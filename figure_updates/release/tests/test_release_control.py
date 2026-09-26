#!/usr/bin/env python3
from __future__ import annotations

import copy
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path


RELEASE = Path(__file__).resolve().parents[1]
PAPER = RELEASE.parent
sys.path.insert(0, str(RELEASE))

import build_archives as archive  # noqa: E402
import release_control as control  # noqa: E402


class ReleaseControlTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="paper-release-control-")
        self.base = Path(self.temporary.name)
        self.paper = self.base / "Paper"
        self.paper.mkdir(mode=0o750)
        self._write("README.txt", b"software/core\n")
        self._write("Fig2/raw_data/acquisition.bin", b"neuro" * 80)
        self._write("FigS3/exports/m26-001/1_L0_rgb.tif", b"histology" * 80)
        self._write("FigS5/channel_tiffs/supplement.bin", b"supplement" * 20)
        self.contract = PAPER / "scripts/utilities/01_validate_bundle.py"
        self.archives = archive.build_archives(
            self.paper,
            self.base / "archive-set",
            contract_script=self.contract,
            payload_limit=1000,
            max_archive_bytes=5000,
            compression_level=6,
        )
        self.plan_path = RELEASE / "release_plan.json"
        self.audit_path = (
            archive.local_input_path(PAPER, "manuscript_sources/PAPER_AUTHOR_PUBMED_AUDIT_20260830.json")
        )
        self.plan = control.load_json(self.plan_path, "release plan")
        self.audit = control.load_json(self.audit_path, "author audit")
        self.summary = control.verify_archive_set(
            self.archives, contract_script=self.contract
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _write(self, relative: str, payload: bytes) -> Path:
        target = self.paper / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        return target

    def approvals(self) -> dict:
        return {
            "schema": control.APPROVALS_SCHEMA,
            "release_date": "2026-08-31",
            "article_doi": None,
            "licenses": {
                "software_core": "other-open",
                "main_figure_data": "cc-by-4.0",
                "supplementary_figure_data": "cc-by-4.0",
                "figs3_wsi": "cc-by-4.0",
            },
            "rights_confirmations": {
                key: True
                for key in self.plan["required_inputs"]["rights_confirmations"]
            },
            "public_visibility_acknowledged": True,
            "clean_room_reproduction_passed": True,
            "credential_scope_confirmations": {
                "zenodo_sandbox": True,
                "zenodo_production": True,
                "github": True,
            },
            "quota_confirmations": {
                "main_figure_data": False,
                "supplementary_figure_data": False,
                "figs3_wsi": False,
            },
            "remote_targets": {
                "github_owner": "release-test-owner",
                "github_repository_name": self.plan["github"]["repository_name"],
                "zenodo_sandbox_record_ids": {
                    "software_core": "101",
                    "main_figure_data": "102",
                    "supplementary_figure_data": "103",
                    "figs3_wsi": "104",
                },
                "zenodo_production_record_ids": {
                    "software_core": "201",
                    "main_figure_data": "202",
                    "supplementary_figure_data": "203",
                    "figs3_wsi": "204",
                },
            },
        }

    @staticmethod
    def environment() -> dict[str, str]:
        return {
            "ZENODO_SANDBOX_TOKEN": "sandbox-secret-never-serialize",
            "ZENODO_TOKEN": "production-secret-never-serialize",
            "GH_TOKEN": "github-secret-never-serialize",
        }

    def test_plan_is_exactly_bound_to_title_authors_split_and_public_visibility(
        self,
    ) -> None:
        result = control.validate_release_plan(self.plan, self.audit, self.audit_path)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["authors"], 18)
        self.assertEqual(result["records"], list(archive.RECORDS))
        self.assertEqual(result["visibility"], control.REQUIRED_VISIBILITY)
        self.assertEqual(
            [row["zenodo_name"] for row in self.plan["manuscript"]["authors"]],
            list(control.EXPECTED_ZENODO_NAMES),
        )

        for mutation in ("title", "author", "visibility", "record"):
            changed = copy.deepcopy(self.plan)
            if mutation == "title":
                changed["manuscript"]["title"] += " changed"
            elif mutation == "author":
                changed["manuscript"]["authors"][0]["display"] = "Wrong Author"
            elif mutation == "visibility":
                changed["visibility"]["github"] = "private"
            else:
                changed["records"].pop("supplementary_figure_data")
            with self.subTest(mutation=mutation):
                with self.assertRaises(control.ReleaseControlError):
                    control.validate_release_plan(changed, self.audit, self.audit_path)

    def test_plan_rejects_missing_malformed_or_reordered_zenodo_names(self) -> None:
        for mutation in ("missing", "malformed", "reordered"):
            changed = copy.deepcopy(self.plan)
            authors = changed["manuscript"]["authors"]
            if mutation == "missing":
                del authors[0]["zenodo_name"]
            elif mutation == "malformed":
                authors[0]["zenodo_name"] = "Nancy Segura"
            else:
                authors[0]["zenodo_name"], authors[1]["zenodo_name"] = (
                    authors[1]["zenodo_name"],
                    authors[0]["zenodo_name"],
                )
            with self.subTest(mutation=mutation):
                with self.assertRaises(control.ReleaseControlError):
                    control.validate_release_plan(changed, self.audit, self.audit_path)

    def test_public_plan_mirror_is_exact(self) -> None:
        public_plan = control.load_json(
            RELEASE / "github_public" / "release_plan.json",
            "public release-plan mirror",
        )
        self.assertEqual(public_plan, self.plan)

    def test_archive_verifier_checks_closed_world_shards_and_members(self) -> None:
        self.assertEqual(self.summary["status"], "PASS")
        self.assertTrue(self.summary["canonical_contract_bound"])
        self.assertEqual(self.summary["tree_files"], 4)
        self.assertEqual(self.summary["verified_zip_members"], 4)
        self.assertEqual(set(self.summary["records"]), set(archive.RECORDS))

        extra = self.archives / "unexpected.txt"
        extra.write_text("not allowed", encoding="utf-8")
        with self.assertRaises(control.ReleaseControlError):
            control.verify_archive_set(self.archives, contract_script=self.contract)

    def test_archive_verifier_rejects_shard_checksum_tampering(self) -> None:
        shard = next(self.archives.glob("*/*.zip"))
        with shard.open("ab") as handle:
            handle.write(b"tampered")
        with self.assertRaisesRegex(
            control.ReleaseControlError, "size differs|SHA-256 differs"
        ):
            control.verify_archive_set(self.archives, contract_script=self.contract)

    def test_required_input_template_is_blocked_and_names_every_gate(self) -> None:
        template = control.load_json(
            RELEASE / "required_inputs.template.json", "approval template"
        )
        result = control.evaluate_readiness(
            self.plan,
            template,
            self.summary,
            "zenodo-production-draft",
            environ={},
        )
        self.assertEqual(result["status"], "BLOCKED")
        joined = "\n".join(result["blockers"])
        for phrase in (
            "license not approved: software_core",
            "rights confirmation is false",
            "public visibility is not acknowledged",
            "clean-room reproduction is not confirmed",
            "release_date is missing",
            "credential scopes not confirmed",
            "credential is absent",
        ):
            self.assertIn(phrase, joined)
        self.assertFalse(result["credential_values_embedded"])
        self.assertFalse(result["remote_execution_enabled"])

    def test_all_operations_can_become_ready_but_never_execute_remotely(self) -> None:
        approved = self.approvals()
        for operation in control.OPERATION_POLICY:
            with self.subTest(operation=operation):
                result = control.evaluate_readiness(
                    self.plan,
                    approved,
                    self.summary,
                    operation,
                    environ=self.environment(),
                )
                self.assertEqual(result["status"], "PASS", result["blockers"])
                self.assertTrue(result["credential_present"])
                self.assertFalse(result["remote_execution_enabled"])
                with self.assertRaisesRegex(
                    control.ReleaseControlError, "intentionally disabled"
                ):
                    control.remote_action(operation)

    def test_quota_is_deferred_for_draft_but_required_for_upload(self) -> None:
        approved = self.approvals()
        large = copy.deepcopy(self.summary)
        large["records"]["main_figure_data"]["compressed_archive_bytes"] = (
            archive.ZENODO_DEFAULT_RECORD_BYTES + 1
        )
        draft = control.evaluate_readiness(
            self.plan,
            approved,
            large,
            "zenodo-production-draft",
            environ=self.environment(),
        )
        self.assertEqual(draft["status"], "PASS", draft["blockers"])
        self.assertNotIn(
            "Zenodo >50 GB quota not confirmed: main_figure_data",
            draft["blockers"],
        )

        request = control.prepare_request_document(
            self.plan,
            approved,
            large,
            "zenodo-production-draft",
            environ=self.environment(),
        )
        self.assertFalse(
            request["approvals"]["quota_confirmations"]["main_figure_data"]
        )

        blocked = control.evaluate_readiness(
            self.plan,
            approved,
            large,
            "zenodo-production-upload",
            environ=self.environment(),
        )
        self.assertIn(
            "Zenodo >50 GB quota not confirmed: main_figure_data",
            blocked["blockers"],
        )
        approved["quota_confirmations"]["main_figure_data"] = True
        wrong_license = copy.deepcopy(approved)
        wrong_license["licenses"]["software_core"] = "mit"
        blocked = control.evaluate_readiness(
            self.plan,
            wrong_license,
            large,
            "zenodo-production-draft",
            environ=self.environment(),
        )
        self.assertIn(
            "license differs from approved plan: software_core", blocked["blockers"]
        )
        approved["remote_targets"]["zenodo_production_record_ids"][
            "software_core"
        ] = None
        blocked = control.evaluate_readiness(
            self.plan,
            approved,
            large,
            "zenodo-production-upload",
            environ=self.environment(),
        )
        self.assertIn("remote record ID is missing: software_core", blocked["blockers"])

    def test_prepared_request_is_public_four_record_and_credential_free(self) -> None:
        environment = self.environment()
        request = control.prepare_request_document(
            self.plan,
            self.approvals(),
            self.summary,
            "zenodo-production-draft",
            environ=environment,
        )
        self.assertEqual(request["schema"], control.REQUEST_SCHEMA)
        self.assertEqual(request["operation"], "zenodo-production-draft")
        self.assertTrue(request["local_only"])
        self.assertFalse(request["remote_execution_enabled"])
        self.assertFalse(request["credentials_embedded"])
        self.assertEqual(set(request["zenodo_records"]), set(archive.RECORDS))
        displays = [
            row["manuscript_display"] for row in request["manuscript"]["authors"]
        ]
        self.assertEqual(displays, [row["display"] for row in self.audit["authors"]])
        zenodo_names = [row["name"] for row in request["manuscript"]["authors"]]
        self.assertEqual(zenodo_names, list(control.EXPECTED_ZENODO_NAMES))
        self.assertEqual(
            [row["service_name"] for row in request["manuscript"]["authors"]],
            [row["service_name"] for row in self.plan["manuscript"]["authors"]],
        )
        for row in request["zenodo_records"].values():
            self.assertEqual(row["access"], control.REQUIRED_VISIBILITY["zenodo"])
            self.assertGreater(len(row["archive_shards"]), 0)
            self.assertEqual(
                [creator["name"] for creator in row["metadata"]["creators"]],
                list(control.EXPECTED_ZENODO_NAMES),
            )
        self.assertEqual(
            request["zenodo_records"]["software_core"]["ancillary_files"],
            [
                {
                    "name": "manifest.json",
                    "size_bytes": self.summary["manifest_size_bytes"],
                    "sha256": self.summary["manifest_sha256"],
                }
            ],
        )
        for record in archive.RECORDS[1:]:
            self.assertEqual(request["zenodo_records"][record]["ancillary_files"], [])
        serialized = control.canonical_json(request)
        for secret in environment.values():
            self.assertNotIn(secret.encode("utf-8"), serialized)

        target = self.base / "request.json"
        control.write_request(target, request)
        self.assertEqual(stat.S_IMODE(target.stat().st_mode), 0o600)
        self.assertEqual(json.loads(target.read_text(encoding="utf-8")), request)
        with self.assertRaises(control.ReleaseControlError):
            control.write_request(target, request)


if __name__ == "__main__":
    unittest.main()
