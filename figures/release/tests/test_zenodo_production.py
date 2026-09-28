#!/usr/bin/env python3
from __future__ import annotations

import copy
import hashlib
import json
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock


RELEASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RELEASE))

import build_archives as archives  # noqa: E402
import zenodo_production as production  # noqa: E402


def digest(payload: bytes, name: str) -> str:
    if name == "sha256":
        return hashlib.sha256(payload).hexdigest()
    return hashlib.md5(payload, usedforsecurity=False).hexdigest()


class FakeZenodo:
    def __init__(self, payloads: dict[str, bytes]) -> None:
        self.payloads = payloads
        self.depositions: dict[int, dict] = {}
        self.creates = 0
        self.updates = 0
        self.uploads = 0
        self.publishes = 0
        self.calls = 0

    def find_drafts(self, release_id: str) -> list[dict]:
        self.calls += 1
        return [
            value
            for value in self.depositions.values()
            if value["submitted"] is False
            and release_id in value["metadata"].get("keywords", [])
        ]

    def create_deposition(self, metadata: dict) -> dict:
        self.calls += 1
        self.creates += 1
        identifier = 100 + self.creates
        stored_metadata = copy.deepcopy(metadata)
        stored_metadata["prereserve_doi"] = {
            "doi": f"10.5281/zenodo.{identifier}",
            "recid": identifier,
        }
        value = {
            "id": identifier,
            "conceptrecid": str(identifier - 1),
            "submitted": False,
            "state": "inprogress",
            "metadata": stored_metadata,
            "files": [],
            "links": {"bucket": f"https://zenodo.org/api/files/bucket-{identifier}"},
        }
        self.depositions[identifier] = value
        return copy.deepcopy(value)

    def get_deposition(self, deposition_id: int) -> dict:
        self.calls += 1
        return copy.deepcopy(self.depositions[int(deposition_id)])

    def update_deposition(self, deposition_id: int, metadata: dict) -> dict:
        self.calls += 1
        self.updates += 1
        stored = copy.deepcopy(metadata)
        stored["prereserve_doi"] = copy.deepcopy(
            self.depositions[int(deposition_id)]["metadata"]["prereserve_doi"]
        )
        self.depositions[int(deposition_id)]["metadata"] = stored
        return self.get_deposition(deposition_id)

    def upload_file(
        self, bucket_url: str, filename: str, path: Path, size_bytes: int
    ) -> dict:
        self.calls += 1
        self.uploads += 1
        identifier = int(bucket_url.rsplit("-", 1)[1])
        payload = path.read_bytes()
        assert len(payload) == size_bytes
        raw = {
            "key": filename,
            "size": len(payload),
            "checksum": "md5:" + digest(payload, "md5"),
            "links": {
                "self": f"https://zenodo.org/api/files/bucket-{identifier}/{filename}"
            },
        }
        self.depositions[identifier]["files"].append(raw)
        return copy.deepcopy(raw)

    def download_file(
        self,
        url: str,
        target: Path,
        *,
        expected_size: int,
        expected_sha256: str,
    ) -> None:
        self.calls += 1
        filename = url.rsplit("/", 1)[1]
        payload = self.payloads[filename]
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)

    def publish(self, deposition_id: int) -> dict:
        self.calls += 1
        self.publishes += 1
        row = self.depositions[int(deposition_id)]
        row["submitted"] = True
        row["state"] = "done"
        row["record_id"] = int(deposition_id)
        row["record_url"] = f"https://zenodo.org/records/{deposition_id}"
        return copy.deepcopy(row)

    def get_public_record(self, record_id: int) -> dict:
        self.calls += 1
        return self.get_deposition(record_id)


class ZenodoProductionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="zenodo-production-test-")
        self.base = Path(self.temporary.name)
        self.archives = self.base / "archives"
        self.archives.mkdir()
        self.payloads = {
            "software.part001.zip": b"software archive",
            "manifest.json": b'{"manifest":"remote"}\n',
            "main.part001.zip": b"main archive",
            "supp.part001.zip": b"supp archive",
            "wsi.part001.zip": b"wsi archive",
        }
        specs = {
            "software_core": [
                ("software.part001.zip", "software_core/software.part001.zip"),
                ("manifest.json", "manifest.json"),
            ],
            "main_figure_data": [
                ("main.part001.zip", "main_figure_data/main.part001.zip")
            ],
            "supplementary_figure_data": [
                ("supp.part001.zip", "supplementary_figure_data/supp.part001.zip")
            ],
            "figs3_wsi": [("wsi.part001.zip", "figs3_wsi/wsi.part001.zip")],
        }
        records: dict[str, dict] = {}
        for index, record in enumerate(production.RECORDS, 1):
            files = []
            for filename, relative in specs[record]:
                payload = self.payloads[filename]
                path = self.archives / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(payload)
                files.append(
                    {
                        "name": filename,
                        "relative_local_path": relative,
                        "size_bytes": len(payload),
                        "sha256": digest(payload, "sha256"),
                        "md5": None,
                        "uploaded": False,
                        "remote_verified": False,
                    }
                )
            records[record] = {
                "title": f"Reviewed title: {record}",
                "description": f"Reviewed description for {record}",
                "upload_type": "software" if record == "software_core" else "dataset",
                "publication_date": "2026-09-02",
                "creators": [
                    {"name": name} for name in production.control.EXPECTED_ZENODO_NAMES
                ],
                "license": production.EXPECTED_LICENSES[record],
                "article_doi": None,
                "files": files,
                "deposition_id": None,
                "reserved_doi": None,
                "metadata_synced": False,
                "remote_verified": False,
                "published": False,
                "public_url": None,
                "quota": {
                    "required_record_bytes": sum(row["size_bytes"] for row in files),
                    "confirmation_required": False,
                    "confirmed_for_deposition_id": None,
                    "allocated_record_bytes": archives.ZENODO_DEFAULT_RECORD_BYTES,
                },
            }
        manifest_payload = self.payloads["manifest.json"]
        self.state = {
            "schema": production.STATE_SCHEMA,
            "environment": "production",
            "api_base": production.PRODUCTION_API,
            "release_id": "apotome-" + "a" * 24,
            "request_sha256": "b" * 64,
            "reviewed_gates": {
                "request_operation": "zenodo-production-draft",
                "rights_confirmations": {
                    key: True for key in production.EXPECTED_RIGHTS_CONFIRMATIONS
                },
                "public_visibility_acknowledged": True,
                "local_clean_room_reproduction_passed": True,
            },
            "archive_root": str(self.archives),
            "manifest": {
                "sha256": digest(manifest_payload, "sha256"),
                "size_bytes": len(manifest_payload),
                "tree_sha256": "c" * 64,
                "tree_files": 4,
            },
            "github_url": "https://github.com/cfarkas/allulose-sucrose-hypothalamus",
            "related_identifier_policy": {
                "mutual_reserved_dois": "isSupplementTo",
                "github_repository": "isSupplementTo",
                "article_doi": "isSupplementTo",
            },
            "records": records,
            "remote_download": {
                "status": "NOT_STARTED",
                "root": None,
                "completed_files": [],
                "reconstructed_root": None,
            },
            "reproduction": {
                "status": "NOT_STARTED",
                "log_path": None,
                "log_sha256": None,
                "command": None,
            },
            "created_at": "2026-09-02T00:00:00+00:00",
            "updated_at": "2026-09-02T00:00:00+00:00",
        }
        self.client = FakeZenodo(self.payloads)
        self.saved: list[dict] = []
        self.persist = lambda state: self.saved.append(copy.deepcopy(state))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def draft_upload(self) -> None:
        production.create_and_link_drafts(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )
        production.upload_all(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )

    def make_large_unconfirmed(self, record: str, *, extra_bytes: int = 1) -> dict:
        row = self.state["records"][record]
        required = archives.ZENODO_DEFAULT_RECORD_BYTES + extra_bytes
        remaining = sum(item["size_bytes"] for item in row["files"][1:])
        row["files"][0]["size_bytes"] = required - remaining
        row["quota"] = {
            "required_record_bytes": required,
            "confirmation_required": True,
            "confirmed_for_deposition_id": None,
            "allocated_record_bytes": None,
        }
        return row

    def seed_explicit_adoption(
        self, record: str, *, wrong_title: bool = False, submitted: bool = False
    ) -> int:
        metadata = production._create_metadata(self.state, record)
        metadata["description"] = "Description from the superseded three-record request"
        metadata["keywords"][-1] = "apotome-" + "d" * 24
        metadata["notes"] = (
            "Three-record reproducibility release. Executor marker: apotome-" + "d" * 24
        )
        metadata["related_identifiers"] = [
            {
                "identifier": "10.5281/zenodo.999999",
                "relation": "isSupplementTo",
                "resource_type": "dataset",
            }
        ]
        if wrong_title:
            metadata["title"] += " wrong"
        deposition = self.client.create_deposition(metadata)
        deposition_id = int(deposition["id"])
        self.client.depositions[deposition_id]["submitted"] = submitted
        self.state["records"][record]["deposition_id"] = deposition_id
        return deposition_id

    def bind_segmented_transport(self) -> dict:
        records = {}
        for record in production.RECORDS:
            rows = []
            for item in self.state["records"][record]["files"]:
                if item["name"] == "manifest.json":
                    continue
                payload = self.payloads[item["name"]]
                segment_bytes = max(1, len(payload) // 2)
                count = (len(payload) + segment_bytes - 1) // segment_bytes
                segments = []
                offset = 0
                for number in range(1, count + 1):
                    chunk = payload[offset : offset + segment_bytes]
                    name = f"{item['name']}.chunk{number:03d}-of-{count:03d}"
                    self.payloads[name] = chunk
                    segments.append(
                        {
                            "name": name,
                            "number": number,
                            "offset_bytes": offset,
                            "size_bytes": len(chunk),
                            "sha256": digest(chunk, "sha256"),
                            "md5": digest(chunk, "md5"),
                        }
                    )
                    offset += len(chunk)
                rows.append(
                    {
                        "name": item["name"],
                        "relative_local_path": item["relative_local_path"],
                        "size_bytes": item["size_bytes"],
                        "sha256": item["sha256"],
                        "md5": digest(payload, "md5"),
                        "segment_bytes": segment_bytes,
                        "segments": segments,
                    }
                )
            records[record] = {"archives": rows}
        document = {
            "schema": production.SEGMENTED_TRANSPORT_SCHEMA,
            "release_id": self.state["release_id"],
            "canonical_manifest": copy.deepcopy(self.state["manifest"]),
            "records": records,
            "reassembly": production.TRANSPORT_REASSEMBLY_INSTRUCTION,
        }
        path = self.base / production.TRANSPORT_MANIFEST_FILENAME
        path.write_bytes(production.canonical_json(document))
        self.payloads[production.TRANSPORT_MANIFEST_FILENAME] = path.read_bytes()
        production.bind_segmented_transport(self.state, path)
        return document

    def seed_segmented_remote(self, *, verify: bool = True) -> None:
        production.create_and_link_drafts(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )
        for record in production.RECORDS:
            deposition_id = self.state["records"][record]["deposition_id"]
            for item in production.expected_remote_items(self.state, record):
                payload = self.payloads[item["name"]]
                self.client.depositions[deposition_id]["files"].append(
                    {
                        "key": item["name"],
                        "size": len(payload),
                        "checksum": "md5:" + digest(payload, "md5"),
                        "links": {
                            "self": (
                                "https://zenodo.org/api/files/"
                                f"bucket-{deposition_id}/{item['name']}"
                            )
                        },
                    }
                )
        if verify:
            production.verify_remote(self.state, self.client, self.persist)

    def test_adopt_draft_parser_requires_unique_known_records_and_positive_ids(
        self,
    ) -> None:
        arguments = production.parser().parse_args(
            [
                "--state",
                str(self.base / "state.json"),
                "init",
                "--request",
                str(self.base / "request.json"),
                "--archive-root",
                str(self.archives),
                "--adopt-draft",
                "software_core=101",
                "--adopt-draft",
                "main_figure_data=102",
            ]
        )
        self.assertEqual(
            production.parse_adopted_drafts(arguments.adopt_draft),
            {"software_core": 101, "main_figure_data": 102},
        )
        for values in (
            ["unknown=101"],
            ["software_core=0"],
            ["software_core=-1"],
            ["software_core=101", "software_core=102"],
            ["software_core=101", "main_figure_data=101"],
        ):
            with self.subTest(values=values):
                with self.assertRaises(production.ZenodoReleaseError):
                    production.parse_adopted_drafts(values)

    def test_explicit_adoption_preserves_empty_draft_but_not_quota_authority(
        self,
    ) -> None:
        row = self.make_large_unconfirmed("main_figure_data")
        deposition_id = self.seed_explicit_adoption("main_figure_data")
        creates_before = self.client.creates
        production.create_and_link_drafts(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )
        self.assertEqual(
            self.client.creates - creates_before, len(production.RECORDS) - 1
        )
        self.assertEqual(row["deposition_id"], deposition_id)
        self.assertEqual(row["reserved_doi"], f"10.5281/zenodo.{deposition_id}")
        self.assertEqual(row["concept_record_id"], str(deposition_id - 1))
        self.assertEqual(
            row["bucket_url"], f"https://zenodo.org/api/files/bucket-{deposition_id}"
        )
        self.assertTrue(row["metadata_synced"])
        self.assertIsNone(row["quota"]["confirmed_for_deposition_id"])
        self.assertIsNone(row["quota"]["allocated_record_bytes"])

    def test_explicit_adoption_preflights_wrong_title_and_submitted_state(self) -> None:
        base_state = copy.deepcopy(self.state)
        for mutation in ("title", "submitted"):
            with self.subTest(mutation=mutation):
                self.state = copy.deepcopy(base_state)
                self.client = FakeZenodo(self.payloads)
                self.saved = []
                self.persist = lambda state: self.saved.append(copy.deepcopy(state))
                self.seed_explicit_adoption(
                    "main_figure_data",
                    wrong_title=mutation == "title",
                    submitted=mutation == "submitted",
                )
                creates_before = self.client.creates
                updates_before = self.client.updates
                with self.assertRaisesRegex(
                    production.ZenodoReleaseError,
                    "identity differs|not an unpublished draft",
                ):
                    production.create_and_link_drafts(
                        self.state,
                        self.client,
                        self.persist,
                        execute=True,
                        production=True,
                    )
                self.assertEqual(self.client.creates, creates_before)
                self.assertEqual(self.client.updates, updates_before)

    def test_init_binds_reviewed_request_to_dynamic_manifest(self) -> None:
        rights_keys = list(production.EXPECTED_RIGHTS_CONFIRMATIONS)
        plan = {
            "required_inputs": {
                "rights_confirmations": rights_keys,
                "quota_confirmations": list(production.RECORDS[1:]),
            },
            "manuscript": {
                "authors": [
                    {"zenodo_name": name}
                    for name in production.control.EXPECTED_ZENODO_NAMES
                ]
            },
            "records": {
                record: {
                    "title": self.state["records"][record]["title"],
                    "description": self.state["records"][record]["description"],
                    "resource_type": self.state["records"][record]["upload_type"],
                }
                for record in production.RECORDS
            },
        }
        shard_names = {
            "software_core": "software.part001.zip",
            "main_figure_data": "main.part001.zip",
            "supplementary_figure_data": "supp.part001.zip",
            "figs3_wsi": "wsi.part001.zip",
        }
        summary_records = {}
        request_records = {}
        for record in production.RECORDS:
            filename = shard_names[record]
            payload = self.payloads[filename]
            shard = {
                "name": filename,
                "size_bytes": len(payload),
                "sha256": digest(payload, "sha256"),
                "member_count": 1,
            }
            summary_records[record] = {
                "compressed_archive_bytes": len(payload),
                "shard_count": 1,
                "shards": [shard],
            }
            row = self.state["records"][record]
            request_records[record] = {
                "metadata": {
                    "title": row["title"],
                    "description": row["description"],
                    "resource_type": row["upload_type"],
                    "publication_date": row["publication_date"],
                    "license": row["license"],
                    "creators": [
                        {"name": name}
                        for name in production.control.EXPECTED_ZENODO_NAMES
                    ],
                },
                "access": production.control.REQUIRED_VISIBILITY["zenodo"],
                "archive_shards": [shard],
                "ancillary_files": (
                    [
                        {
                            "name": "manifest.json",
                            "size_bytes": self.state["manifest"]["size_bytes"],
                            "sha256": self.state["manifest"]["sha256"],
                        }
                    ]
                    if record == "software_core"
                    else []
                ),
                "article_doi": None,
            }
        summary = {
            "status": "PASS",
            "manifest_sha256": self.state["manifest"]["sha256"],
            "manifest_size_bytes": self.state["manifest"]["size_bytes"],
            "tree_sha256": self.state["manifest"]["tree_sha256"],
            "tree_files": self.state["manifest"]["tree_files"],
            "records": summary_records,
        }
        request = {
            "schema": production.REQUEST_SCHEMA,
            "operation": "zenodo-production-draft",
            "status": "READY_FOR_SEPARATELY_REVIEWED_EXECUTOR",
            "local_only": True,
            "remote_execution_enabled": False,
            "credentials_embedded": False,
            "visibility": production.control.REQUIRED_VISIBILITY,
            "approvals": {
                "public_visibility_acknowledged": True,
                "clean_room_reproduction_passed": True,
                "rights_confirmations": {key: True for key in rights_keys},
                "quota_confirmations": {
                    "main_figure_data": False,
                    "supplementary_figure_data": False,
                    "figs3_wsi": False,
                },
            },
            "zenodo_records": request_records,
            "github": {
                "owner": "cfarkas",
                "repository_name": "allulose-sucrose-hypothalamus",
            },
            "archive_manifest": {
                "sha256": summary["manifest_sha256"],
                "tree_sha256": summary["tree_sha256"],
                "tree_files": summary["tree_files"],
            },
        }
        request_path = self.base / "request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        with mock.patch.object(
            production.control, "load_validated_plan", return_value=(plan, {})
        ), mock.patch.object(
            production.control, "verify_archive_set", return_value=summary
        ), mock.patch.object(
            production.archives, "ZENODO_DEFAULT_RECORD_BYTES", 1
        ):
            observed = production.initialize_state(
                request_path,
                self.archives,
                plan_path=self.base / "plan.json",
                audit_path=self.base / "audit.json",
                contract_script=self.base / "contract.py",
                adopt_drafts={"figs3_wsi": 104},
            )
        self.assertEqual(observed["records"]["figs3_wsi"]["deposition_id"], 104)
        self.assertIsNone(
            observed["records"]["figs3_wsi"]["quota"]["confirmed_for_deposition_id"]
        )
        observed["records"]["figs3_wsi"]["deposition_id"] = None
        self.assertEqual(
            observed["manifest"]["sha256"], self.state["manifest"]["sha256"]
        )
        self.assertEqual(observed["records"]["software_core"]["license"], "other-open")
        self.assertEqual(
            observed["records"]["main_figure_data"]["license"], "cc-by-4.0"
        )
        self.assertEqual(
            [row["name"] for row in observed["records"]["software_core"]["creators"]],
            list(production.control.EXPECTED_ZENODO_NAMES),
        )
        self.assertTrue(
            observed["records"]["main_figure_data"]["quota"]["confirmation_required"]
        )
        self.assertIsNone(
            observed["records"]["main_figure_data"]["quota"][
                "confirmed_for_deposition_id"
            ]
        )
        self.assertIsNone(
            observed["records"]["main_figure_data"]["quota"]["allocated_record_bytes"]
        )

        client = FakeZenodo(self.payloads)
        saved: list[dict] = []
        persist = lambda state: saved.append(copy.deepcopy(state))
        with mock.patch.object(production.archives, "ZENODO_DEFAULT_RECORD_BYTES", 1):
            production.create_and_link_drafts(
                observed, client, persist, execute=True, production=True
            )
            self.assertEqual(client.creates, len(production.RECORDS))
            with self.assertRaisesRegex(production.ZenodoReleaseError, "quota"):
                production.upload_all(
                    observed, client, persist, execute=True, production=True
                )
        self.assertEqual(client.uploads, 0)

    def test_request_rejects_missing_malformed_or_reordered_zenodo_names(self) -> None:
        rights_keys = list(production.EXPECTED_RIGHTS_CONFIRMATIONS)
        base_plan = {
            "required_inputs": {
                "rights_confirmations": rights_keys,
                "quota_confirmations": list(production.RECORDS[1:]),
            },
            "manuscript": {
                "authors": [
                    {"zenodo_name": name}
                    for name in production.control.EXPECTED_ZENODO_NAMES
                ]
            },
            "records": {
                record: {
                    "title": self.state["records"][record]["title"],
                    "description": self.state["records"][record]["description"],
                    "resource_type": self.state["records"][record]["upload_type"],
                }
                for record in production.RECORDS
            },
        }
        request_records = {
            record: {
                "metadata": {
                    "title": self.state["records"][record]["title"],
                    "description": self.state["records"][record]["description"],
                    "resource_type": self.state["records"][record]["upload_type"],
                    "license": self.state["records"][record]["license"],
                    "creators": copy.deepcopy(
                        self.state["records"][record]["creators"]
                    ),
                },
                "access": production.control.REQUIRED_VISIBILITY["zenodo"],
            }
            for record in production.RECORDS
        }
        base_request = {
            "schema": production.REQUEST_SCHEMA,
            "operation": "zenodo-production-draft",
            "status": "READY_FOR_SEPARATELY_REVIEWED_EXECUTOR",
            "local_only": True,
            "remote_execution_enabled": False,
            "credentials_embedded": False,
            "visibility": production.control.REQUIRED_VISIBILITY,
            "approvals": {
                "public_visibility_acknowledged": True,
                "clean_room_reproduction_passed": True,
                "rights_confirmations": {key: True for key in rights_keys},
                "quota_confirmations": {},
            },
            "zenodo_records": request_records,
        }
        for mutation in ("missing", "malformed", "reordered"):
            plan = copy.deepcopy(base_plan)
            request = copy.deepcopy(base_request)
            if mutation == "missing":
                del plan["manuscript"]["authors"][0]["zenodo_name"]
            elif mutation == "malformed":
                plan["manuscript"]["authors"][0]["zenodo_name"] = "Nancy Segura"
            else:
                creators = request["zenodo_records"]["software_core"]["metadata"][
                    "creators"
                ]
                creators[0], creators[1] = creators[1], creators[0]
            with self.subTest(mutation=mutation):
                with self.assertRaises(production.ZenodoReleaseError):
                    production._validate_request(request, plan)

    def test_every_remote_mutation_requires_explicit_gate(self) -> None:
        with self.assertRaises(production.ZenodoReleaseError):
            production.create_and_link_drafts(
                self.state,
                self.client,
                self.persist,
                execute=False,
                production=True,
            )
        with self.assertRaises(production.ZenodoReleaseError):
            production.upload_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=False,
            )
        with self.assertRaises(production.ZenodoReleaseError):
            production.publish_all(
                self.state,
                self.client,
                self.persist,
                execute=False,
                production=True,
                confirmation=production.PUBLISH_CONFIRMATION,
                expected_state_sha256=production.state_sha256(self.state),
            )
        self.assertEqual(self.client.calls, 0)

    def test_state_requires_exact_complete_rights_gate_set(self) -> None:
        production.validate_state(self.state)

        missing = copy.deepcopy(self.state)
        del missing["reviewed_gates"]["rights_confirmations"][
            production.EXPECTED_RIGHTS_CONFIRMATIONS[-1]
        ]
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "rights/public/local-clean-room"
        ):
            production.validate_state(missing)

        extra = copy.deepcopy(self.state)
        extra["reviewed_gates"]["rights_confirmations"]["unreviewed_extra"] = True
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "rights/public/local-clean-room"
        ):
            production.validate_state(extra)

        reordered = copy.deepcopy(self.state)
        creators = reordered["records"]["software_core"]["creators"]
        creators[0], creators[1] = creators[1], creators[0]
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "Creator names/order"
        ):
            production.validate_state(reordered)

        malformed = copy.deepcopy(self.state)
        malformed["records"]["software_core"]["creators"][0]["name"] = "Nancy Segura"
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "Creator names/order"
        ):
            production.validate_state(malformed)

    def test_drafts_reserve_then_crosslink_all_dois_and_mixed_licenses(self) -> None:
        production.create_and_link_drafts(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )
        self.assertEqual(self.client.creates, len(production.RECORDS))
        self.assertEqual(self.client.updates, len(production.RECORDS))
        dois = {
            self.state["records"][record]["reserved_doi"]
            for record in production.RECORDS
        }
        self.assertEqual(len(dois), len(production.RECORDS))
        for record in production.RECORDS:
            metadata = production.effective_metadata(self.state, record)
            linked_dois = {
                item["identifier"]
                for item in metadata["related_identifiers"]
                if item["identifier"].startswith("10.")
            }
            self.assertEqual(
                linked_dois, dois - {self.state["records"][record]["reserved_doi"]}
            )
            self.assertEqual(metadata["license"], production.EXPECTED_LICENSES[record])
            self.assertEqual(metadata["access_right"], "open")
            self.assertTrue(self.state["records"][record]["metadata_synced"])

    def test_drafts_resynchronize_reviewed_description_changes(self) -> None:
        production.create_and_link_drafts(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )
        creates = self.client.creates
        updates = self.client.updates
        for record in production.RECORDS:
            self.state["records"][record][
                "description"
            ] += " Segments concatenate to the manifest-bound canonical archive."

        production.create_and_link_drafts(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )

        self.assertEqual(self.client.creates, creates)
        self.assertEqual(self.client.updates, updates + len(production.RECORDS))
        for record in production.RECORDS:
            deposition_id = self.state["records"][record]["deposition_id"]
            production.verify_metadata(
                self.client.depositions[deposition_id],
                production.effective_metadata(self.state, record),
                require_related=True,
            )

    def test_metadata_accepts_only_zenodo_ampersand_serialization(self) -> None:
        expected = production._create_metadata(self.state, "figs3_wsi")
        expected["description"] = "Reviewed mouse H&E whole-slide exports"
        deposition = {"metadata": copy.deepcopy(expected)}
        production.verify_metadata(deposition, expected, require_related=False)
        deposition["metadata"][
            "description"
        ] = "Reviewed mouse H&amp;E whole-slide exports"
        production.verify_metadata(deposition, expected, require_related=False)

        for changed in (
            "Reviewed mouse H and E whole-slide exports",
            "Reviewed mouse H&amp;amp;E whole-slide exports",
            "Reviewed mouse H&#38;E whole-slide exports",
            "Reviewed mouse H&amp;E altered exports",
            "Reviewed mouse H&amp;E whole-slide exports &lt;script&gt;",
        ):
            deposition["metadata"]["description"] = changed
            with self.assertRaisesRegex(
                production.ZenodoReleaseError, "metadata differs"
            ):
                production.verify_metadata(deposition, expected, require_related=False)

        deposition["metadata"] = copy.deepcopy(expected)
        deposition["metadata"]["title"] = expected["title"].replace("&", "&amp;")
        if deposition["metadata"]["title"] == expected["title"]:
            deposition["metadata"]["title"] += " &amp; encoded"
        with self.assertRaisesRegex(production.ZenodoReleaseError, "metadata differs"):
            production.verify_metadata(deposition, expected, require_related=False)

        expected["description"] = "Reviewed literal H&amp;E label"
        deposition["metadata"] = copy.deepcopy(expected)
        production.verify_metadata(deposition, expected, require_related=False)
        deposition["metadata"]["description"] = "Reviewed literal H&amp;amp;E label"
        with self.assertRaisesRegex(production.ZenodoReleaseError, "metadata differs"):
            production.verify_metadata(deposition, expected, require_related=False)

    def test_quota_confirmation_is_bound_to_exact_draft_and_bytes(self) -> None:
        row = self.make_large_unconfirmed("main_figure_data", extra_bytes=123)
        row["deposition_id"] = 102
        quota = row["quota"]
        with self.assertRaises(production.ZenodoReleaseError):
            production.confirm_quota(
                self.state,
                "main_figure_data",
                deposition_id=999,
                allocated_record_bytes=quota["required_record_bytes"],
                confirmed=True,
            )
        production.confirm_quota(
            self.state,
            "main_figure_data",
            deposition_id=102,
            allocated_record_bytes=quota["required_record_bytes"],
            confirmed=True,
        )
        production._quota_ready(
            "main_figure_data", self.state["records"]["main_figure_data"]
        )
        row["deposition_id"] = 103
        with self.assertRaisesRegex(production.ZenodoReleaseError, "quota"):
            production.validate_state(self.state)
        with self.assertRaisesRegex(production.ZenodoReleaseError, "quota"):
            production._quota_ready(
                "main_figure_data", self.state["records"]["main_figure_data"]
            )

    def test_quota_state_rejects_lowered_derived_required_bytes(self) -> None:
        row = self.state["records"]["main_figure_data"]
        row["quota"]["required_record_bytes"] -= 1
        with self.assertRaisesRegex(production.ZenodoReleaseError, "required bytes"):
            production.validate_state(self.state)
        with self.assertRaisesRegex(production.ZenodoReleaseError, "required bytes"):
            production._quota_ready("main_figure_data", row)

    def test_quota_state_rejects_flipped_large_record_threshold(self) -> None:
        row = self.make_large_unconfirmed("main_figure_data")
        row["quota"]["confirmation_required"] = False
        row["quota"]["allocated_record_bytes"] = archives.ZENODO_DEFAULT_RECORD_BYTES
        with self.assertRaisesRegex(production.ZenodoReleaseError, "threshold"):
            production.validate_state(self.state)
        with self.assertRaisesRegex(production.ZenodoReleaseError, "threshold"):
            production._quota_ready("main_figure_data", row)

    def test_quota_state_rejects_inconsistent_small_record_allocation(self) -> None:
        row = self.state["records"]["software_core"]
        row["quota"]["allocated_record_bytes"] -= 1
        with self.assertRaisesRegex(production.ZenodoReleaseError, "Small"):
            production.validate_state(self.state)
        with self.assertRaisesRegex(production.ZenodoReleaseError, "Small"):
            production._quota_ready("software_core", row)

    def test_quota_state_rejects_partial_large_confirmation(self) -> None:
        row = self.make_large_unconfirmed("main_figure_data")
        row["quota"]["allocated_record_bytes"] = row["quota"]["required_record_bytes"]
        with self.assertRaisesRegex(production.ZenodoReleaseError, "current draft"):
            production.validate_state(self.state)
        with self.assertRaisesRegex(production.ZenodoReleaseError, "current draft"):
            production._quota_ready("main_figure_data", row)

    def test_upload_preflights_all_quotas_before_any_network_action(self) -> None:
        production.create_and_link_drafts(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )
        self.make_large_unconfirmed("figs3_wsi")
        calls_before = self.client.calls
        with self.assertRaisesRegex(production.ZenodoReleaseError, "quota"):
            production.upload_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=True,
            )
        self.assertEqual(self.client.calls, calls_before)
        self.assertEqual(self.client.uploads, 0)

    def test_publish_rejects_unconfirmed_draft_bound_quota(self) -> None:
        self.prepare_publish_gate()
        self.make_large_unconfirmed("main_figure_data")
        with self.assertRaisesRegex(production.ZenodoReleaseError, "quota"):
            production.publish_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=True,
                confirmation=production.PUBLISH_CONFIRMATION,
                expected_state_sha256=production.state_sha256(self.state),
            )
        self.assertEqual(self.client.publishes, 0)

    def test_upload_is_file_resumable_and_conflicts_fail_closed(self) -> None:
        self.draft_upload()
        first_uploads = self.client.uploads
        self.assertEqual(first_uploads, 5)
        production.upload_all(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
        )
        self.assertEqual(self.client.uploads, first_uploads)
        first = self.client.depositions[101]["files"][0]
        first["checksum"] = "md5:" + "0" * 32
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "conflicts|size/MD5"
        ):
            production.upload_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=True,
            )

    def test_segmented_transport_binding_and_expected_remote_inventory(self) -> None:
        document = self.bind_segmented_transport()
        production.validate_state(self.state)
        software = production.expected_remote_items(self.state, "software_core")
        self.assertEqual(
            {item["transport_kind"] for item in software},
            {"segment", "canonical_manifest", "transport_manifest"},
        )
        self.assertNotIn("software.part001.zip", {item["name"] for item in software})
        self.assertIn("manifest.json", {item["name"] for item in software})
        self.assertIn(
            production.TRANSPORT_MANIFEST_FILENAME,
            {item["name"] for item in software},
        )

        document["records"]["software_core"]["archives"][0]["segments"][0][
            "offset_bytes"
        ] = 1
        path = Path(self.state["segmented_transport"]["manifest_path"])
        path.write_bytes(production.canonical_json(document))
        sha256, md5, size = production.hashes_file(path)
        self.state["segmented_transport"].update(
            {"sha256": sha256, "md5": md5, "size_bytes": size}
        )
        with self.assertRaisesRegex(production.ZenodoReleaseError, "segment differs"):
            production.validate_state(self.state)

    def test_segmented_transport_storage_override_routes_physical_inventory(
        self,
    ) -> None:
        document = self.bind_segmented_transport()
        logical_record = "figs3_wsi"
        storage_record = "supplementary_figure_data"
        segment = document["records"][logical_record]["archives"][0]["segments"][0]
        name = segment["name"]
        segment["storage_record"] = storage_record
        path = Path(self.state["segmented_transport"]["manifest_path"])
        path.write_bytes(production.canonical_json(document))
        self.payloads[production.TRANSPORT_MANIFEST_FILENAME] = path.read_bytes()
        sha256, md5, size = production.hashes_file(path)
        self.state["segmented_transport"].update(
            {"sha256": sha256, "md5": md5, "size_bytes": size}
        )

        production.validate_state(self.state)
        self.assertNotIn(
            name,
            {
                item["name"]
                for item in production.expected_remote_items(self.state, logical_record)
            },
        )
        self.assertIn(
            name,
            {
                item["name"]
                for item in production.expected_remote_items(self.state, storage_record)
            },
        )
        self.seed_segmented_remote()
        self.assertTrue(self.state["records"][logical_record]["remote_verified"])
        self.assertTrue(self.state["records"][storage_record]["remote_verified"])

    def test_segmented_verify_download_reassembles_and_discards_segments(self) -> None:
        document = self.bind_segmented_transport()
        self.seed_segmented_remote()
        self.assertTrue(
            all(
                self.state["records"][record]["remote_verified"]
                for record in production.RECORDS
            )
        )
        download = self.base / "segmented-download"
        reconstructed = self.base / "segmented-paper"

        def fake_reconstruct(manifest, output, *, archive_root):
            self.assertEqual(Path(manifest), download / "manifest.json")
            self.assertEqual(Path(archive_root), download)
            for record in production.RECORDS:
                for archive in document["records"][record]["archives"]:
                    self.assertEqual(
                        (download / archive["relative_local_path"]).read_bytes(),
                        self.payloads[archive["name"]],
                    )
            Path(output).mkdir()
            return Path(output), Path(archive_root)

        with mock.patch.object(
            production.reconstruction, "reconstruct", side_effect=fake_reconstruct
        ):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=download,
                reconstructed_root=reconstructed,
            )
        self.assertEqual(self.state["remote_download"]["status"], "PASS")
        self.assertTrue((download / production.TRANSPORT_MANIFEST_FILENAME).is_file())
        self.assertFalse((download / ".transport-segments").exists())
        self.assertFalse(any(".chunk" in path.name for path in download.rglob("*")))
        expected_segments = sum(
            len(archive["segments"])
            for record in production.RECORDS
            for archive in document["records"][record]["archives"]
        )
        self.assertEqual(
            len(self.state["remote_download"]["completed_segments"]),
            expected_segments,
        )

    def test_segmented_downloads_archives_in_parallel_and_persists_on_main(
        self,
    ) -> None:
        self.bind_segmented_transport()
        self.seed_segmented_remote()
        download = self.base / "parallel-segmented-download"
        reconstructed = self.base / "parallel-segmented-paper"
        original_download = self.client.download_file
        release = threading.Event()
        lock = threading.Lock()
        active = 0
        maximum_active = 0

        def tracked_download(url, target, *, expected_size, expected_sha256) -> None:
            nonlocal active, maximum_active
            filename = url.rsplit("/", 1)[1]
            if ".chunk" not in filename:
                original_download(
                    url,
                    target,
                    expected_size=expected_size,
                    expected_sha256=expected_sha256,
                )
                return
            with lock:
                active += 1
                maximum_active = max(maximum_active, active)
                if active >= 2:
                    release.set()
            try:
                if not release.wait(timeout=2):
                    raise AssertionError("A second archive download never overlapped")
                original_download(
                    url,
                    target,
                    expected_size=expected_size,
                    expected_sha256=expected_sha256,
                )
            finally:
                with lock:
                    active -= 1

        self.client.download_file = tracked_download
        main_thread = threading.get_ident()
        persisted_threads = []
        factory_threads = []
        worker_clients = []
        original_persist = self.persist

        def guarded_persist(state) -> None:
            persisted_threads.append(threading.get_ident())
            original_persist(state)

        def client_factory():
            with lock:
                factory_threads.append(threading.get_ident())
                worker = mock.Mock()
                worker.download_file.side_effect = tracked_download
                worker_clients.append(worker)
            return worker

        with mock.patch.object(
            production.reconstruction,
            "reconstruct",
            side_effect=lambda manifest, output, *, archive_root: (
                Path(output).mkdir(),
                Path(archive_root),
            ),
        ):
            production.download_and_reconstruct(
                self.state,
                self.client,
                guarded_persist,
                download_root=download,
                reconstructed_root=reconstructed,
                download_workers=4,
                download_client_factory=client_factory,
            )
        self.assertGreaterEqual(maximum_active, 2)
        self.assertGreaterEqual(len(factory_threads), 2)
        self.assertNotIn(main_thread, factory_threads)
        self.assertEqual(
            len({id(client) for client in worker_clients}), len(worker_clients)
        )
        self.assertTrue(persisted_threads)
        self.assertEqual(set(persisted_threads), {main_thread})
        self.assertEqual(self.state["remote_download"]["status"], "PASS")
        self.assertFalse((download / ".transport-segments").exists())

    def test_parallel_segment_failure_persists_progress_and_resumes(self) -> None:
        transport = self.bind_segmented_transport()
        self.seed_segmented_remote()
        download = self.base / "parallel-resume-download"
        reconstructed = self.base / "parallel-resume-paper"
        archive = transport["records"]["software_core"]["archives"][0]
        first_name = archive["segments"][0]["name"]
        failure_name = archive["segments"][1]["name"]
        original_download = self.client.download_file
        lock = threading.Lock()
        failed = False
        failure_raised = threading.Event()
        calls = {}

        def fail_once(url, target, *, expected_size, expected_sha256) -> None:
            nonlocal failed
            filename = url.rsplit("/", 1)[1]
            with lock:
                calls[filename] = calls.get(filename, 0) + 1
                should_fail = filename == failure_name and not failed
                if should_fail:
                    failed = True
            if should_fail:
                failure_raised.set()
                raise production.ZenodoReleaseError("injected segment failure")
            if filename.startswith("main.part") and not failure_raised.wait(timeout=2):
                raise AssertionError("Software failure was not reached")
            original_download(
                url,
                target,
                expected_size=expected_size,
                expected_sha256=expected_sha256,
            )

        self.client.download_file = fail_once
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "injected segment failure"
        ):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=download,
                reconstructed_root=reconstructed,
                download_workers=2,
            )
        first_id = f"software_core/{first_name}"
        self.assertIn(first_id, self.state["remote_download"]["completed_segments"])
        part = download / "software_core" / f".{archive['name']}.part"
        self.assertEqual(part.stat().st_size, archive["segments"][0]["size_bytes"])
        self.assertEqual(calls[first_name], 1)
        self.assertGreaterEqual(
            len(self.state["remote_download"]["completed_files"]), 2
        )
        deferred_names = {
            segment["name"]
            for record in ("supplementary_figure_data", "figs3_wsi")
            for row in transport["records"][record]["archives"]
            for segment in row["segments"]
        }
        self.assertFalse(deferred_names & set(calls))

        with mock.patch.object(
            production.reconstruction,
            "reconstruct",
            side_effect=lambda manifest, output, *, archive_root: (
                Path(output).mkdir(),
                Path(archive_root),
            ),
        ):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=download,
                reconstructed_root=reconstructed,
                download_workers=4,
            )
        self.assertEqual(calls[first_name], 1)
        self.assertEqual(calls[failure_name], 2)
        self.assertEqual(self.state["remote_download"]["status"], "PASS")
        self.assertFalse((download / ".transport-segments").exists())

    def test_download_worker_bounds_fail_before_filesystem_changes(self) -> None:
        self.assertEqual(production.DEFAULT_DOWNLOAD_WORKERS, 4)
        parsed = production.parser().parse_args(
            [
                "--state",
                str(self.base / "state.json"),
                "download-roundtrip",
                "--download-root",
                str(self.base / "cli-download"),
                "--reconstructed-root",
                str(self.base / "cli-paper"),
            ]
        )
        self.assertEqual(parsed.download_workers, 4)
        for workers in (0, production.MAX_DOWNLOAD_WORKERS + 1, True):
            download = self.base / f"invalid-workers-{workers}"
            with self.assertRaisesRegex(
                production.ZenodoReleaseError, "Download workers must be between"
            ):
                production.download_and_reconstruct(
                    self.state,
                    self.client,
                    self.persist,
                    download_root=download,
                    reconstructed_root=self.base / f"invalid-paper-{workers}",
                    download_workers=workers,
                )
            self.assertFalse(download.exists())

    def test_state_rejects_archive_path_reused_across_records(self) -> None:
        duplicate = self.state["records"]["main_figure_data"]["files"][0][
            "relative_local_path"
        ]
        self.state["records"]["supplementary_figure_data"]["files"][0][
            "relative_local_path"
        ] = duplicate
        with self.assertRaisesRegex(
            production.ZenodoReleaseError,
            "Duplicate local archive path across records",
        ):
            production.validate_state(self.state)

    def test_segmented_verify_and_download_reject_corruption(self) -> None:
        self.bind_segmented_transport()
        self.seed_segmented_remote(verify=False)
        first = self.client.depositions[101]["files"][0]
        first["checksum"] = "md5:" + "0" * 32
        with self.assertRaisesRegex(production.ZenodoReleaseError, "size/MD5"):
            production.verify_remote(self.state, self.client, self.persist)

        first["checksum"] = "md5:" + digest(self.payloads[first["key"]], "md5")
        production.verify_remote(self.state, self.client, self.persist)
        self.payloads[first["key"]] = b"x" * len(self.payloads[first["key"]])
        with self.assertRaisesRegex(production.ZenodoReleaseError, "SHA-256"):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=self.base / "corrupt-segment-download",
                reconstructed_root=self.base / "corrupt-segment-paper",
            )

    def test_segmented_download_reconciles_interrupted_append(self) -> None:
        transport = self.bind_segmented_transport()
        self.seed_segmented_remote()
        download = self.base / "interrupted-segment-download"
        download.mkdir()
        self.state["remote_download"].update(
            {
                "status": "IN_PROGRESS",
                "root": str(download),
                "completed_files": [],
                "completed_segments": [],
            }
        )
        archive = transport["records"]["software_core"]["archives"][0]
        first = archive["segments"][0]
        part = download / "software_core" / f".{archive['name']}.part"
        part.parent.mkdir()
        part.write_bytes(self.payloads[first["name"]])
        cached = download / ".transport-segments" / "software_core" / first["name"]
        cached.parent.mkdir(parents=True)
        cached.write_bytes(self.payloads[first["name"]])

        with mock.patch.object(
            production.reconstruction,
            "reconstruct",
            side_effect=lambda manifest, output, *, archive_root: (
                Path(output).mkdir(),
                Path(archive_root),
            ),
        ):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=download,
                reconstructed_root=self.base / "interrupted-segment-paper",
            )
        self.assertEqual(
            (download / archive["relative_local_path"]).read_bytes(),
            self.payloads[archive["name"]],
        )
        self.assertFalse((download / ".transport-segments").exists())

    def test_segmented_download_rolls_back_mid_segment_append(self) -> None:
        transport = self.bind_segmented_transport()
        self.seed_segmented_remote()
        download = self.base / "mid-segment-download"
        download.mkdir()
        self.state["remote_download"].update(
            {
                "status": "IN_PROGRESS",
                "root": str(download),
                "completed_files": [],
                "completed_segments": [],
            }
        )
        archive = transport["records"]["software_core"]["archives"][0]
        first, second = archive["segments"][:2]
        part = download / "software_core" / f".{archive['name']}.part"
        part.parent.mkdir()
        second_payload = self.payloads[second["name"]]
        part.write_bytes(
            self.payloads[first["name"]]
            + second_payload[: max(1, len(second_payload) // 2)]
        )
        cached = download / ".transport-segments" / "software_core" / second["name"]
        cached.parent.mkdir(parents=True)
        cached.write_bytes(second_payload)
        downloaded_names = []
        original_download = self.client.download_file

        def track_download(url, target, *, expected_size, expected_sha256) -> None:
            downloaded_names.append(url.rsplit("/", 1)[1])
            original_download(
                url,
                target,
                expected_size=expected_size,
                expected_sha256=expected_sha256,
            )

        self.client.download_file = track_download
        with mock.patch.object(
            production.reconstruction,
            "reconstruct",
            side_effect=lambda manifest, output, *, archive_root: (
                Path(output).mkdir(),
                Path(archive_root),
            ),
        ):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=download,
                reconstructed_root=self.base / "mid-segment-paper",
            )
        self.assertNotIn(second["name"], downloaded_names)
        self.assertEqual(
            (download / archive["relative_local_path"]).read_bytes(),
            self.payloads[archive["name"]],
        )
        self.assertFalse((download / ".transport-segments").exists())

    def test_authenticated_download_checks_sha_and_reconstructs(self) -> None:
        self.draft_upload()
        download = self.base / "remote-download"
        reconstructed = self.base / "remote-paper"

        def fake_reconstruct(manifest, output, *, archive_root):
            self.assertEqual(Path(manifest), download / "manifest.json")
            self.assertEqual(Path(archive_root), download)
            Path(output).mkdir()
            return Path(output), Path(archive_root)

        with mock.patch.object(
            production.reconstruction, "reconstruct", side_effect=fake_reconstruct
        ):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=download,
                reconstructed_root=reconstructed,
            )
        self.assertEqual(self.state["remote_download"]["status"], "PASS")
        self.assertEqual(
            self.state["remote_download"]["manifest_sha256"],
            self.state["manifest"]["sha256"],
        )
        for record in production.RECORDS:
            for item in self.state["records"][record]["files"]:
                observed = download / item["relative_local_path"]
                self.assertEqual(observed.stat().st_size, item["size_bytes"])
                self.assertEqual(production.sha256_file(observed), item["sha256"])

    def test_download_rejects_corrupt_client_output(self) -> None:
        self.draft_upload()
        original = self.client.download_file

        def corrupt(url, target, *, expected_size, expected_sha256):
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b"x" * expected_size)

        self.client.download_file = corrupt
        with self.assertRaisesRegex(production.ZenodoReleaseError, "SHA-256"):
            production.download_and_reconstruct(
                self.state,
                self.client,
                self.persist,
                download_root=self.base / "bad-download",
                reconstructed_root=self.base / "bad-paper",
            )
        self.client.download_file = original

    def prepare_publish_gate(self) -> None:
        self.draft_upload()
        figure_log = self.base / "figures.log"
        figure_log.write_bytes(b"ALL FIGURES PASS\n")
        self.state["remote_download"].update(
            {
                "status": "PASS",
                "manifest_sha256": self.state["manifest"]["sha256"],
                "tree_sha256": self.state["manifest"]["tree_sha256"],
                "root": str(self.base / "download"),
                "reconstructed_root": str(self.base / "paper"),
            }
        )
        self.state["reproduction"] = {
            "status": "PASS",
            "returncode": 0,
            "log_path": str(figure_log),
            "log_sha256": production.sha256_file(figure_log),
            "command": ["reproduce_all_figures.sh", "--force"],
            "remote_manifest_sha256": self.state["manifest"]["sha256"],
            "remote_tree_sha256": self.state["manifest"]["tree_sha256"],
        }

    def test_reproduction_must_be_launched_from_remote_tree(self) -> None:
        root = self.base / "paper"
        root.mkdir()
        script = root / "reproduce_all_figures.sh"
        script.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        script.chmod(0o755)
        self.state["remote_download"].update(
            {
                "status": "PASS",
                "manifest_sha256": self.state["manifest"]["sha256"],
                "tree_sha256": self.state["manifest"]["tree_sha256"],
                "reconstructed_root": str(root),
            }
        )
        log = self.base / "figures.log"
        with self.assertRaises(production.ZenodoReleaseError):
            production.reproduce_roundtrip_figures(
                self.state, self.persist, log_path=log, execute=False
            )

        def runner(command, *, cwd, log):
            self.assertEqual(cwd, root)
            self.assertEqual(command[0], str(script))
            log.write(b"ALL FIGURES PASS\n")
            return 0

        production.reproduce_roundtrip_figures(
            self.state,
            self.persist,
            log_path=log,
            execute=True,
            runner=runner,
        )
        self.assertEqual(self.state["reproduction"]["status"], "PASS")
        self.assertEqual(
            self.state["reproduction"]["remote_manifest_sha256"],
            self.state["manifest"]["sha256"],
        )

    def test_publish_requires_roundtrip_phrase_and_current_state_digest(self) -> None:
        self.prepare_publish_gate()
        with self.assertRaisesRegex(production.ZenodoReleaseError, "confirm"):
            production.publish_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=True,
                confirmation="yes",
                expected_state_sha256=production.state_sha256(self.state),
            )
        with self.assertRaisesRegex(production.ZenodoReleaseError, "SHA-256"):
            production.publish_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=True,
                confirmation=production.PUBLISH_CONFIRMATION,
                expected_state_sha256="0" * 64,
            )
        expected = production.state_sha256(self.state)
        production.publish_all(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
            confirmation=production.PUBLISH_CONFIRMATION,
            expected_state_sha256=expected,
        )
        self.assertEqual(self.client.publishes, len(production.RECORDS))
        self.assertTrue(
            all(self.state["records"][name]["published"] for name in production.RECORDS)
        )

    def test_publish_is_blocked_without_remote_reproduction(self) -> None:
        self.draft_upload()
        before = self.client.publishes
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "roundtrip|downloads"
        ):
            production.publish_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=True,
                confirmation=production.PUBLISH_CONFIRMATION,
                expected_state_sha256=production.state_sha256(self.state),
            )
        self.assertEqual(self.client.publishes, before)

    def test_public_url_map_has_only_manifest_shards(self) -> None:
        self.prepare_publish_gate()
        production.publish_all(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
            confirmation=production.PUBLISH_CONFIRMATION,
            expected_state_sha256=production.state_sha256(self.state),
        )
        document = production.build_public_url_map(self.state, self.client)
        self.assertEqual(document["schema"], production.reconstruction.URL_MAP_SCHEMA)
        self.assertNotIn("manifest.json", document["archives"])
        self.assertEqual(len(document["archives"]), len(production.RECORDS))
        self.assertTrue(
            all(
                url.startswith("https://zenodo.org/records/")
                and url.endswith("?download=1")
                for url in document["archives"].values()
            )
        )

    def test_segmented_publish_preflight_and_v2_public_url_map(self) -> None:
        transport = self.bind_segmented_transport()
        self.seed_segmented_remote()
        figure_log = self.base / "segmented-figures.log"
        figure_log.write_bytes(b"ALL FIGURES PASS\n")
        self.state["remote_download"].update(
            {
                "status": "PASS",
                "manifest_sha256": self.state["manifest"]["sha256"],
                "tree_sha256": self.state["manifest"]["tree_sha256"],
                "transport_manifest_sha256": self.state["segmented_transport"][
                    "sha256"
                ],
                "root": str(self.base / "segmented-download"),
                "reconstructed_root": str(self.base / "segmented-paper"),
            }
        )
        self.state["reproduction"] = {
            "status": "PASS",
            "returncode": 0,
            "log_path": str(figure_log),
            "log_sha256": production.sha256_file(figure_log),
            "command": ["reproduce_all_figures.sh", "--force"],
            "remote_manifest_sha256": self.state["manifest"]["sha256"],
            "remote_tree_sha256": self.state["manifest"]["tree_sha256"],
            "remote_transport_manifest_sha256": self.state["segmented_transport"][
                "sha256"
            ],
        }
        production.publish_all(
            self.state,
            self.client,
            self.persist,
            execute=True,
            production=True,
            confirmation=production.PUBLISH_CONFIRMATION,
            expected_state_sha256=production.state_sha256(self.state),
        )
        document = production.build_public_url_map(self.state, self.client)
        self.assertEqual(
            document["schema"], production.reconstruction.URL_MAP_SCHEMA_V2
        )
        self.assertEqual(
            document["manifest"],
            {
                "sha256": self.state["manifest"]["sha256"],
                "size_bytes": self.state["manifest"]["size_bytes"],
            },
        )
        expected_archives = {
            archive["relative_local_path"]
            for record in production.RECORDS
            for archive in transport["records"][record]["archives"]
        }
        self.assertEqual(set(document["archives"]), expected_archives)
        for archive in document["archives"].values():
            self.assertEqual(set(archive), {"size_bytes", "sha256", "segments"})
            self.assertTrue(archive["segments"])
            for segment in archive["segments"]:
                self.assertEqual(set(segment), {"name", "size_bytes", "sha256", "url"})
                self.assertTrue(segment["url"].endswith("?download=1"))

    def test_segmented_publish_preflight_rejects_remote_extra(self) -> None:
        self.bind_segmented_transport()
        self.seed_segmented_remote()
        figure_log = self.base / "segmented-extra-figures.log"
        figure_log.write_bytes(b"ALL FIGURES PASS\n")
        self.state["remote_download"].update(
            {
                "status": "PASS",
                "manifest_sha256": self.state["manifest"]["sha256"],
                "tree_sha256": self.state["manifest"]["tree_sha256"],
                "transport_manifest_sha256": self.state["segmented_transport"][
                    "sha256"
                ],
            }
        )
        self.state["reproduction"] = {
            "status": "PASS",
            "returncode": 0,
            "log_path": str(figure_log),
            "log_sha256": production.sha256_file(figure_log),
            "remote_manifest_sha256": self.state["manifest"]["sha256"],
            "remote_tree_sha256": self.state["manifest"]["tree_sha256"],
            "remote_transport_manifest_sha256": self.state["segmented_transport"][
                "sha256"
            ],
        }
        self.client.depositions[101]["files"].append(
            {
                "key": "unexpected.bin",
                "size": 1,
                "checksum": "md5:" + digest(b"x", "md5"),
                "links": {
                    "self": "https://zenodo.org/api/files/bucket-101/unexpected.bin"
                },
            }
        )
        with self.assertRaisesRegex(
            production.ZenodoReleaseError, "file inventory differs"
        ):
            production.publish_all(
                self.state,
                self.client,
                self.persist,
                execute=True,
                production=True,
                confirmation=production.PUBLISH_CONFIRMATION,
                expected_state_sha256=production.state_sha256(self.state),
            )
        self.assertEqual(self.client.publishes, 0)

    def test_state_and_token_files_are_private_and_never_embed_credentials(
        self,
    ) -> None:
        state_path = self.base / "state.json"
        production.write_state(state_path, self.state)
        self.assertEqual(state_path.stat().st_mode & 0o777, 0o600)
        self.assertNotIn("example-production-token", state_path.read_text())
        state_path.chmod(0o644)
        with self.assertRaises(production.ZenodoReleaseError):
            production.load_state(state_path)
        token_path = self.base / "token"
        token_path.write_text("example-production-token\n", encoding="utf-8")
        token_path.chmod(0o600)
        self.assertEqual(
            production.load_token(token_path, environ={}),
            "example-production-token",
        )
        token_path.chmod(0o644)
        with self.assertRaises(production.ZenodoReleaseError):
            production.load_token(token_path, environ={})

    def test_production_host_validation_rejects_token_exfiltration_urls(self) -> None:
        with self.assertRaises(production.ZenodoReleaseError):
            production._safe_api_url("https://evil.example/api/files/x")
        with self.assertRaises(production.ZenodoReleaseError):
            production._safe_api_url("https://zenodo.org/api/files/x?access_token=leak")
        with self.assertRaises(production.ZenodoReleaseError):
            production._safe_public_url(
                "https://zenodo.org/records/1/files/x?token=leak"
            )

    def test_public_url_prefers_zenodo_html_over_external_doi_url(self) -> None:
        deposition = {
            "record_url": None,
            "doi_url": "https://doi.org/10.5281/zenodo.22265239",
            "links": {
                "record": "https://zenodo.org/api/records/22265239",
                "html": "https://zenodo.org/records/22265239",
            },
        }
        self.assertEqual(
            production._public_url(deposition),
            "https://zenodo.org/records/22265239",
        )


if __name__ == "__main__":
    unittest.main()
