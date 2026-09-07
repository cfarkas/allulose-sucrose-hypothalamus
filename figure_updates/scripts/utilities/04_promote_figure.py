#!/usr/bin/env python3
"""Promote a finished render into its live FigN directory.

Renderers are fresh-output only and refuse to write into a live figure
directory, because publication there is an atomic directory replacement and a
figure directory also holds raw_data, the numbered scripts, README.txt and
legacy. So a rebuild lands in a working directory and this brings it home: the
master graphics under the name the paper cites them by, plus panels, legends,
source data and provenance.

Figure 2's renderer places legends, source_data and receipts beside its
``outputs`` directory. When that documented outputs directory is promoted,
those siblings are mapped to legends, source_data and provenance respectively.

It only ever adds or replaces the specific files it is promoting. By default,
anything it would replace whose bytes differ is copied to
FigN/legacy/superseded_promotion_* first, files that already match are skipped,
and nothing is ever deleted. The root all-figure launcher passes
``--replace-current`` after a complete validated rebuild. That mode atomically
rewrites every planned canonical artifact, including byte-identical files, and
stamps the installed copy with the current installation time. It deliberately
does not create another in-tree legacy copy because the canonical tree is being
sanitized and the pre-rebuild archive is maintained separately. Every mode
finishes by verifying source and destination SHA-256 values byte for byte.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve()
PAPER = next((parent for parent in HERE.parents
              if (parent.name == "Paper" or ((parent / "scripts" / "setup").is_dir()
                  and (parent / "README.txt").is_file()))), None)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")

# The basename each figure is promoted and cited under.
PROMOTED_MASTER = {
    "Fig1": "Figure_1_behavior_experiment_1",
    "Fig2": "Figure_2",
    "Fig3": "Figure_3_cFos_NPY",
    "Fig4": "Figure_4",
    "Fig5": "Figure_5",
    "FigS5": "Figure_S5",
    "FigS6": "Figure_S6",
    "FigS1": "Figure_S1_single_bottle_male",
    "FigS2": "Figure_S2_single_bottle_female",
    "FigS3": "Figure_S3",
    "FigS4": "Figure_S4_behavior_preference",
}
# What the renderer calls the master when that differs from the promoted name.
RENDER_MASTER = {
    "Fig4": ("Figure_4", "Figure_ARC_ME_multipanel"),
    "Fig5": ("Figure_5", "Figure_GFAP_Iba1_microglia_ARC_ME"),
}
SUBTREES = ("panels", "legends", "source_data", "provenance")
# Renderer-specific sibling layouts used by the documented promotion command.
# The source side is relative to ``source.parent`` and the target side is
# relative to FigN/. Figure 2 calls its provenance directory ``receipts``.
SIBLING_SUBTREES = {
    "Fig2": (
        ("legends", "legends"),
        ("source_data", "source_data"),
        ("receipts", "provenance"),
    ),
}
# Figure 3's promoted bytes are asserted from inside its own frozen reproduction
# bundle, by Fig3/reproduce_final_20260823/reference/STAGE_POSTFLIGHT.json, which
# is itself covered by that bundle's SHA256_MANIFEST.tsv. Its driver prints
# [PASS] when a rerun matches, and that is the verification. Overwriting the
# promoted files with a rerun's would leave the shipped figure disagreeing with
# the reference the bundle asserts, and Figure 3's PDF is not byte-stable across
# runs anyway; only its PNG is pixel-bound.
CERTIFIED_AGAINST_FROZEN_REFERENCE = {"Fig3"}
MASTER_SUFFIXES = (".pdf", ".png")
# Multipanel masters are English-only. Spanish is retained only for isolated
# panels under ``panels`` and those files are still promoted by SUBTREES.
MASTER_VARIANTS = {
    "FigS5": ("",),
}
# Never touched, whatever a render directory happens to contain.
PROTECTED = ("raw_data", "raw", "scripts", "hil_review", "legacy", "README.txt")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_install(src: Path, dst: Path, *, stamp_current: bool) -> None:
    """Install ``src`` at ``dst`` atomically, optionally refreshing its mtime."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=dst.parent, prefix=f".{dst.name}.promotion-", delete=False
        ) as handle:
            temporary = Path(handle.name)
        shutil.copy2(src, temporary)
        if stamp_current:
            os.utime(temporary, None)
        os.replace(temporary, dst)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def find_master(source: Path, figure: str) -> str:
    """The basename the renderer actually used for its master graphics."""
    for candidate in RENDER_MASTER.get(figure, (PROMOTED_MASTER[figure],)):
        if any((source / f"{candidate}{suffix}").is_file() for suffix in MASTER_SUFFIXES):
            return candidate
    raise SystemExit(
        f"No master graphic in {source}. Expected one of "
        + ", ".join(f"{c}.pdf/.png" for c in RENDER_MASTER.get(
            figure, (PROMOTED_MASTER[figure],)))
    )


def planned_copies(source: Path, target: Path, figure: str) -> list[tuple[Path, Path]]:
    render_name = find_master(source, figure)
    promoted_name = PROMOTED_MASTER[figure]
    pairs: list[tuple[Path, Path]] = []
    for variant in MASTER_VARIANTS.get(figure, ("",)):
        for suffix in MASTER_SUFFIXES:
            src = source / f"{render_name}{variant}{suffix}"
            if src.is_file():
                pairs.append((src, target / f"{promoted_name}{variant}{suffix}"))
    for extra in sorted(source.glob(f"{render_name}_caption.txt")):
        pairs.append((extra, target / f"{promoted_name}_caption.txt"))
    for name in SUBTREES:
        subdir = source / name
        if not subdir.is_dir():
            continue
        for src in sorted(subdir.rglob("*")):
            if src.is_file():
                pairs.append((src, target / name / src.relative_to(subdir)))
    if source.name == "outputs":
        for source_name, target_name in SIBLING_SUBTREES.get(figure, ()):
            subdir = source.parent / source_name
            if not subdir.is_dir():
                continue
            for src in sorted(subdir.rglob("*")):
                if src.is_file():
                    pairs.append((src, target / target_name / src.relative_to(subdir)))
    return pairs


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", required=True, choices=sorted(PROMOTED_MASTER))
    parser.add_argument("--from", dest="source", type=Path, required=True,
                        help="Finished render directory to promote from.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Report what would change and write nothing.")
    parser.add_argument(
        "--force", "--replace-certified", dest="replace_certified",
        action="store_true",
        help=("Promote a figure whose promoted bytes are asserted by a frozen "
              "reproduction bundle. Only Figure 3 is in that state."),
    )
    parser.add_argument(
        "--replace-current",
        action="store_true",
        help=("Atomically rewrite every current publication file, including "
              "byte-identical files, without copying old bytes into FigN/legacy. "
              "Intended only for the root launcher after every staged figure "
              "has passed validation."),
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help=("Write nothing; require every planned staged artifact to exist in "
              "the canonical figure package with the same SHA-256."),
    )
    parser.add_argument(
        "--verify-reference-images",
        action="store_true",
        help=("Write nothing; require regenerated master/panel PNG files to "
              "match the canonical repository references byte for byte."),
    )
    args = parser.parse_args()

    if (args.figure in CERTIFIED_AGAINST_FROZEN_REFERENCE
            and not args.replace_certified and not args.verify_only
            and not args.verify_reference_images):
        raise SystemExit(
            f"{args.figure} is certified against a frozen reproduction bundle and is "
            f"not promoted from a rerun. Its own driver, "
            f"Fig3/reproduce_final_20260823/RUN_REPRODUCE_FINAL.sh, verifies a rerun "
            f"against the frozen references and prints [PASS]; that is the check, and "
            f"the promoted files stay as they are. Pass --force only if "
            f"you intend to supersede the certified artifact."
        )

    source = args.source.expanduser().resolve()
    target = (PAPER / args.figure).resolve()
    if not source.is_dir():
        raise SystemExit(f"Render directory not found: {source}")
    if source == target or target in source.parents:
        raise SystemExit(f"Refusing to promote a figure directory onto itself: {source}")
    if not target.is_dir():
        raise SystemExit(f"Figure directory not found: {target}")

    pairs = planned_copies(source, target, args.figure)
    if args.verify_reference_images:
        image_pairs: list[tuple[Path, Path]] = []
        for src, dst in pairs:
            relative = dst.relative_to(target)
            if (dst.suffix.lower() == ".png"
                    and (len(relative.parts) == 1 or relative.parts[0] == "panels")):
                image_pairs.append((src, dst))
        if not image_pairs:
            raise SystemExit(f"No master or panel PNG artifacts planned for {args.figure}")
        failures: list[str] = []
        for src, dst in image_pairs:
            if not dst.is_file():
                failures.append(f"missing reference: {dst.relative_to(PAPER)}")
            elif sha256(src) != sha256(dst):
                failures.append(f"reference SHA-256 mismatch: {dst.relative_to(PAPER)}")
        if failures:
            for failure in failures:
                print(f"  {failure}")
            raise SystemExit(
                f"Repository image reproducibility failed for {args.figure}: "
                f"{len(failures)}/{len(image_pairs)} PNG artifacts"
            )
        print(
            f"[PASS] {args.figure} shipped-image byte reproduction: "
            f"{len(image_pairs)}/{len(image_pairs)} PNG artifacts"
        )
        return 0
    if args.verify_only:
        failures: list[str] = []
        for src, dst in pairs:
            if not dst.is_file():
                failures.append(f"missing: {dst.relative_to(PAPER)}")
            elif sha256(src) != sha256(dst):
                failures.append(f"SHA-256 mismatch: {dst.relative_to(PAPER)}")
        if failures:
            for failure in failures:
                print(f"  {failure}")
            raise SystemExit(
                f"Byte-for-byte canonical verification failed for {args.figure}: "
                f"{len(failures)}/{len(pairs)} artifacts"
            )
        print(
            f"[PASS] {args.figure} byte-for-byte canonical match: "
            f"{len(pairs)}/{len(pairs)} artifacts"
        )
        return 0

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = target / "legacy" / f"superseded_promotion_{stamp}"

    added, replaced, refreshed, unchanged, preserved = [], [], [], [], []
    for src, dst in pairs:
        relative = dst.relative_to(target)
        if relative.parts[0] in PROTECTED:
            raise SystemExit(f"Refusing to promote into a protected path: {relative}")
        if not dst.exists():
            added.append((src, dst))
        elif sha256(src) == sha256(dst):
            if args.replace_current:
                refreshed.append((src, dst))
            else:
                unchanged.append((src, dst))
        else:
            replaced.append((src, dst))
            preserved.append((dst, archive / relative))

    print(f"Promoting {args.figure} from {source}")
    print(
        f"  {len(added)} new, {len(replaced)} content-replaced, "
        f"{len(refreshed)} byte-identical rewritten, "
        f"{len(unchanged)} already identical skipped"
    )
    for _src, dst in added + replaced:
        print(f"    {dst.relative_to(PAPER)}")
    if preserved and not args.replace_current:
        print(f"  {len(preserved)} existing files first copied to "
              f"{archive.relative_to(PAPER)}")
    elif preserved:
        print(f"  {len(preserved)} existing files replaced in canonical-current mode; "
              "no in-tree legacy copy created")
    if args.dry_run:
        print("\n[DRY RUN] nothing was written")
        return 0
    if not added and not replaced and not refreshed:
        print("\n[OK] already promoted; nothing to do")
        return 0

    if not args.replace_current:
        for old, keep in preserved:
            keep.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(old, keep)
    written = added + replaced + refreshed
    for src, dst in written:
        atomic_install(src, dst, stamp_current=args.replace_current)

    verification_failures = [
        dst for src, dst in pairs
        if not dst.is_file() or sha256(src) != sha256(dst)
    ]
    if verification_failures:
        for dst in verification_failures:
            print(f"  byte mismatch after install: {dst.relative_to(PAPER)}")
        raise SystemExit(
            f"Post-install byte verification failed for {args.figure}: "
            f"{len(verification_failures)}/{len(pairs)} artifacts"
        )
    print(
        f"  [PASS] byte-for-byte post-install match: "
        f"{len(pairs)}/{len(pairs)} artifacts"
    )

    receipt = {
        "schema": "figure_promotion_receipt_v1",
        "figure": args.figure,
        "promoted_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_render_dir": str(source),
        "master_name_in_render": find_master(source, args.figure),
        "master_name_promoted": PROMOTED_MASTER[args.figure],
        "replacement_policy": (
            "validated_canonical_current_no_in_tree_legacy"
            if args.replace_current else "copy_superseded_to_figure_legacy"
        ),
        "superseded_copied_to": (
            str(archive.relative_to(PAPER))
            if preserved and not args.replace_current else None
        ),
        "counts": {"added": len(added), "replaced": len(replaced),
                   "rewritten_identical": len(refreshed),
                   "already_identical_skipped": len(unchanged)},
        "byte_for_byte_verified": True,
        "verified_artifacts": len(pairs),
        "files": [
            {"source": str(src), "promoted": str(dst.relative_to(PAPER)),
             "bytes": dst.stat().st_size, "sha256": sha256(dst)}
            for src, dst in sorted(written, key=lambda pair: str(pair[1]))
        ],
    }
    receipt_path = target / "provenance" / f"PROMOTION_RECEIPT_{stamp}.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    print(f"\n[OK] {args.figure} promoted; receipt: {receipt_path.relative_to(PAPER)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
