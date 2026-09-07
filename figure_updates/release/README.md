# Public release control

This directory prepares a public, checksum-bound release without contacting
Zenodo or GitHub. Remote actions are deliberately absent from the local control
program and `remote-action` always fails closed.

## Release objects

The public Paper tree is one reproducible object reconstructed from four
scientifically distinct, mutually linked records. ZIP parts are transport
shards within a record, not separate records. Main-figure, supplementary-
figure, and full-resolution Figure S3 WSI inputs are independently downloadable
reproduction subsets.

| Record key | Public object | Reason for the split |
| --- | --- | --- |
| `software_core` | Scripts, environments, source tables, HIL-reviewed masks, figures, legends, thesis/manuscript deliverables, provenance, and reconstruction tooling | Versioned software/reproduction package; expected below the default 50 GB quota |
| `main_figure_data` | Raw and acquisition-derived inputs for Figures 1–5 | Independently downloadable source-data object for the main figures |
| `supplementary_figure_data` | Raw and acquisition-derived inputs for Figures S1–S6, including Figure S3 L2 working-resolution exports, conversion receipts, SHA-256 manifests, and downstream analysis/QC outputs | Supplementary source-data object; explicitly excludes the 24 Figure S3 full-resolution L0 BigTIFF exports |
| `figs3_wsi` | The 24 Figure S3 full-resolution scanner-native L0 H&E BigTIFF exports | Dedicated CC BY 4.0 WSI object; matched L2 exports and receipts remain in `supplementary_figure_data` |

Every Zenodo record and its files are set to public/open. The proposed GitHub
repository is public. Figure 5 must retain its documented classifier-provisional
caveat, and Figure S3 HistoPLUS phenotypes must remain described as exploratory
model-predicted outputs rather than validated mouse cell identities.

The archive builder keeps each ZIP below 45 GB, safely below the current 50 GB
individual-file limit, and permits at most 100 uploaded files per logical
record. A record above the default 50 GB quota requires an approved allocation.
The builder also enforces the account-wide 150 GB additional allowance and
never treats transport shards as extra records to evade quota policy. Limit
references checked on 2026-09-01:

- <https://support.zenodo.org/help/en-gb/1-upload-deposit/80-what-are-the-size-limitations-of-zenodo>
- <https://help.zenodo.org/docs/deposit/manage-quota/>
- <https://help.zenodo.org/docs/deposit/manage-files/>
- <https://support.zenodo.org/help/en-gb/1-upload-deposit/68-should-i-upload-data-and-software-in-one-or-two-records>

## Files

- `build_archives.py` inventories the canonical public tree, assigns every file
  exactly once, produces deterministic ZIP64 shards, and records paths, modes,
  sizes, member names, and SHA-256 values in `manifest.json`.
- The generated `manifest.json` must be uploaded with `software_core` and
  mirrored in the public GitHub reconstruction repository; it is required
  before any shard can be verified or downloaded safely.
- `reconstruct_paper.py` verifies local shards or downloads an exact public-
  HTTPS shard map, then transactionally recreates the manifest-bound Paper tree
  without extracting unlisted or path-traversing ZIP members.
- `release_plan.json` binds the exact reviewed manuscript title and 18-author
  display list to the four records and to public visibility.
- `required_inputs.template.json` is a deliberately incomplete approval file.
  Copy it outside the committed release metadata and resolve every applicable
  field. Never put credentials in it.
- `release_control.py` validates metadata, decompresses and hashes every ZIP
  member, evaluates release gates, and can emit a credential-free local request
  document for a separately reviewed executor.
- `tests/` exercises deterministic reconstruction metadata, checksum tampering,
  public visibility, author identity, quota gates, credential non-disclosure,
  and remote-action refusal.

The canonical selection excludes development-only `.mypy_cache/`,
`__pycache__/`, `*.orig`, and `*.rej` artifacts in addition to the documented
machine-local working directories. These files are not scientific inputs or
reproduction outputs. It also excludes downloaded third-party files named
`full_text_open_access.pdf`: an index-provided open-access URL is not treated as
sufficient evidence for redistribution under this package's licenses. The
entire `literature_webscrap_30_08_2026/` working scrape is local-only because it
is not a figure input and includes third-party metadata. Other unused
exploratory inputs are also omitted by the exact-prefix archive contract; the
active paper analyses do not read them. Figure S3 uses the original within-slide
organ segmentation masks for tile selection, HistoPLUS inference, classical
features, and animal-level statistics. Required microscope tile/channel
stitching receipts for Figures 2–4 remain because they reconstruct the
acquisitions and are distinct from biological cross-sample alignment.

## Required decisions and authority

No public draft or upload is ready until all of these are supplied or confirmed:

1. Confirmation of the approved root license matrix: MIT for code and
   software-support files, and CC BY 4.0 for the authors' original documents
   and data. The mixed `software_core` record uses Zenodo's `other-open`
   identifier and points to `LICENSE`; all three data records use `cc-by-4.0`.
   The private approvals file must repeat these exact identifiers.
2. Approval from all authors and redistribution rights for the manuscript,
   figures, HIL-reviewed masks, main-figure data, supplementary-figure data,
   and the dedicated Figure S3 full-resolution histology WSI record.
3. Review of third-party material/model-output redistribution and confirmation
   that no sensitive or restricted data are present.
4. A release date, public-visibility acknowledgement, and a completed clean-room
   reproduction test.
5. Zenodo quota increases for each verified record above 50 GB. Create the
   drafts first, allocate storage in each draft's **Manage storage** screen,
   and bind the reviewed allocation to that exact deposition ID before upload.
   Upfront support approval is still required if the account-wide additional
   allowance is insufficient.
6. A GitHub owner and any Zenodo record IDs created by the later draft step.
7. Credentials with appropriate scopes, supplied only through environment
   variables: `ZENODO_SANDBOX_TOKEN`, `ZENODO_TOKEN`, and either `GH_TOKEN` or
   `GITHUB_TOKEN`. Sandbox and production Zenodo credentials are separate.

ORCIDs and affiliations may be added after author confirmation, but must never
be guessed. `Carlos Farkas*` is preserved exactly in the manuscript display
list; the service-facing creator name removes only the corresponding-author
asterisk.


## Local workflow

Validate the immutable title/author/public-visibility plan:

```bash
/home/server/anaconda3/bin/python release/release_control.py validate-plan
```

Run the release unit tests:

```bash
/home/server/anaconda3/bin/python -m unittest discover -s release/tests -p 'test_*.py' -v
```

Build the archive set outside the Paper tree, using an absent output path:

```bash
/home/server/anaconda3/bin/python release/build_archives.py \
  --paper-root "$PWD" \
  --output /ABSOLUTE/ABSENT/PATH/paper-public-archive-set \
  --jobs 10
```

Verify the completed set, including every decompressed member checksum:

```bash
/home/server/anaconda3/bin/python release/release_control.py verify-archives \
  --archive-root /ABSOLUTE/PATH/paper-public-archive-set
```

Reconstruct from a local verified archive set into an absent destination:

```bash
/home/server/anaconda3/bin/python release/reconstruct_paper.py \
  --manifest /ABSOLUTE/PATH/paper-public-archive-set/manifest.json \
  --archive-root /ABSOLUTE/PATH/paper-public-archive-set \
  --output /ABSOLUTE/ABSENT/PATH/Paper
```

To download public Zenodo files first, provide `zenodo-urls.json` with schema
`apotome-paper-archive-url-map-v1` and an `archives` object mapping every exact
manifest `relative_archive_path` to its public HTTPS file URL. The key set must
match the manifest exactly; placeholders, embedded credentials, and non-HTTPS
URLs are rejected.

```bash
/home/server/anaconda3/bin/python release/reconstruct_paper.py \
  --manifest /ABSOLUTE/PATH/manifest.json \
  --url-map /ABSOLUTE/PATH/zenodo-urls.json \
  --download-dir /ABSOLUTE/ABSENT/PATH/verified-zenodo-shards \
  --output /ABSOLUTE/ABSENT/PATH/Paper
```

The unchanged template is expected to report `BLOCKED`:

```bash
/home/server/anaconda3/bin/python release/release_control.py readiness \
  --archive-root /ABSOLUTE/PATH/paper-public-archive-set \
  --approvals release/required_inputs.template.json \
  --operation zenodo-production-draft
```

After a separately stored approvals file and credential environment pass every
gate, `prepare-request` may write a local JSON request. It still performs no
network action:

```bash
/home/server/anaconda3/bin/python release/release_control.py prepare-request \
  --archive-root /ABSOLUTE/PATH/paper-public-archive-set \
  --approvals /ABSOLUTE/PRIVATE/PATH/release-approvals.json \
  --operation zenodo-production-draft \
  --output /ABSOLUTE/ABSENT/PATH/zenodo-draft.request.json
```

The three data-record quota confirmations may remain `false` in this pre-draft request: no
deposition exists yet to receive an allocation. After draft creation, record
each actual allocation with `zenodo_production.py confirm-quota`; upload and
publication remain blocked until every record above 50 GB is bound to its
exact draft ID with sufficient allocated bytes.

`remote-action` cannot create a draft, upload a file, publish a record, create a
repository, or create a GitHub release. That separation is intentional: the
local archive, clean-room, rights, quota, metadata, and checksum gates must be
reviewed before any external state is changed.

## GitHub source update, 2026-09-07

The published v1.0.0 manifests and release_plan.json describe the archived ten
figures. Their DOI, chunk and archive hashes are kept fixed. The current working
source adds Figure S6 and corrects Figure 1 terminology.

`export_figure_updates.py --repository /path/to/current/github/clone` exports
actual Fig1/Fig5/FigS6 scripts and mirrors, shared utilities, small Fig1/S6 inputs,
and reference graphics into `figure_updates/`, with `figure-updates.json`. It
requires an absent update payload so an earlier export is preserved. Keep the
current GitHub bootstrap and guides when exporting; do not replace them with
an older local template. `figure_updates.py` verifies this payload and supports
a checked installation over the archived Paper tree with backups.

Future Zenodo file releases should version the software/core record and the
supplementary-data record (the latter now also holds FigS6/raw_data). The main
raw-data and full-resolution WSI payloads have no changes from this figure update.
The reviewed Fig5 campaign and portable training inputs are now frozen under
Fig5/microglial_review_data/ and exported to GitHub (about 74 MB in ZIPs). The
explicit freeze script is build_microglial_review_data.py; local analyses/
working artifacts remain excluded by the canonical archive contract. No new Zenodo version is published by the source
export. See https://help.zenodo.org/docs/deposit/manage-versions/.
