FIGURE 4 — HIL FINAL
======================

Status
------
Figure 4 is final and HIL-accepted. All 76 sections of the strict ARC/ME/VMN
HIL review are accepted: 75 experimental sections plus one
technical negative control. Every quantified section uses the HIL mask; no
automated or mixed ARC/ME/VMN region is present in any inferential value.

  Figure_4.pdf  68c5449efda64442e87e5df99bfaebf34fd502505fd183f314e98eba026a033c
  Figure_4.png  adcf46808ecf259c302df67fe37f44575d99199d6c20b740c85391a7308fb3e1
  provenance/PROMOTION_RECEIPT_20260827T131441Z.json records the current
  promoted files and the automatic pre-replacement archive.

The immutable HUMAN_FINAL_BUILD_RECEIPT.json and its 2026-08-25 addendum
document the preceding HIL-final generation. On 2026-08-27 all ten isolated
panels gained paired English/Spanish 600-dpi PDF/PNG outputs. In the current
generation, A-C identify FR722 S02 Allulose, D corrects the upper ventricular
lumen from the two accepted HIL ARC lobes, and E remains the exact A+B+C
microscopy merge with its registered-cell magnification. F keeps the complete
ungated transgenic NPY-GFP microscopy signal and shows natural POMC intensity
strictly inside accepted POMC ROIs. Its smaller POMC/NPY-GFP ROI cartoon keeps
marker ROIs outline-free and adds white dashed tissue/ventricle guides. I/J give
their heatmaps more horizontal room. The immediately preceding generation is
preserved in
legacy/pre_real_roi_masked_pomc_signal_20260827T131354Z/.

The earlier machine-provisional preview is NOT the figure. It stays in place,
untouched, as Figure_ARC_ME_multipanel.png, its recovered PDF, and
figure_bundle_v3, and it must not be cited as the HIL final.


Reproducing it
--------------
Everything Figure 4 needs ships in this directory, so the three documented
commands work in a freshly unpacked copy of the tree, with no arguments:

  cd <wherever you unpacked Paper>
  ./Fig4/00_create_conda_envs.sh
  export PAPER_PYTHON=<path it prints>
  PY="$PAPER_PYTHON"

  $PY Fig4/02_analyze_pomc_cfos.py
  $PY Fig4/05_analyze_spatial_distributions.py \
    --analysis-dir analyses/Fig4/results/human_final_run \
    --output-dir analyses/Fig4/results/human_final_run/spatial
  $PY Fig4/03_make_figure_4_pomc_cfos.py

The analyzer reads raw_data, validates the 76 accepted HIL reviews in
hil_review/human_final_20260824_v1 against the prepared review images in
hil_review/prepared_review_input_20260824, and writes
analyses/Fig4/results/human_final_run. The renderer reads that directory and
writes analyses/Fig4/figure/human_final_render. Both directories are generated,
so they are absent in a fresh copy; add --force to an analyzer rerun to replace
only its own output directory. Promotion into this directory is a separate,
deliberate copy.

The historical explicit analyzer/renderer form is in
provenance/HUMAN_FINAL_BUILD_RECEIPT.json. The analysis takes roughly fifteen
minutes and writes about 2.2 GB.

Fig4/04_run_pomc_cfos_staged.sh is the fail-closed staging driver. It refuses
--mode build for Figure 4 by design,
because the analyzer validates its receipts against reconstructed files inside
its own already-nonempty output and so requires replacement mode, which a
no-replacement driver will not grant. Its Figure 4 modes are prepare-hil,
review, check and render, all writing only under a --stage-root outside the
tree. The three commands above are how this figure is rebuilt.


Doing the HIL review again
----------------------------
The shipped HIL review is complete, so a reproduction never has to redo it. To
repeat it anyway, launch exactly one HIL reviewer process from the Paper root — never
two writers against the same output directory:

  $PY scripts/shared/01_annotate_regions.py \
    --dataset pomc \
    --input-root Fig4/hil_review/prepared_review_input_20260824/reconstructed_sections \
    --analysis-dir Fig4/hil_review/prepared_review_input_20260824 \
    --output-dir analyses/Fig4/hil/accepted_annotations \
    --host 127.0.0.1 --port 0

Point --output-dir at a new directory to keep the shipped HIL review intact. The
analyzer picks up analyses/Fig4/hil/accepted_annotations automatically when it
exists, and otherwise falls back to the shipped HIL review.

Only --output-dir has to be given. --input-root and --analysis-dir now resolve on
their own to analyses/Fig4/results/final_run when that working copy exists and to
the shipped Fig4/hil_review/prepared_review_input_20260824 otherwise, so the
reviewer starts in a freshly unpacked tree, which previously it could not do.

Replacement is the reviewer's default: loading an accepted section, drawing new
polygons and pressing "Save drawn subset + accept" rewrites that section's
JSON/TIFF receipt pair atomically under the current reviewer, session and source
hashes. Nothing is deleted. --lock-accepted restores the old refusal.


Make Figure 4, from the browser
-------------------------------
The reviewer's sidebar carries a "Make Figure 4" button, and it finishes the job.
When the analysis directory it was given is already complete it renders from it;
otherwise it runs Fig4/02_analyze_pomc_cfos.py over the accepted receipts first
and then renders, streaming both logs back into the page. That is the whole
chain, so a completed review reaches a figure without leaving the browser.

The render lands in a new timestamped directory under analyses/Fig4/figure and
its master graphics are Figure_4.pdf and Figure_4.png, the names this figure is
promoted and cited under. When the build finishes the page prints their full
paths. Earlier it rendered into a random directory under /tmp under the renderer's
own default name, Figure_ARC_ME_multipanel, so a build that had actually
succeeded produced nothing the reviewer could recognise, and pressing the button
looked like nothing happening at all. Nothing is written into this directory
either way: promotion stays a separate, deliberate copy.

A server keeps the code it was started with. If scripts/shared/01_annotate_regions.py
changes while a reviewer is open, stop that process and relaunch it; the page
prints STALE SERVER in the figure build panel when it detects this.

  --figure-analysis-dir   completed analysis to render from, or an absent
                          directory for the button's own analysis run to write.
                          Omitted, the button chooses
                          analyses/Fig4/results/human_final_run when that is
                          absent, else a fresh directory beside it under
                          analyses/Fig4/results.
  --figure-prepared-dir   prepared review images the accepted receipts are
                          validated against; defaults to --analysis-dir.
  --figure-input-root     raw input the analysis run reads; defaults to the
                          analyzer's own default, Fig4/raw_data.
  --figure-output-root    parent of the fresh, absent render directory; defaults
                          to analyses/Fig4/figure.
  --no-figure-analyze     refuse to analyze; render only, from a completed
                          --figure-analysis-dir.

The checkbox "Accept the machine default segmentation of every remaining section
first" turns pending native automated masks into accepted receipts in one pass.
It requires typing the exact attestation
"I REVIEWED THE MACHINE DEFAULT SEGMENTATIONS", records
machine_default_batch=true, review_granularity=batch_machine_default and that
attestation string in every receipt it writes, and refuses any composite whose
medial seam still needs image-by-image confirmation. Receipts written that way
are always distinguishable from image-by-image review in the durable provenance,
and they are not what produced the shipped figure.

On a terminal, Fig4/02_analyze_pomc_cfos.py asks the same question: reuse the
segmentations it found, open the reviewer and segment again, or make_figure_4
from them. --make-figure requests the render non-interactively and
--figure-output-root chooses where it lands.


How a review receipt is bound
-----------------------------
A receipt is accepted only when the reviewed reconstructed DAPI SHA-256, the
stitch-registration output DAPI hash, the accepted mask TIFF hash, the image
dimensions, and the sample/section/animal identities all match the section being
analysed. Those hashes are the identity.

The absolute directory a receipt was created in is not. A review is validated
from a different root whenever the tree is unpacked somewhere else, whenever the
analyzer writes a fresh output directory beside the immutable prepared review
input, and always for anyone who downloads this archive. Requiring the recorded
root turned a complete 76/76 review into 0/76 at analysis time, so the analyzer
and the renderer now compare the <animal>/<section>/<file> identity and log a
"relocated stage" warning naming both paths when the roots differ. Ten of the 76
receipts were drawn against a /tmp copy of the identical prepared images and are
reused this way; the other 66 name the in-tree copy.


Panel set
---------
The figure is A-J. Panels A, B and C are the isolated DAPI, c-FOS and POMC
channels of the corrected DAPI-registered representative field, each drawn in
the exact tint it contributes to the merge, so their clipped sum over the
identical pixels is Panel E. Panel D is the channel-separated Cellpose cartoon,
placed directly under the channels it depicts and ordered DAPI, c-FOS, POMC to
match A, B and C. A-C are labeled FR722 S02, Allulose; paired Spanish isolates
say Alulosa. Panel D's upper ventricular lumen is delineated row by row between
the medial boundaries of the two accepted HIL ARC lobes, extended to the
cropped image top and stopped before the first accepted HIL ME row. It uses no
automated wall CSV or marker intensity. Panel E is the exact A+B+C clipped-sum
microscopy composite with the registered POMC-positive/c-FOS-positive cell
magnification. Panels G-H are the animal-level c-FOS/DAPI and POMC-activation
results.

Panel E names its animal, section and condition above the A/B/C channel key. At
the lower left is the deterministic POMC-positive-cell magnification; the dashed
box and two leaders mark exactly the region it magnifies. No NPY-GFP intensity
or segmentation layer is present in E.

Panel F names NPY-GFP with GFP as the suffix. Its DAPI and c-FOS contributions
use the exact registered microscopy planes, and the real NPY-GFP TIFF uses a
deterministic median/99.7th-percentile display window and gamma 0.72 followed by
a monotone screen blend. The NPY-GFP microscopy is not ROI-gated: its complete
transgenic green signal is retained, including signal outside the NPY-GFP ROI
mask. The real registered POMC TIFF uses the same full-field normalization as
panels C/E, is reduced to display gain 0.72, and is then multiplied by the exact
accepted POMC mask. It therefore retains 254 natural intensity levels and
contributes at 11,759 signal-bearing pixels inside the 11,790-pixel mask, with
zero contribution outside it. No flat POMC paint is used. The 40%-wide black
cartoon inset fills the exact 62 POMC and 202 NPY-GFP ROI interiors without
marker-ROI outlines. Thin white dashed lines delimit the DAPI-supported tissue
exterior and HIL-ARC-derived ventricular lumen; no tissue/ARC/ME fill, ARC/ME
contour, halo, dilation, smoothing, interpolation, or resampling is present.

Channel colours are fixed paper-wide: DAPI blue, c-FOS magenta, NPY-GFP green. POMC
is amber, and was green until 2026-08-27, which collided with NPY. The medial
eminence contour is brick rather than orange so it cannot be read as amber POMC
signal, and the c-FOS/POMC double-positive class in the segmentation overlay is
white rather than yellow for the same reason. POMC uses a bright amber on the
black microscopy field and a darker amber for titles on white; one value cannot
be seen on both.

The promoted 2026-08-27 figure uses FR-722 and FR6-1 as the two animals
carrying both POMC and NPY-GFP, ten sections in all; the representative field is
one of them. On 2026-08-28 a separate, not-yet-integrated extension added
FR7-5 and FR8-1 so a future solution-blind NPY/POMC overlap endpoint can use
animal n=4 after their new masks pass QC. No nucleus in the promoted dataset is
positive for both POMC and NPY: they are two different arcuate populations, so
the cartoon maps where each sits rather than treating edge overlap as cell-level
co-expression; the two promoted marker masks share 137 edge pixels, which retain
both color components in the cartoon. Panels I and J are the exact spatial
occurrence analysis Figure 3 runs, applied to these reconstructed ARC/ME
sections by 05_analyze_spatial_distributions.py. Their heatmap columns receive
more internal width without stretching the fixed 6.15 x 3.12 panel canvas. The
representative field is FR722 S02, Allulose, selected deterministically.

Panel files follow the letters:
  A dapi, B cfos, C pomc, D cellpose_cartoon, E microscopy,
  F microscopy_npy_gfp,
  G cfos_dapi, H pomc_activation, I spatial_cfos, J spatial_cfos_pomc
Each has English and _spanish panels/*.pdf and panels/*.png at 600 dpi, plus an
English legends/*_LEGEND.txt. The superseded standalone-cartoon F files are under
legacy/pre_npy_gfp_EF_20260827T064813Z/retired_from_active/; the still earlier F
anatomy files are under legacy/pre_user_requested_EF_bilingual_20260827T060204Z/.

On 2026-08-26 the cartoon row moved under A/B/C and the panels were re-lettered
in reading order, so D, E and F changed meaning. Because panel files are named
after their content, the old-lettering files were left behind rather than
overwritten; they are archived in
legacy/superseded_panel_lettering_20260826/ and were removed from panels/ and
legends/ so that each panel letter names exactly one panel. Addendum 17 of
incident.txt records the change.

The still older Figure_ARC_ME_panel_C_cellpose_cartoon.png and the prior F
anatomy files were moved from the active panel directory into
legacy/pre_user_requested_EF_bilingual_20260827T060204Z/ without deletion. The
current cartoon is panel D, and the active panels directory now contains exactly
the ten English/ten Spanish A-J PNGs (with matching PDFs).


Raw and scientific status
-------------------------
raw_data holds the complete 20x acquisition and its Cellpose companions.
RAW_DATA_MANIFEST.csv SHA-256:
  97cab53d18b684970c5b2253d8659e1df2b503e02b08aa1a92bacd9ae43ab586
FR7-4 was published twice: the nested acquisition copy FR7-3/FR7-4 carried its
*_seg.npy masks and the flat FR7-4 copy carried only the TIFFs, so five FR7-4
sections were not quantifiable from the shipped tree at all. The 30 masks were
copied byte-exact into the flat path and verified against both the manifest and
the prepared overlay; FR7-4_SEG_SELF_CONTAINMENT_RECEIPT.json records every one.
Nothing was moved or deleted.

Segmentation is never rerun: every Cellpose mask the analysis needs ships beside
its image, and --segment-missing is off. The POMC and c-FOS models ship under
raw_data/FR6-1; the DAPI model does not, so --segment-missing is the one option
this tree cannot satisfy on its own.

The HIL-final run quantifies 231,625 per-cell rows over 76 reconstructed
sections and 15 animals, one of which is the excluded technical control. Cells
and sections are never n. Animal-level inference is allowed only where each
condition has at least three independent animals; the analysed conditions carry
Water 3, Sucrose 5 and Allulose 6 animals, and ME remains descriptive where that
gate fails. Historical HIL-reviewed statistics must never be hardcoded into a new
run.


NPY/POMC overlap extension — 2026-08-28
---------------------------------------
Two independently copied NPY-transgenic CZI acquisitions now live under:

  raw_data/npy_pomc_overlap_extension_20260706

They are FR7-5 Sucrose and FR8-1 Allulose. Eight channel-separated OME-TIFF
Z-stacks and eight native-resolution Cellpose MIPs are complete. Figure 4 names
AF546 POMC throughout. The source CZI files remain intact because the
reconstruction incident requires separate confirmation before destructive
moves.

The inputs are ready for Cellpose, but no new masks or n=4 overlap statistic are
claimed yet. After segmentation, the solution-blind overlap analysis will use
one value per animal across FR722, FR6-1, FR7-5 and FR8-1. Reproduce the MIPs
into an absent directory with:

  "$PAPER_PYTHON" Fig4/06_prepare_npy_pomc_overlap_cellpose_inputs.py \
    --output-dir /tmp/fig4_npy_pomc_cellpose_inputs

See raw_data/npy_pomc_overlap_extension_20260706/README.txt and
provenance/NPY_POMC_OVERLAP_EXTENSION_RECEIPT_20260828.json. The promoted
Figure_4.pdf/png have not been changed by this extension.


Safety
------
Never run the unsafe Figure 4 spatial/deletion test or any quarantined
live-writing workflow. Paper scripts must never access MinKNOW/Minion_Data or
storage_disk2. Superseded copies of the two numbered scripts are preserved
byte-exact under legacy/superseded_scripts_20260824_portability/.

POMC size QC (2026-09-09): all 76 reconstructed sections use a uniform
post-segmentation minimum area equal to mean DAPI nuclear area in that section
and region. Original images, source masks and accepted anatomy remain intact.
Derived POMC assignments and statistics are recomputed. This cohort underwent
familiarization and a 4-hour fast; the separate 16-hour ACTH/CLIP cohort is S5.
Run 07_audit_pomc_size_sensitivity.py after the analyzer for the 0.5/1/1.5-factor
sensitivity. See provenance/POMC_SIZE_QC_20260909.md for counts and interpretation.

10 SEPTEMBER 2026 — SPATIAL RING EXPLANATION
Panel K shows the six covariance-normalized DAPI-count quantile rings, with inner
rings 1–2 highlighted around a schematic third ventricle. Only the within-ring
percentage formula accompanies the drawing; methodological detail is in the
legend. The nuclei are illustrative, not experimental observations. Rings are constructed separately within each reconstructed section;
corresponding-ring numerator and denominator counts are summed over an animal
before calculating the six occurrence percentages. Complete English master and
individual A–K panels in English and Spanish are rebuilt by the standard renderer.
See scripts/shared/spatial_ring_cartoon.py and its mathematical verification tests.
