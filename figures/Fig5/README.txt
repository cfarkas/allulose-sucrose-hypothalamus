FIGURE 5 — REGIONAL GLIAL SIGNALS AND MICROGLIAL MORPHOLOGY

Figure 5 measures Iba1 and GFAP signals and describes predicted microglial
morphology in reviewed hypothalamic regions. Anatomy, DAPI segmentation and
DAPI-seeded Iba1 masks have been reviewed across 61 sections.

Panels
  A–B  Iba1 and GFAP microscopy.
  C    Microglial morphology classification with TinyMorphCNN.
  D    Predicted morphology-class composition.
  E–F  GFAP workflow and DAPI/Iba1/GFAP mask illustration.
  G    Regional Iba1 signal.
  H    Fraction assigned to activated/amoeboid morphology classes.
  I    Regional GFAP signal.

Cells and sections are subsamples; animal-level values enter the inferential
comparisons. The morphology classes are provisional shape assignments, not
independent evidence of functional glial activation. TinyMorphCNN is a
79,780-parameter study-specific network developed using morphometric
pseudo-labels. Figure S6 evaluates morphometric proposals and frozen DINOv2
features against a selected human-review sample. Details are in
source_data/deep_learning_implementation.md and its JSON companion.

Files
  Figure_5.pdf/.png  Complete English figure, with a 600-dpi PNG.
  panels/           Individual A–I panels in English and Spanish, PNG and PDF.
  legends/          Figure and panel captions in both languages.
  source_data/      Animal values, statistical tests and model methods.
  raw_data/         Section microscopy and segmentation inputs.
  hil_review/       Accepted anatomy and reviewed DAPI/Iba1 masks.
  microglial_review_data/  Portable crops, masks, labels and a saved candidate.

Reproduce from Paper/ using the paper environment and a fresh output directory:
  python Fig5/01_analyze_gfap_iba1_microglia.py \
    --reviewed-dapi-root Fig5/hil_review/dapi_final_20260830_v1 \
    --reviewed-iba1-root Fig5/hil_review/iba1_final_20260830_v1 \
    --output /tmp/fig5_build/analysis </dev/null
  python Fig5/02_make_figure_5_gfap_iba1_microglia.py \
    --analysis-dir /tmp/fig5_build/analysis \
    --output-dir /tmp/fig5_build/figure --dpi 600

The full microscopy and accepted masks come from the Zenodo reconstruction.
The portable review bundle can be used directly after cloning GitHub. In a
reconstructed Paper/ tree, verify it with:
  python Fig5/12_microglial_choices.py --validate-only
Restore the 203 saved human choices and their matching candidate with:
  python Fig5/12_microglial_choices.py --output-dir /tmp/microglial_choices

Optional independent review and training
  python Fig5/12_microglial_choices.py --do_microglial_choices \
    --output-dir /tmp/microglial_independent

Open the MICROGLIAL_REVIEW_URL printed in the terminal. The session provides
300 targets with no prefilled choices. Correct cell/soma masks as needed and
assign one of four morphology classes. Training requires 200–300 labels with
class and animal coverage; it starts automatically at 300 qualifying labels.
Use --review-dir to resume a session. Training and validation use separate
animals. Review candidates are separate outputs from the figure's fixed inputs.

From the GitHub root, prefix these direct Fig5 commands with figures/.
The full workflow also accepts python3 reproduce.py --do_microglial_choices.
Further review, CNN training, DINOv2 preparation and transfer evaluation are
implemented by the numbered scripts in this directory. Each input bundle is
verified against its size and SHA-256 manifest.
