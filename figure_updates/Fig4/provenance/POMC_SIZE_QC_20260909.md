POMC size quality control, 9 September 2026

Following visual review of small POMC objects, a post-segmentation size heuristic
was applied uniformly across all conditions and anatomical regions. Objects with
full area below the mean DAPI nuclear area in the same reconstructed section and
region were excluded. Assignment to region uses the object centroid. No c-FOS
status, treatment identity or statistical outcome enters the rule. The factor
1.0 was fixed before rerunning treatment statistics.

All 76 reconstructed sections were assessed: 3,403 processed POMC objects,
780 exclusions (552 ARC, 57 ME, 171 outside accepted regions). Original TIFFs,
segmentation files, accepted HIL anatomy and registration hashes are retained.
Filtered labels are separate derived TIFFs, with object-level JSON decisions.
POMC nuclear assignments and dependent statistics and spatial plots were rerun.
Mask objects are not equivalent to the number of DAPI nuclei assigned POMC.

ARC mean c-FOS/POMC fractions were 18.4%, 22.0%, 25.6% in water, sucrose,
allulose. The exact water–allulose contrast is nominal p=0.047619; global ANOVA
p=0.176454. Cage-mean water–allulose p=0.200. Minimum-area factors 0.5, 1.0
and 1.5 retain the direction, with nominal p=0.261905, 0.047619, 0.047619.
The 1.0-factor decision is a quality-control heuristic, not a validated somatic
identity classifier or an independently confirmed treatment effect.

The ten stored POMC training fields were additionally checked: 67 of 346 objects
were excluded in derived masks. Historical training masks and model weights
were not overwritten, models were not retrained, and loss curves remain those
of the historical training runs. Figure S7 explicitly identifies the later QC.

Reproduce analysis and sensitivity:
  python Fig4/02_analyze_pomc_cfos.py --output /tmp/pomc_new --no-auto-conda
  python Fig4/07_audit_pomc_size_sensitivity.py --analysis-dir /tmp/pomc_new
  python Fig4/05_analyze_spatial_distributions.py --analysis-dir /tmp/pomc_new --output-dir /tmp/pomc_new/spatial --dpi 600
  python Fig4/03_make_figure_4_pomc_cfos.py --analysis-dir /tmp/pomc_new --outdir /tmp/pomc_figure_new --figure-name Figure_4 --dpi 600
