# POMC: current statistical inputs and results

These files reproduce the numerical statistics used in the 26 September 2026
manuscript revision. Coordinates, reviewed cell labels, tissue-side maps and
ROI maps are sufficient to recompute conditional pair moments and the spatial
envelope tests. Columns unused by these computations were omitted from the
portable cell table; `input_provenance.json` records the original hashes.
Geometry files are copied byte for byte, and their paths are relative here.

From the Paper or `figure_updates` directory, run:

```sh
python scripts/shared/recheck_current_statistics.py
```

Dependencies: NumPy, pandas, SciPy. The command compares every row and column of
the recomputed numerical tables with the supplied results, including all radii,
animal/cage sample sizes and multiplicity corrections. It writes to a temporary
directory by default. Use `--output-dir NEW_DIRECTORY` to retain results.

This begins after segmentation and cell classification. It does not recreate
raw microscopy, retrain models or establish the independence/exchangeability of
experimental units. Review `scale_resolved_plan.json` for inferential limits.
NPY uses animals (cage membership unknown); POMC uses equally weighted cages.
Both NPY spatial endpoints use a two-test Holm family; the marker-associated
result is nominal p=0.0285714, adjusted p=0.0571429. The POMC endpoints are both
nonsignificant. NPY abundance p values are unadjusted and exploratory.
