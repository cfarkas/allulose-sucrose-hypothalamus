# Acquisition comparability and interpretation

Magnification alone does not invalidate measurements made within regions of
interest (ROIs). Comparability depends on the physical tissue sampled, optical
sectioning, signal detection and the endpoint being measured.

## What ROI quantification controls

A fraction such as c-FOS-positive/NPY-positive nuclei normalizes for the number
of sampled NPY-positive nuclei. Coordinates calibrated in micrometres allow
physical distances to be compared across pixel sizes. A tissue ROI limits the
sampled region; an individual-cell ROI defines the object being measured.
Neither operation alone makes marker detection or anatomical sampling identical.

If the same cells are resolved, classified and sampled comparably, changing
magnification need not change a labeling fraction. Conversely, a faint positive
cell missed by acquisition or segmentation stays missing after ROI normalization.
Numerator and denominator may be affected differently. The same number of pixels
also need not represent the same physical area across acquisitions.

| Endpoint | Useful normalization | Remaining acquisition dependence |
|---|---|---|
| Positive-cell fractions | Positive cells divided by the relevant sampled population | Staining sensitivity, false positives/negatives, overlapping objects and segmentation |
| Spatial clustering | Calibrated coordinates and a null retaining observed positions and label counts | Missed/merged cells, projection depth, sampled anatomy and spatially varying detection |
| Fluorescence intensity | Background correction and an explicitly defined area/cell summary | Illumination, detector response, collection efficiency, staining and saturation |
| Fine cell morphology | Physical scale and reviewed object boundaries | Resolution, optical sectioning and visibility of thin processes |

Optical resolution depends on numerical aperture and wavelength as well as
sampling; magnification labels alone are insufficient. Instrument and detector
settings also affect quantitative fluorescence. See [Waters, 2009](https://doi.org/10.1083/jcb.200903097)
and [Montero Llopis et al., 2021](https://doi.org/10.1038/s41592-021-01156-w).

## The NPY cohort

The documented objectives are 20×, 25× and 40×. No 35× objective is identified in
the NPY acquisition records.

| Animals | Acquisition | Analysis pixel pitch | NPY-associated spatial eligibility |
|---|---|---|---|
| WATER_NPY1–2 | Leica SP8, 20×/0.75 | Approximately 0.565–0.568 µm | Both included |
| WATER_NPY3 | Leica SP8, 40×, 21-plane maximum projection | Approximately 0.569 µm | Excluded: zero double-positive nuclei |
| WATER_NPY4–5 | Zeiss LSM 780, 25×/0.8, nine-plane maximum projections | 0.664213 µm after 2×2 mean reduction from 0.332106 µm | WATER_NPY4 excluded: one double-positive nucleus; WATER_NPY5 included: two |
| Three sucrose and three allulose animals | Zeiss LSM 780, 25×/0.8 | 0.664213 µm | All six included |

All five water animals contribute abundance measurements and the total c-FOS
spatial endpoint. NPY-associated clustering requires at least two double-positive
nuclei, giving three estimable animals in each condition. Thus, the 40× animal
cannot account directly for the reported NPY-associated spatial test: it does
not enter that test. It does contribute to abundance and total c-FOS analyses.

The 20× and 25× numerical apertures are close (0.75 and 0.8). This makes it
inappropriate to infer a large resolution discrepancy from magnification alone.
Pixel pitches are also much smaller than the tested 20–150 µm distances. Those
facts support using physically calibrated coordinates, but do not validate the
equivalence of label detection, optical depth or sampling across instruments.

WATER_NPY4–5 share the sugar images' instrument/objective and analysis pixel pitch. Resampling is scale harmonization, not an
optical or staining calibration. Maximum projections through different depths
can also change overlap and apparent co-localization. The conditional spatial
null controls observed positions and counts; it cannot recover unobserved cells.

The allulose–sucrose comparison shares the documented acquisition setup and is
less affected by instrument differences than comparisons with the mixed-platform
water group. This does not establish an absence of batch or biological-design
effects. Acquisition differences are a plausible confound, not demonstrated
proof that the observed pattern is an imaging artifact. Paired imaging of the
same tissue with each setup, or independent matched acquisitions with detection
validation, would quantify that bias more directly.

Sources: [acquisition scale manifest](https://github.com/cfarkas/allulose-sucrose-hypothalamus/blob/main/figures/Fig3/source_data/current_statistics/inputs/NPY_acquisition_manifest.csv),
[spatial eligibility](https://github.com/cfarkas/allulose-sucrose-hypothalamus/blob/main/figures/Fig3/source_data/current_statistics/spatial_eligibility.csv),
and the acquisition descriptions accompanying the manuscript. The same NPY
instrument concern should not automatically be transferred to the separate
POMC/glial acquisitions, which used the Keyence platform.

## Behavioral units and biological claims

- Figure 1 uses single-bottle **volume changes, mL/cage**, with observations on
  days 1, 3 and 6 and day 6 as the primary endpoint. Figure S4 uses weighed
  bottles and reports **g/cage/day**. These are different recorded quantities;
  no assumption of identical solution density is needed.
- Small cage counts restrict precision and attainable permutation p values.
  Exhaustive enumeration is still computationally exact; exchangeability and
  other design assumptions determine inferential validity.
- c-FOS staining describes nuclear protein expression associated with cellular
  stimulation; it does not measure firing rate, peptide release or satiety.
  ACTH/CLIP staining identifies POMC-related immunoreactivity through a shared
  precursor/product epitope, not peptide processing or release separately.
- No significant Iba1/GFAP intensity difference does not establish unchanged
  glial numbers or function. Predicted morphology classes describe shapes;
  agreement with human shape labels does not validate inflammatory activity.
- The NPY-associated spatial comparison does not meet the significance
  threshold after correction for multiple comparisons. Hepatic epithelial-like
  effects describe relative predicted phenotype composition, not absolute
  hepatocyte abundance.

These limits apply when interpreting the measured markers and model outputs.
