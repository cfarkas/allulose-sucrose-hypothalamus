FIGURE 3 — COMPLETE THREE-CONDITION FIGURE (8 SEPTEMBER 2026)

Figure_3_cFos_NPY.pdf and its 600-dpi PNG are the complete English master.
Individual A–I panels and their legends are available in English and Spanish.
There is no complete Spanish master in the final figure directory.

A: original water/sucrose/allulose microscopy and magnified merges.
B/C/D: water/sucrose/allulose marker-positive DAPI nuclei.
E/F: complete-field c-FOS/DAPI and c-FOS-positive/NPY-positive summaries.
H/I: c-FOS/DAPI and double-positive/DAPI spatial occurrence profiles.
G: schematic DAPI nuclei around the third ventricle, with six quantile rings.

Sucrose E7_FR7-5 uses its original native masks and registered display transform.
Carlos Farkas accepted the ventricular polygon in the local HIL session:
hil_review/sucrose_ventricle_20260908_v1/.
The gray external tissue border uses the same DAPI-supported display-envelope
procedure as Water and Allulose. Neither border defines an analysis ROI.
Native counts: DAPI 3241; c-FOS 81; NPY 302; double-positive 19.
Display-crop counts: NPY 292; c-FOS 66; double-positive 19.

Rebuild this presentation from the current audited analytical outputs:
  FIG3_RENDER_PYTHON=/path/to/scientific/python \
  FIG3_PDF_PYTHON=/path/to/python/with/PyMuPDF \
  ./Fig3/07_rebuild_complete_figure3.sh

The renderer requires numpy, scipy, cv2, skimage and matplotlib; PDF assembly
requires PyMuPDF, Pillow, Matplotlib, NumPy, pandas, SciPy and scikit-image. To finish a new certified seven-panel analytical
build, pass --base-panels-dir /path/to/certified/final/panels. The script retains
its D/E quantitative sources and relabels them E/F in the complete presentation.
The main reproduce_all_figures.sh now runs this finishing step after installing
fresh analyses and before final validation. It never edits raw data or accepted
HIL masks. The older sealed analytical build remains unchanged.

Current presentation receipt:
  provenance/complete_conditions_20260908.json
Current legends:
  legends/Figure_3_cFos_NPY_LEGEND.txt (also Figure_3_cFos_NPY_caption.txt)
  legends/Figure_3_cFos_NPY_LEGEND_spanish.txt

Animal-level statistics and source_data measurements remain unchanged by this
presentation update. Historical promotion receipts document earlier A–G layouts
and remain under provenance/ as superseded records.

10 SEPTEMBER 2026 — SPATIAL RING EXPLANATION
Panel G shows illustrative DAPI nuclei around the third ventricle, six quantile
rings with the inner two highlighted, and one within-ring occurrence equation identifying H and I.
E and F occupy the left of the same row; spatial profiles are H and I below.
The illustration is not experimental data. Inner rings 1–2 contain approximately
one-third of the eligible DAPI nuclei. One field contributes per animal.
The shared renderer is scripts/shared/spatial_ring_cartoon.py. Figure 3 completion
regenerates this panel and corrects the historical density wording in its plot axes.
The certified historical analytical pipeline and its masks remain unchanged.

The ring schematic in G uses a narrow superior third-ventricle neck, an inferior
flare and rounded floor, saturated blue DAPI points and stronger ring colors.
It is the only ring-definition cartoon in the main figures; Figure 4 refers here.

Ring schematic placement, 11 September 2026
------------------------------------------
Figure 3G positions the illustrative nuclear distribution near the floor of the
third ventricle, placing rings 1–2 at that level. All six rings remain concentric
around the mean of the displayed DAPI nuclei and contain 60 illustrative nuclei
each. The ventricular floor is a drawing landmark, not an analytical origin.
The shared spatial_ring_cartoon.py renderer supplies the standalone schematic,
the complete Figure 3 builder, bilingual panels, and manuscript legends. Figure 4
refers to Figure 3G and does not repeat the cartoon.

Schematic relative size, 11 September 2026
----------------------------------------
The third-ventricle contour in Figure 3G is enlarged twofold about its floor.
The rings appear approximately half as large relative to the ventricle, with
ring numbers spread across the lower arc for readability. Illustrative DAPI
nuclei are regenerated outside the enlarged lumen and still define six equal-
count rings around their mean. This is a schematic scaling change only.

Connection to spatial profiles, 11 September 2026
The schematic identifies the DAPI mean with + and continues each covariance
ellipse across the empty illustrative lumen using dashed lines. The ventricle
does not define the analysis origin or ring boundaries. A single equation explicitly
maps H to all c-FOS-positive nuclei and I to c-FOS/NPY double-positive nuclei;
both use all DAPI nuclei in the same ring. The numerator identifies the endpoint (H or I); the DAPI denominator is
shown once. The legend explains individual-animal percentages and arithmetic means.
All 108 stored percentages agree with their counts; H/I use identical ring
boundaries and denominators. The nucleus counts of each animal's six rings
differ by at most one. These checks preserve the experimental results.

Final middle-row layout, 11 September 2026
E and F each receive 34% of the available row width and G receives 32%
(previously 26%, 26%, 48%). The schematic uses a compact 7.2 x 5.4 inch
vector canvas with one formula and a two-line numerator. H/I labels within
the schematic are retained by the assembler. Data panels are not cropped
or stretched; their standalone outputs and underlying measurements are intact.
