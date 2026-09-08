FIGURE 3 — COMPLETE THREE-CONDITION FIGURE (8 SEPTEMBER 2026)

Figure_3_cFos_NPY.pdf and its 600-dpi PNG are the complete English master.
Individual A–H panels and their legends are available in English and Spanish.
There is no complete Spanish master in the final figure directory.

A: original water/sucrose/allulose microscopy and magnified merges.
B/C/D: water/sucrose/allulose marker-positive DAPI nuclei.
E/F: complete-field c-FOS/DAPI and c-FOS-positive/NPY-positive summaries.
G/H: c-FOS/DAPI and double-positive/DAPI spatial occurrence profiles.

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
requires PyMuPDF and Pillow. To finish a new certified seven-panel analytical
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
