FIGURE S1 - MALE SINGLE-BOTTLE EXPERIMENT

The combined publication PDF and PNG are at the parent FigS1 directory's root.
Exact live A-D crops are in panels/, legends and the caption are in legends/, plotted
values/statistics/confidence intervals/source hashes/audit CSVs are in
source_data/, and the deterministic manifest, scientific audit, and this
README are in provenance/. Every publication artifact is a real file.

Panels A-D show the male consumption trajectory/endpoint and
body-weight trajectory/endpoint. Behavior-specific exact cage-label
permutation tests are retained because this is a small cage-randomized design;
pairwise tests are Holm-adjusted within sex and endpoint. No treatment-by-sex
interaction is claimed.

Reproduce from the Paper directory:

  cd /media/server/STORAGE/Apotome/Paper
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/FigS1/01_analyze_single_bottle_male.py
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/FigS1/02_make_figure_s1_male.py

The renderer validates a complete sibling staging tree and replaces only an empty or
previously marker-owned tree. It refuses symlinks, modified files, and untracked
content, so the recovered live FigS1 directory must not be used as
--output-dir; render to a fresh directory and promote reviewed leaf artifacts.
