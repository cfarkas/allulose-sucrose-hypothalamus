FIGURE 1 - SINGLE BOTTLE EXPERIMENT

The combined publication PDF and PNG are at the parent Fig1 directory's root. Exact live
panel crops are in panels/, legends and the caption are in legends/, plotted
values/statistics/confidence intervals/source hashes are in source_data/, and
the deterministic manifest, statistical audit, and this README are in
provenance/. Every publication artifact is a real file; symlinks are forbidden.

Reproduce from the Paper directory:

  cd /media/server/STORAGE/Apotome/Paper
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/Fig1/01_analyze_behavior_experiment1.py
  /home/server/anaconda3/envs/paper_apotome_repro/bin/python scripts/Fig1/02_make_figure_1_behavior.py

The renderer builds a complete sibling staging tree and replaces only an empty or
previously marker-owned tree. It refuses symlinks, modified files, and untracked
content, so the recovered live Fig1 directory must not be used as --output-dir;
render to a fresh directory and promote reviewed leaf artifacts.
