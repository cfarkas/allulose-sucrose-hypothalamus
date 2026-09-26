# Figure S3 method sources

Checked 2026-08-28. Primary documentation or primary research implementations
were used for the technical workflow.

- [TumorQuantAI](https://github.com/cfarkas/tumorquantai): receipt-driven MDS
  conversion and WSI analysis design; MIT-licensed reference implementation.
- [LazySlide workflow cheatsheet](https://lazyslide.readthedocs.io/en/stable/reference/workflow-cheatsheet.html)
  and [first analysis](https://lazyslide.readthedocs.io/en/latest/getting-started/first-analysis.html):
  tissue shapes, physical-resolution tiling, and GPU segmentation workflow.
- [LazySlide paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13076205/): scalable
  whole-slide analysis framework and data model.
- [HistoPLUS](https://github.com/owkin/histoplus/) and the
  [LazySlide model implementation](https://github.com/rendeirolab/lazyslide-models/blob/main/src/lazyslide_models/segmentation/cellvit_family/histoplus.py):
  gated 20x model, 840-pixel tile constraint, and 14 non-background output
  classes. The weight license is CC-BY-NC-ND-4.0 and non-commercial.
- [HoVer-Net implementation](https://github.com/vqdang/hover_net) and
  [paper](https://www.sciencedirect.com/science/article/pii/S1361841519301045):
  horizontal/vertical-map nuclear instance postprocessing lineage.

Inference note: HistoPLUS is not validated here as a mouse normal-organ cell
classifier. The tumor-specific raw class is excluded from all Figure S3
classification/composition displays and inference. The 13 retained outputs are
treated as exploratory model-predicted phenotype components and paired with
classical H&E computer-vision sensitivity features. Cell-example insets are
shown only when the corresponding animal-level contrast passes global BH
q<0.05; their polygons are the actual model-output geometries, not illustrative
annotations.
