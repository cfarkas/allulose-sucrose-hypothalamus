# Project-trained model archive

Release: [models-2026-09-14](https://github.com/cfarkas/allulose-sucrose-hypothalamus/releases/tag/models-2026-09-14)

This archive contains **21 unique project-trained model binaries** and **13 saved Cellpose training-loss arrays** (7,924,959,764 bytes in total). It includes the checkpoints used for the manuscript analyses and the earlier training iterations and development candidates found in the project directories. Each binary retains its original bytes. The release stores the large files outside Git history; `manifest.json` records their sizes, SHA-256 checksums, roles, download URLs and suggested restoration paths.

The associated analysis code is identified by commit [`1962f75bbbaed8760b0cd981fd3c064678555e8e`](https://github.com/cfarkas/allulose-sucrose-hypothalamus/tree/1962f75bbbaed8760b0cd981fd3c064678555e8e). This model archive supplements the existing reproduction package and does not alter its figures or analysis outputs.

## Which models belong to which analysis

- **Custom Cellpose segmentation:** `dapi_channel`, `cfos_channel3`, `POMC_2`, `cfos_Keyence1` and `dapi_Keyence4` are the five manuscript checkpoints. Eight additional checkpoints preserve earlier project training iterations. The `cellpose__` prefix is a release naming convention; the manifest restores the original native model names under `models/cellpose/`. Each checkpoint has an associated saved training-loss array. The saved arrays are training records, not independent validation scores.
- **TinyMorphCNN treatment analysis:** `treatment_tinymorphcnn.pt` is the retained Figure 5 checkpoint. It was trained using high-confidence morphology-derived pseudo-labels. Its four predicted class labels are Ramified, Rod-like, Activated and Amoeboid. These are the class names stored in the model.
- **Reviewed TinyMorphCNN development:** `reviewed_initial_tinymorphcnn.pt`, `reviewed_v2_tinymorphcnn.pt` and `reviewed_v3_tinymorphcnn.pt` preserve successive supervised development runs using 203 reviewed cells from 15 animals, with 163 training cells and 40 held-out development cells. They are successive runs, not cross-validation folds. The v3 reviewed candidate did not replace the Figure 5 treatment checkpoint.
- **Morphometry/DINOv2 classifier candidates:** the four `.joblib` bundles contain final candidates fitted on **all 203 human-reviewed cells from 15 animals** after grouped model selection. They are not the fixed unsupervised morphometric comparator and are not serialized held-out fold models. `morphology_forest.joblib` and the project's `selected_classifier.joblib` are byte-identical; only the former is uploaded, and the manifest records the latter as a restoration alias.

The separate benchmark's held-out performance must be read from its saved out-of-fold predictions and validation design, rather than obtained by evaluating these final full-data fits on their training cells. In the reconstructed Paper project, the relevant records are under `analyses/Fig5/results/transfer_reviewed_20260907_v1/validation/`: `morphometric_proposals_out_of_fold.csv`, `dino_logistic_out_of_fold.csv`, `nested_selected_out_of_fold.csv`, `validation_design.json` and `model_selection_audit.json`. The manifest records these paths. No serialized inner or outer cross-validation fold checkpoints were found. Publishing these candidates neither changes the treatment results nor converts the separate benchmark into direct validation of TinyMorphCNN.

## Download and verify

The download helper lists available assets by default. Run it from the repository root:

```bash
python3 models/download_models.py --list
python3 models/download_models.py --asset cellpose__cfos_channel3 --destination ./model-downloads
python3 models/download_models.py --asset morphology_forest.joblib --destination ./model-downloads
python3 models/download_models.py --all --destination ./model-downloads
```

`--all` explicitly downloads all 34 assets, approximately 7.925 GB. Downloads require a POSIX system such as Linux. The helper verifies file sizes and SHA-256 checksums, restores each asset at its manifest path beneath the destination directory, and creates the `selected_classifier.joblib` alias when downloading `morphology_forest.joblib`. It refuses symbolic links and existing files with different contents. It downloads and verifies files without loading models or running inference.

Alternatively, download individual files from the release page, or retrieve all release assets with GitHub CLI:

```bash
gh release download models-2026-09-14   --repo cfarkas/allulose-sucrose-hypothalamus   --dir model-downloads
```

The `SHA256SUMS.txt` file in this directory lists checksums for all 34 model and loss assets. For the flat directory produced by the GitHub CLI command above, verify the files from the repository root with:

```bash
cd model-downloads
sha256sum -c ../models/SHA256SUMS.txt
```

`manifest.json` uses `restore_path` for a suggested path relative to the repository root and `project_path`, where present, for the original location relative to the reconstructed Paper project. These are different path bases. Downloading model files does not automatically change the analysis configuration or replace generated project outputs. The TinyMorphCNN analyzer accepts a preserved checkpoint through `--cnn-model-in`; its implementation is in `figure_updates/Fig5/01_analyze_gfap_iba1_microglia.py`. The classifier-bundle serialization and loading code is in `figure_updates/Fig5/09_improve_microglia_classifier.py` and `figure_updates/Fig5/10_audit_microglia_transfer.py`.

When downloading manually or with GitHub CLI, restore the selected classifier's original filename by copying `morphology_forest.joblib` to `selected_classifier.joblib` in the same restored classifier directory. The Python download helper creates this alias automatically. Both names must have the SHA-256 value recorded for `morphology_forest.joblib`. No second download is needed.

## External pretrained dependency

The DINO and hybrid classifier candidates require the frozen external DINOv2 ViT-S/14 feature encoder and the feature-processing code. This encoder was **not trained by this project** and is not included among the release assets.

- Upstream file: [dinov2_vits14_pretrain.pth](https://dl.fbaipublicfiles.com/dinov2/dinov2_vits14/dinov2_vits14_pretrain.pth)
- Size: 88,283,115 bytes
- SHA-256: `b938bf1bc15cd2ec0feacfe3a1bb553fe8ea9ca46a7e1d8d00217f29aef60cd9`
- Original project path: `analyses/Fig5/model_assets/dinov2/dinov2_vits14_pretrain.pth`

Generic Cellpose and other externally pretrained caches are likewise excluded. Use the repository's documented software environments and model implementations when loading these artifacts.

## License and attribution

The repository's [license matrix](../LICENSE) governs the authors' contributions: software and supporting configuration use MIT, and the authors' original non-code materials use CC BY 4.0. File-specific notices and applicable third-party terms take precedence. This archive does not assign new rights to upstream pretrained components or relicense third-party material. Refer to the analysis-code citations and the original model providers when attributing reused software and pretrained components.

## Model inventory

| Release asset | Role | Bytes | Suggested restore path |
|---|---|---:|---|
| `cellpose__cfos_channel` | Earlier segmentation iteration | 609,384,247 | `models/cellpose/cfos_channel` |
| `cellpose__cfos_channel2` | Earlier segmentation iteration | 609,384,603 | `models/cellpose/cfos_channel2` |
| `cellpose__cfos_channel3` | Manuscript segmentation | 609,384,603 | `models/cellpose/cfos_channel3` |
| `cellpose__cpsam_cfos1` | Earlier segmentation iteration | 609,383,891 | `models/cellpose/cpsam_cfos1` |
| `cellpose__dapi_channel` | Manuscript segmentation | 609,384,247 | `models/cellpose/dapi_channel` |
| `cellpose__POMC_1` | Earlier segmentation iteration | 609,375,839 | `models/cellpose/POMC_1` |
| `cellpose__POMC_2` | Manuscript segmentation | 609,375,839 | `models/cellpose/POMC_2` |
| `cellpose__cfos_Keyence1` | Manuscript segmentation | 609,384,603 | `models/cellpose/cfos_Keyence1` |
| `cellpose__cpsam_20260107_110116` | Earlier segmentation iteration | 609,387,451 | `models/cellpose/cpsam_20260107_110116` |
| `cellpose__dapi_Keyence` | Earlier segmentation iteration | 609,384,247 | `models/cellpose/dapi_Keyence` |
| `cellpose__dapi_Keyence3` | Earlier segmentation iteration | 609,384,603 | `models/cellpose/dapi_Keyence3` |
| `cellpose__dapi_Keyence4` | Manuscript segmentation | 609,384,603 | `models/cellpose/dapi_Keyence4` |
| `cellpose__dapi_Keyence2` | Earlier segmentation iteration | 609,384,603 | `models/cellpose/dapi_Keyence2` |
| `treatment_tinymorphcnn.pt` | Manuscript treatment analysis | 331,757 | `models/tinymorphcnn/treatment_tinymorphcnn.pt` |
| `reviewed_initial_tinymorphcnn.pt` | Reviewed CNN candidate | 331,757 | `models/tinymorphcnn/reviewed_initial_tinymorphcnn.pt` |
| `reviewed_v2_tinymorphcnn.pt` | Reviewed CNN candidate | 331,885 | `models/tinymorphcnn/reviewed_v2_tinymorphcnn.pt` |
| `reviewed_v3_tinymorphcnn.pt` | Reviewed CNN candidate | 331,885 | `models/tinymorphcnn/reviewed_v3_tinymorphcnn.pt` |
| `dino_logistic.joblib` | Full-reviewed-set classifier candidate | 463,216 | `models/reviewed_classifiers/dino_logistic.joblib` |
| `hybrid_logistic.joblib` | Full-reviewed-set classifier candidate | 465,522 | `models/reviewed_classifiers/hybrid_logistic.joblib` |
| `morphology_forest.joblib` | Full-reviewed-set classifier candidate | 609,998 | `models/reviewed_classifiers/morphology_forest.joblib` |
| `morphology_logistic.joblib` | Full-reviewed-set classifier candidate | 15,901 | `models/reviewed_classifiers/morphology_logistic.joblib` |

The 13 supporting loss arrays and all individual SHA-256 values are listed in `manifest.json` and `SHA256SUMS.txt`.
