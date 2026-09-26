# Model implementation and computational resources

The accompanying JSON records inspected source-file hashes and distinguishes
implementation defaults from settings confirmed by historical run records.

| Component | Manuscript implementation |
|---|---|
| TinyMorphCNN | Study-specific; 79,780 trainable parameters; saved checkpoint loads into the implemented architecture |
| Layers | Four 3×3 convolutions, 16/32/64/96 channels; each has batch normalization and ReLU; first three have 2×2 max pooling; global average pooling; dropout 0.25; linear 96→4; softmax probabilities |
| Inputs | 96×96 pixels; Iba1 intensity, cell mask and skeleton/distance morphology |
| Released optimization defaults | Class-weighted cross-entropy; AdamW; learning rate 0.001, weight decay 0.0001, batch size 64, seed 20260620, patience six; no learning-rate scheduler |
| Confirmed Fig. 5 run | CPU; 17,112 pseudo-labeled candidate cells from 15 animals; whole-animal 80/20 split; 25-epoch cap, stopped at 11; minimum validation loss at epoch 5 |
| DINOv2 | Externally pretrained ViT-S/14 visual foundation model; frozen encoder; no microscopy fine-tuning |
| DINOv2 inputs/features | 224×224 target/context crops; pretrained normalization; four right-angle rotations; class and mask-weighted patch tokens; 1,536 features per cell |
| S6 downstream fitting | Training-only scaling and 32-component whitened PCA; class-balanced multinomial logistic regression; C=0.1/1/10; L-BFGS, max_iter=2000, tolerance=1e-6 |
| S6 validation | 15 leave-one-animal-out folds; three stratified group inner folds; seed 20260907 |
| Recorded DINOv2 execution | CUDA feature extraction, PyTorch 2.8.0+cu128 |
| Author-confirmed workstation | AMD EPYC 7713, 64 cores/128 threads; approximately 1 TiB RAM (995.51 GiB OS-visible); NVIDIA RTX A5000, 24 GB |

The Fig. 5 checkpoint lacks a complete historical optimizer configuration.
Defaults above come from released code; device and epoch history are recorded
independently. Later manually reviewed TinyMorphCNN development runs are
separate from the Fig. 5 pseudo-label model. S6 compares morphometric proposals
against DINOv2 plus a fitted logistic classifier; it does not directly validate
Fig. 5 TinyMorphCNN. A foundation model supplies externally learned general
representations. Freezing its encoder does not preclude fitting a downstream
classifier. Morphology classes do not independently establish inflammation.

Sources: [TinyMorphCNN](../01_analyze_gfap_iba1_microglia.py),
[DINOv2 workflow](../09_improve_microglia_classifier.py),
[model downloads](../../../models/README.md).
References: [DINOv2](https://arxiv.org/abs/2304.07193),
[PyTorch](https://papers.neurips.cc/paper_files/paper/2019/hash/bdbca288fee7f92f2bfa9f7012727740-Abstract.html),
[AdamW](https://arxiv.org/abs/1711.05101).
