# MapCell

<img src="./mapcell_logo.PNG" alt="logo" title="mapcell" style="zoom:150%;" />

#  <a href="https://f1000.com/prime/736854052?bd=1" target="_blank"><img src="https://s3.amazonaws.com/cdn.f1000.com/images/badges/badgef1000.gif" alt="Access the recommendation on F1000Prime" id="bg" /></a>  ![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)

MapCell: Learning a comparative cell type distance metric with *Siamese neural nets* with applications towards cell-type identification across experimental datasets [(biorxiv)](https://www.biorxiv.org/content/10.1101/828699v1) [(pdf)](https://www.biorxiv.org/content/biorxiv/early/2019/11/04/828699.full.pdf)

## Features  

- Requires few training examples per cell type (20 cells per cell type was used in this work) 

- Accurately map cells across different scRNA-seq platforms at single-cell level (rather than cluster level)

- The ability to identify novel cell types not seen in our training set and we demonstrated a process of retraining to incorporate new cell types into the model

- Generalizable distance metric to map unseen cell types without re-training

- Ability to map cell annotations across different species

- Highly scalable from small (training from cells of a single patient can be used to predict across 6 other patients) to large models (training from the Human Cell Landscape resource that comprises a wide survey of cell types derived from 50 different tissues.) 

## Usage

This repository is meant to host the trained models from the above paper and to serve as a easy launching point for you to start using these models and begin annotating your own scRNA-seq dataset.  

R notebooks examples are provided in the repository as to illustrate the ease of modding and using the trained models. Current examples includes:

1. [10X count peripheral blood annotation](./deployment_nb.md) and;
2. [3K PBMC seurat object annotation](./seurat_annotation.md).

To test out the annotation speed, we have also included a [notebook](./annotation_speed_stress_test) that records the annotation speed.

#### Environment Notes

Main requirement for running the R notebook examples is installation of [keras]( https://github.com/rstudio/keras) library for R. The given examples have been tested to run on a 16 core linux server with 32gb ram where 10 000 cells are annotated in 20 seconds. Higher RAM and better hardware such as availability of GPU will result in better performance. However, the models are intentionally made small so as to be able to run fast without GPUs.

## Citations

Please use the following bibtex entry:

```
@article {Koh828699,
	author = {Koh, Winston and Hoon, Shawn},
	title = {MapCell: Learning a comparative cell type distance metric with Siamese  neural nets with applications towards cell-types identification across experimental datasets},
	year = {2019},
	doi = {10.1101/828699},
	journal = {bioRxiv}
}
```

## Future Work

We are working on the release of a series of models featuring major tissue types pending the feedback and demand we received from the community.

Docker images to abstract away the need for notebooks are being built. 







