![logo](./mapcell_logo.PNG "mapcell")
# mapcell
MapCell: Learning a comparative cell type distance metric with Siamese neural nets with applications towards cell-type identification across experimental datasets [(biorxiv)](https://www.biorxiv.org/content/10.1101/828699v1) [(pdf)](https://www.biorxiv.org/content/biorxiv/early/2019/11/04/828699.full.pdf)

Features of MapCell are: 
-	Requires few training examples per cell type (20 cells per cell type was used in this work) 
-	Accurately map cells across different scRNA-seq platforms at single-cell level (rather than cluster level)
-	The ability to identify novel cell types not seen in our training set and we demonstrated a process of retraining to incorporate new cell types into the model
-	Generalizable distance metric to map unseen cell types without re-training
-	Ability to map cell annotations across different species
-	Highly scalable from small (training from cells of a single patient can be used to predict across 6 other patients) to large models (training from the Human Cell Landscape resource that comprises a wide survey of cell types derived from 50 different tissues.) 
