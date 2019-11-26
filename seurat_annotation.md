### Deployment of siamese SNN for transferring annotation from one annotated experiment to another

Case Study 1: Annotating 10X data from manuscript
-------------------------------------------------

``` r
# Load Functions
source("./annotation_functions/function_lib_annotation.R")

# Read in External counts data matrix
library(Seurat)
seurat_obj   = readRDS("./external_data/seurat_objs/pbmc3k_final.rds")
```

``` r
# Load Trained Model to use and the reference data for transfer of annotation
gene_names    = "./models_reference/human_cell_landscape/peri-blood/gene_names.rds"%>%  readRDS() 
model_sel     = "./models_reference/human_cell_landscape/peri-blood/weights.h5"     %>% build_embedding(rna_length = length(gene_names),weights_h5 = .)
refCells_sel  = "./models_reference/human_cell_landscape/peri-blood/embed_ref.rds"  %>% readRDS()
metaCells_sel = "./models_reference/human_cell_landscape/peri-blood/meta_ref.rds"   %>% readRDS()

common_genes  = seurat_obj@raw.data@Dimnames[[1]] %>% intersect(gene_names)
```

``` r
sel_cols = 1:ncol(seurat_obj@raw.data)
chunks   = split(sel_cols, ceiling(seq_along(sel_cols)/4000)) 

annot_out=list()
for(i in 1:length(chunks)){
    
  common.mat      =
    apply(seurat_obj@raw.data[,chunks[[i]]],2,function(x){
      api.vec               = rep(0L,length(gene_names))
      names(api.vec)        = gene_names
      api.vec[common_genes] = x[common_genes]
      cell_exp              = as.numeric(api.vec/max(api.vec))
      return(cell_exp)
    })
  rownames(common.mat) = gene_names
  
  annotate_cell.list = annotate_fn_v2(common.mat, 
                                      model_sel,
                                      refCells_sel,
                                      metaCells_sel$Celltype_clean,
                                      neighbors = 20)
  
annot_out[[i]] = 
    annotate_cell.list %>%
    arrange(-annot_cts) %>% slice(1) %>% ungroup() %>% arrange(annot_cts)

}
```

    ## embedding start

    ## annotating

``` r
annot_all = 
  annot_out %>% do.call("rbind",.) %>%
  group_by(cell_bc) %>% top_n(n=1,wt=-knn_dist)
```

``` r
combined_annot.df = seurat_obj@meta.data %>% rownames_to_column("cell_bc") %>% left_join(annot_all)
```

    ## Joining, by = "cell_bc"
