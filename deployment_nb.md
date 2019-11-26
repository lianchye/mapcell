### Deployment of siamese SNN for transferring annotation from one annotated experiment to another

Case Study 1: Annotating 10X data from manuscript
-------------------------------------------------

``` r
# Load Functions
source("./annotation_functions/function_lib_annotation.R")

# Read in External counts data matrix
ext_data.ct=readRDS("./external_data/tenX_data/ext_data_ct.rds")
ext_data.ct %>% dim() # Check dimensions
```

    ## [1] 33694  7865

``` r
# Load Trained Model to use and the reference data for transfer of annotation
gene_names    = "./models_reference/cell_NSCLC/P1_annot_reference/input_genes.csv"%>%  read_csv() %>% pull(symbol)
model_sel     = "./models_reference/cell_NSCLC/P1_annot_reference/weights.h5"     %>% build_embedding(rna_length = length(gene_names),weights_h5 = .)
refCells_sel  = "./models_reference/cell_NSCLC/P1_annot_reference/embed_ref.rds"  %>% readRDS()
metaCells_sel = "./models_reference/cell_NSCLC/P1_annot_reference/meta_ref.rds"   %>% readRDS()
```

``` r
# Peek at the trained model
model_sel
```

    ## Model
    ## ___________________________________________________________________________
    ## Layer (type)                     Output Shape                  Param #     
    ## ===========================================================================
    ## input_1 (InputLayer)             (None, 33694)                 0           
    ## ___________________________________________________________________________
    ## embedding (Sequential)           (None, 32)                    17530912    
    ## ===========================================================================
    ## Total params: 17,530,912
    ## Trainable params: 17,530,912
    ## Non-trainable params: 0
    ## ___________________________________________________________________________

``` r
# Perform Annotation
annotate_out.df = 
  annotate_fn_v2(
    ext_data.ct,                  # External 10X data
    model_sel,                    # Model trained on P1 data and human annotations
    refCells_sel,                 # Reference cells that were used for training
    metaCells_sel$Major_cell_type # Labels of the reference cells
  )
```

    ## embedding start

    ## annotating

``` r
# Output contains:
# distance to 21 closest cell [knn_dist] ;
# the label of the cells [knn_annot] ; 
# Majority agreement of cell type across the 21 cells [annot_call];
# Number of reference cell supporting the majority call [annot_cts]

annotate_out_perCell.df = annotate_out.df %>% group_by(cell_bc) %>% group_by(cell_bc) %>% arrange(-annot_cts) %>% slice(1) %>% ungroup()
annotate_out_perCell.df
```

    ## # A tibble: 7,865 x 5
    ##    cell_bc          knn_annot     knn_dist annot_call annot_cts
    ##    <chr>            <chr>            <dbl> <chr>          <int>
    ##  1 AAACCCAAGATTGTGA bMonocytes 0.000000561 bMonocytes        21
    ##  2 AAACCCACATCGGTTA bMonocytes 0.000166    bMonocytes        21
    ##  3 AAACCCAGTACCGCGT tMoMacDC   0.0123      tMoMacDC          21
    ##  4 AAACCCAGTATCGAAA tB cells   0.0163      bNK cells          7
    ##  5 AAACCCAGTCGTCATA bNK cells  0.000114    bNK cells         12
    ##  6 AAACCCAGTCTACACA tMoMacDC   0.00991     bMonocytes        20
    ##  7 AAACCCAGTGCAAGAC bMonocytes 0.0000169   bMonocytes        21
    ##  8 AAACCCAGTGCATTTG tMoMacDC   0.0225      bMonocytes        20
    ##  9 AAACCCATCCGATGTA bT cells   0.00607     bT cells          16
    ## 10 AAACCCATCTCAACGA bNK cells  0.00000552  bT cells          16
    ## # ... with 7,855 more rows
