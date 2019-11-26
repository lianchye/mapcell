Deployment of siamese SNN for transferring annotation from one annotated experiment to another
==============================================================================================

Load 10X reference 10k pbmc data
--------------------------------

``` r
# Load Functions
source("./annotation_functions/function_lib_annotation.R")

# Read in External counts data matrix
ext_data.ct=readRDS("./external_data/tenX_data/ext_data_ct.rds")
ext_data.ct %>% dim() # Check dimensions
```

    ## [1] 33694  7865

load model reference
--------------------

``` r
# Load Trained Model to use and the reference data for transfer of annotation
gene_names    = "./models_reference/cell_NSCLC/P1_annot_reference/input_genes.csv" %>%  read_csv() %>% pull(symbol)
model_sel     = "./models_reference/cell_NSCLC/P1_annot_reference/weights.h5"      %>% build_embedding(rna_length = length(gene_names),weights_h5 = .)
refCells_sel  = "./models_reference/cell_NSCLC/P1_annot_reference/embed_ref.rds"   %>% readRDS()
metaCells_sel = "./models_reference/cell_NSCLC/P1_annot_reference/meta_ref.rds"    %>% readRDS()
```

Verify model structure
----------------------

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

Sample across different cell number and collect the annotation time to run the model on them
--------------------------------------------------------------------------------------------

``` r
# Perform Annotation
cell_count_range =c(10,50,100,500,1000,2500,5000,7500,10000)
time.out =list()

for(i in 1L:length(cell_count_range)){ # i=1
  sample_cells = sample(1:ncol(ext_data.ct),cell_count_range[i],replace=TRUE)
  
  time.out[[i]]=
    system.time({
      annotate_out.df = 
        annotate_fn_v2(
          ext_data.ct[,sample_cells],   # External 10X data
          model_sel,                    # Model trained on P1 data and human annotations
          refCells_sel,                 # Reference cells that were used for training
          metaCells_sel$Major_cell_type # Labels of the reference cells
        )
    })
  
  gc()
}
```

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

    ## embedding start

    ## annotating

Plot the time changes with number of cells
------------------------------------------

``` r
speed.plot=
  data.frame(nos_cells_annot = cell_count_range,
             elapsed_time    = lapply(time.out,function(x){x[3]}) %>% unlist() %>% as.numeric()) %>% 
  ggplot(aes(x=nos_cells_annot,y=elapsed_time))+
  geom_point()+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  # scale_y_log10(
  #  breaks = scales::trans_breaks("log10", function(x) 10^x),
  #  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  geom_line(alpha=0.7)+
  ylab("Elapsed Time (seconds)")+
  xlab("Number of cell annotated in a batch on 32gb ram") +
  theme_pubclean() 

speed.plot
```

![](annotation_speed_stress_test_files/figure-markdown_github/unnamed-chunk-4-1.png)
