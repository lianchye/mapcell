# Deployment of siamese SNN for transferring annotation from an annotated experiment to a non-annotated one

Case Study 1: Annotating 10X data from manuscript
-------------------------------------------------

### Source function and load the unannotated external data

``` r
# Load Functions
source("./annotation_functions/function_lib_annotation.R")

# Read in External counts data matrix
ext_data.ct=readRDS("./external_data/tenX_data/ext_data_ct.rds")
ext_data.ct %>% dim() # Check dimensions
```

    ## [1] 33694  7865

### Load up the trained models 

``` r
# Load Trained Model to use and the reference data for transfer of annotation
gene_names    = "./models_reference/cell_NSCLC/P1_annot_reference/input_genes.csv"%>%  read_csv() %>% pull(symbol)
model_sel     = "./models_reference/cell_NSCLC/P1_annot_reference/weights.h5"     %>% build_embedding(rna_length = length(gene_names),weights_h5 = .)
refCells_sel  = "./models_reference/cell_NSCLC/P1_annot_reference/embed_ref.rds"  %>% readRDS()
metaCells_sel = "./models_reference/cell_NSCLC/P1_annot_reference/meta_ref.rds"   %>% readRDS()
```

### Check model structure

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

### Run annotation

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

### Understanding the output

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

### Read in meta-data to evaluate the results of the annotation using mapcell

```{r}
# Read in external meta data about cells
tsne_coord= readRDS("./external_data/tenX_data/ext10x_30_tsne.rds")
ext10_pro = readRDS("./external_data/tenX_data/pbmc10k_protein_v3_adt.rds")
plot.df = tsne_coord %>% left_join(annotate_out_perCell.df)
```

### Unsupervised clustering with no annotation 

```{r}
# plotting just the plane t-SNE with no annotation
before.plot = 
  plot.df %>% 
  ggplot(aes(x=X1,y=X2))+
  geom_point(size=0.7)+
  theme_pubclean()+
  guides(color = guide_legend(override.aes = list(size=5,alpha=1)))
before.plot %>% ggsave(plot = .,filename = "./plots/ext10_tsne.png")
before.plot
```

![ext10_tsne](.\plots\ext10_tsne.png)

### Now layer on the transferred annotation using SNN model as color

```{r}
# Plot the annotated data as colors
after.plot =plot.df %>% 
  ggplot(aes(x=X1,y=X2,color=annot_call))+
  geom_point(size=0.4,alpha=0.7)+
  theme_pubclean()+
  guides(color = guide_legend(override.aes = list(size=5,alpha=1)))
after.plot %>% ggsave(plot = .,filename = "./plots/ext10_annot.png")
after.plot
```

![ext10_annot](.\plots\ext10_annot.png)

### Verify the transferred annotation using independent protein meta-data

```{r}
# Check the annotation of the 
poi=rownames(ext10_pro)

plot_pro.df=plot.df %>% left_join(data.frame(cell_bc=colnames(ext10_pro)[-1],ext10_pro[,-1] %>% apply(.,2,function(x){log2((x/sum(x)*10^6)+1)}) %>% t()))

plot.list=lapply(poi,function(x){
  minor_pro.plot=
  plot_pro.df %>%  
  ggplot(aes_string(x="X1",y="X2",color=x))+
  geom_point(size=0.4)+
  scale_color_viridis_c(name=x,option="A")+
  coord_equal()+
  theme_pubclean()
  return(minor_pro.plot)
})
names(plot.list)=poi

plot_sub.list = plot.list[c("CD3","CD14","CD19","CD56")]
protein.plot=ggarrange(plotlist = plot_sub.list,ncol=length(plot_sub.list),nrow=1)
ggsave(protein.plot,file="./plots/ext10_protein_verify.png",width = 12,height = 4)
protein.plot
```

![ext10_protein_verify](C:\Users\Administrator\Documents\GitHub\mapcell\plots\ext10_protein_verify.png)