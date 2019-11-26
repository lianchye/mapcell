# Required Libraries
library(class)
library(ggsci)
library(ggpubr)
library(keras)
library(tidyverse)
library(KernelKnn)
library(doParallel)

# Optional libraries for plotting
# library(viridis)
# library(viridisLite)
# library(Rtsne.multicore)
# library(ggalluvial)

# Annotation function [Vers 1]
annotate_fn = function(
  ext_data.ct,     # External Data for annotataion
  model_embed,     # Embedding Neural Network Keras Model
  model_embed_ref, # Reference Cells Count vectors 
  ref_annotation,  # Annotation/Cell type labels of reference cells
  neighbors=  21,  # Annotation parameter 1
  nos_threads=10,  # Multicore threads
  method=""){
  ext_data.embed = model_embed %>% predict(ext_data.ct %>% apply(., 2, function(x){x/max(x)}) %>% as.matrix() %>% t(),verbose=1)
  snn_pred       = knn(model_embed_ref,ext_data.embed,ref_annotation,k=neighbors,prob=FALSE) 
  knn_index_out  = KernelKnn::knn.index.dist(model_embed_ref,ext_data.embed,k=neighbors,method = "pearson_correlation",threads = nos_threads)
  snn_annot.df   = data.frame(cell_bc     = colnames(ext_data.ct),
                              snn_pred    = snn_pred,
                              knn1_dist   = knn_index_out$test_knn_dist[,1],
                              knnAll_dist = knn_index_out$test_knn_dist[,1:neighbors] %>% rowSums(),
                              stringsAsFactors = FALSE)
  snn_annot.df
}

# Annotation function [Vers 2]
annotate_fn_v2 = function(
  ext_data.ct,     # External Data for annotataion
  model_embed,     # Embedding Neural Network Keras Model
  model_embed_ref, # Reference Cells Count vectors
  ref_annotation,  # Annotation/Cell type labels of reference cells
  neighbors   = 21,
  method      = "pearson_correlation",
  nos_threads = 10 # Multicore threads
){
  message("embedding start")
  ext_data.embed = model_embed %>% predict(ext_data.ct %>% apply(., 2, function(x){x/max(x)}) %>% as.matrix() %>% t(),verbose=1)
  knn_index_out  = KernelKnn::knn.index.dist(model_embed_ref,ext_data.embed,k=neighbors,method = method,threads = nos_threads)
  
  message("annotating")
  snn_annot.df   = lapply(1:nrow(ext_data.embed),function(x){
    data.frame(cell_bc     = colnames(ext_data.ct)[x],
               knn_annot   = knn_index_out$test_knn_idx[x,] %>% ref_annotation[.],
               knn_dist    = knn_index_out$test_knn_dist[x,],
               stringsAsFactors = FALSE)
    
  }) %>% do.call("rbind",.) %>%
    group_by(cell_bc) %>%
    mutate(annot_call= table(knn_annot) %>% sort(decreasing = TRUE) %>% .[1] %>% names(),
           annot_cts= table(knn_annot) %>% sort(decreasing = TRUE) %>% .[1] )
  snn_annot.df
}


# Sparse PCA
sparse_pca=function(x,n_pcs=30){
  svd_res=irlba::irlba(x,n_pcs)
  n <- dim(x)[1]
  s <- apply(x, 2, sd, na.rm=TRUE)
  variance_sum <- sum(apply(x,2,var,na.rm=TRUE)/(s^2)) # sample variance sum
  var_pcs <- svd_res$d^2/(n-1)/variance_sum
  return(list(x=svd_res$u %*% diag(svd_res$d),tot_var=variance_sum))
}

# Build embedding layer
build_embedding  = function(layer1_nodes = 512, 
                            layer2_nodes = 512,
                            layer3_nodes = 32,
                            rna_length   = 33694,
                            weights_h5   = "na"){
  
  original_dim <- rna_length %>% as.integer()
  input_a      <- layer_input(shape = c(original_dim),name="input_1")
  
  # Embedding Neural Network with dropout layer to prevent overfitting
  embedding_nn <-
    keras_model_sequential(name = "embedding") %>%
    layer_dropout (rate = 0.5) %>%
    layer_dense   (layer1_nodes, kernel_regularizer =regularizer_l1(l = 1e-4)) %>%
    layer_dropout (rate = 0.5) %>%
    layer_dense   (layer2_nodes, activation = "relu") %>%
    layer_dropout (rate = 0.5) %>%
    layer_dense   (layer3_nodes , activation = "relu")
  
  out.mod = keras_model(inputs = input_a, input_a %>% embedding_nn)
  
  if(!(weights_h5=="na")){
    out.mod %>% load_model_weights_hdf5(weights_h5)
  }
  return(out.mod)
}

