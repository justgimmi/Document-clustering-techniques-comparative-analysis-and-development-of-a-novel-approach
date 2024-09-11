setwd("C:/Users/apalm/Documents/GitHub/Document-clustering-techniques-comparative-analysis-and-development-of-a-novel-approach")

##### Package used during the thesis ####
library(future.apply);library(tm);library(topicmodels)
library(parallel);library(doParallel);library(quanteda);library(seededlda)
library(foreach)
library(dbscan)
library(mclust)
library(kernlab)
library(cluster)
library(proxy)
library(lsa)
library(datasets)
library(caret)
library(OneR)
library(pROC)
library(factoextra)
library(clusterSim)










##### P_Value estimation via Bootstrap ####

# The code is a slightly modification of the original paper in order to take into
# account computational complexity (it is parallelized)

# Original_Paper: "https://link.springer.com/chapter/10.1007/978-3-031-15509-3_33"

replica_function = function(dtm, term_doc12_vector, p, ndoc1, ndoc2, k, index1, index2) {
  #' @param dtm description: Document-Term Matrix
  #' @param term_doc_12_vector description: Vector of size 1xp obtained multiplying the corresponding row of theta and beta
  #' @param p description: Dimension of the vocabulary
  #' @param ndoc1 description: Number of words in the first document
  #' @param ndoc2 description: Number of words in the second document
  #' @param k description: Number of topics in LDA
  #' @param index1 description: Index of the first document
  #' @param index2 description: Index of the second document
  # Resample the two documents according to the new estimated probabilities
  
  new_terms_doc1 = data.frame(table(sample(1:p, size = ndoc1, replace = TRUE,
                                           prob = term_doc12_vector)))
  new_terms_doc1$Var1 = as.numeric(as.character(new_terms_doc1$Var1))
  new_terms_doc2 = data.frame(table(sample(1:p, size = ndoc2, replace = TRUE,
                                           prob = term_doc12_vector)))
  new_terms_doc2$Var1 = as.numeric(as.character(new_terms_doc2$Var1))
  newdoc1_dtm = numeric(p)
  for (j in 1:nrow(new_terms_doc1)) {newdoc1_dtm[new_terms_doc1[j,1]] = new_terms_doc1[j,2]}
  newdoc2_dtm = numeric(p)
  for (j in 1:nrow(new_terms_doc2)) {newdoc2_dtm[new_terms_doc2[j,1]] = new_terms_doc2[j,2]}
  
  # Substitute the two documents inside the original matrix
  dtm_new = dtm
  dtm_new[index1,] = newdoc1_dtm
  dtm_new[index2,] = newdoc2_dtm
  dtm_new = as.dfm(dtm)
  # Fit again LDA
  z_new = LDA(dtm_new, k = k, method = "Gibbs", control = list(alpha = 1/k, iter = 500, burnin = 100))
  theta_new = z_new@gamma    
  phi_new = exp(z_new@beta)
  # Compute the bootstrap replicate value of the Kullback-Leibler divergence
  sim_doc1_new = data.frame(theta_new[index1,])
  sim_doc2_new = data.frame(theta_new[index2,])
  KLD_star = sum(sim_doc1_new*log(sim_doc1_new/sim_doc2_new))
  
  return(KLD_star)
  
}

boottest_similarity<-function(dtm,index1,index2,k, B){
  #' @param dtm description: Document-Term matrix
  #' @param index1 description: Index of the first document
  #' @param index2 description: Index of the second document
  #' @param k description: Number of topics in LDA
  #' @param B description: Number of bootstrap replicates

  #'The following function is pretty smooth. Given two indices, it computes
  #'the bootstrap p-value linked with an homogeneity test.
  #' 
  if (index1 >= index2) {
    print("Not valid values")
    return(NaN)
    
  }
  print("Iniziamo")
  ndoc1<-sum(dtm[index1,]) # length of the first document
  ndoc2<-sum(dtm[index2,]) # length of the second document 
  p<-ncol(dtm) # number of words in the vocabulary

  z = LDA(dtm, k = k, method = "Gibbs",
          control = list(alpha = 1/k, iter = 500, burnin = 100))

  theta = z@gamma 
  
  doc1 = data.frame(theta[index1,])
  doc2 = data.frame(theta[index2,])
  
  KLD = sum(doc1*log(doc1/doc2))

  #Merge the two documents of interest
  
  dtm_merged = dtm
  
  dtm_doc12 = apply(dtm_merged[c(index1,index2), ], 2, sum)
  
  #Replace the first doc with the merged document
  dtm_merged[index1, ] = dtm_doc12

  dtm_merged <- dtm_merged[-c(index2), ]
  dtm_merged = as.dfm(dtm_merged)

  
  print("merging the rows")
  #Run LDA with merged abstract to create population
  
  z_merged = LDA(dtm_merged, k = k, method = "Gibbs", control = list(alpha = 1/k, iter = 500, burnin = 100))

  theta_merged = z_merged@gamma    
  phi_merged = exp(z_merged@beta)
  #Compute term probability of terms in merged document
  theta_doc12 = data.frame(t(theta_merged[index1,]))
  #term_doc12_vector<-as.matrix(theta_doc12) %*% as.matrix(t(phi_merged))
  term_doc12_vector = as.matrix(theta_doc12) %*% as.matrix(phi_merged)

  
  print("starting the bootstrap loop")
  
  
  results = future_sapply(seq_len(B), function(i, dtm, term_doc12_vector, p, ndoc1, ndoc2, k, index1, index2) replica_function(dtm, term_doc12_vector, p, ndoc1, ndoc2, k, index1, index2),
                          dtm = dtm, term_doc12_vector = term_doc12_vector,
                          p = p, ndoc1 = ndoc1, ndoc2 = ndoc2, k = k, index1 = index1, index2 = index2,
                          future.seed = T)
  
  
  print("Done this cycle")
  result = list(pvalue = mean(results >= KLD), stats = results, KLD = KLD)
  return(result)
}

##### P_Value_Matrix ####
# P_value Matrix
p_value_matrix = function(ind, m, dtm, k, B, theta){
  #' @param ind description: indices to consider in the document-term matrix
  #' @param m description: empty matrix of size M x M
  #' @param k description: number of topics for LDA
  #' @param B description: number of Bootstrap replicates
  #' @param theta description: document topic matrix obtained fitting LDA
  
  num_cores = detectCores(logical = T)-1
  plan(multisession, workers = num_cores)
  message("Number of parallel workers: ", nbrOfWorkers())
  # Start the parallelized session
  iter = 1
  for (i in seq_len(nrow(ind))) {
    cat("We are at the iteration number", iter, "over", nrow(ind), "\n")
    m[ind[i,1], ind[i,2]] = boottest_similarity(dtm = dtm, index1 = ind[i, 1], index2 = ind[i,2], k = k,  B = B, theta = theta)$pvalue
    # compute the Bootstrap pvalue
    m[ind[i, 2], ind[i,1]] = m[ind[i,1], ind[i,2]] # assign the computed value to the element in  the lower triangular matrix
    iter = iter + 1
  }
  return(m)
}


#### Test_Based_Clustering_Method ####

clustering_method = function(mat, dtm, alpha, B, k){
  #' @param mat description: starting p-values matrix
  #' @param dtm description: initially document term matrix
  #' @param alpha description: significance level
  #' @param B description: Bootstrap replicates
  #' @param K description: number of topics
  
  M = ncol(mat)
  partizione = lapply(1:M, function(x) list(x)) # save the partition for each iteration
  colnames(mat) = paste("Cluster", sep = " ", seq_len(M))
  storage = matrix(0, nrow = 1, ncol = M) # store clusters
  colnames(storage) = colnames(mat)
  p_values = numeric(M) # save merging height
  au = colnames(mat)
  cluster = 1
  iter = 1
  partizione[[iter]] = storage  # save the trivial partition
  
  while (any(mat[upper.tri(mat) == T] >= alpha) &  (iter < M-1)) {
    # while there exists a p_value greater than alpha
    iter = iter + 1
    cat("Iteration number", iter, "/", M, "\n")
    inde = which(upper.tri(mat) == T & mat == max(mat), arr.ind = T)[1,] # find the documents corresponding to the highest pvalue
    p_values[iter] = max(mat) # save this value
    
    # Storage procedure:
    # - find the names of the clusters;
    # - merge them in the storage value according to the rule below
    r = au[inde[1]]
    r2 = au[inde[2]]
    if (storage[1, r] == 0 & storage[1, r2] == 0){
      storage[1, r] = cluster
      storage[1, r2] = cluster
      cluster = cluster + 1}
    else if (storage[1, r] > 0 & storage[1, r2] == 0) {
      storage[1, r2] = storage[1, r]}
    else if (storage[1, r] ==  0 & storage[1, r2] > 0) {
      storage[1, r] = storage[1, r2] }   
    else if (storage[1, r] > 0 & storage[1, r2] > 0) {
      storage[storage == storage[1, r2]] = storage[1, r] }
    
    partizione[[iter]] = storage # save the obtained partition
    
    if (iter == M-1) {break}
    
    dtm[inde[1], ] = apply(dtm[c(inde[1], inde[2]), ], 2, sum) # merge the two documents in the Document Term matrix
    
    dtm = dtm[-c(inde[2]), ] # delete the second document from the Document Term matrix
    
    au = au[au != r2] # delete from the list of all the clusters the one that has been merged
    
    dtm = as.dfm(dtm)
    mat = mat[-c(inde[2]), -c(inde[2])] #delete from the matrix of pvalue the row and the column of the cluster that has been merged
    # We need to compute the pvalues between the new document and the remaining ones
    
    s = seq(from = 1, to = ncol(mat)) # define the sequence of all the remaining clusters excluding the newly formed
    
    s = s[s != inde[1]]
    
    index = data.frame(rep(inde[1], length(s)), s)
    # This part of code is needed due to the fact the function pvalue matrix relies on another function that work exclusively when index1 > index2
    for (i in 1:nrow(index)) {
      if (index[i, 1] > index[i, 2]) {
        aux = index[i, 1]
        index[i, 1] = index[i, 2]
        index[i, 2] = aux
      }      
    }
    z = LDA(dtm, k = k, method = "Gibbs",
            control = list(alpha = 1/k, iter = 1000)) # fit LDA with the new Document Term matrix
    theta = z@gamma
    mat = p_value_matrix(ind = index, m = mat, dtm = dtm, k = k, B = B, theta = theta) # compute the pvalue between the newly formed cluster and the others
  }
  plan(sequential) # stop the parallelized procedure
  print("Resting my workers")
  results = list(clustering = storage, height = p_values, partition = partizione)
  return(results)
}