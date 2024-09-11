# Document-clustering-techniques-comparative-analysis-and-development-of-a-novel-approach
Code linked to my master thesis

TThe idea for this project, as well as the initial part of the code, is inspired by the work presented in https://link.springer.com/chapter/10.1007/978-3-031-15509-3_33. The approach is straightforward: given a corpus, we first compute the p-value associated with the homogeneity test described in the referenced study. Specifically, after applying Latent Dirichlet Allocation (LDA), we measure the extent to which two documents share similar topics using the Kullback-Leibler divergence. 

Subsequently, a test-based clustering method is employed. In this method, the two documents with the highest p-value are merged, and then the p-value between this newly formed document and the remaining documents is computed. This process is repeated iteratively until no p-value exceeds the chosen significance level,  $\alpha$.
