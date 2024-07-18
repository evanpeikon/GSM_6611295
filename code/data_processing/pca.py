import scanpy as sc

# perform z-transformation
sc.pp.scale(adata_combined, zero_center=True)

# perform PCA and plot elbow plot
sc.tl.pca(adata_combined, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_combined, log=True) 

# perform PCA, retaining 7 components
sc.tl.pca(adata_combined, svd_solver='arpack', n_comps=7)
