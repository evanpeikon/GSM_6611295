n_pcs_available = adata_combined.obsm['X_pca'].shape[1]
sc.pp.neighbors(adata_combined, n_pcs=n_pcs_available, use_rep='X_pca', metric='cosine', random_state=0)
sc.tl.leiden(adata_combined, resolution=1.0, flavor="igraph", n_iterations=2, directed=False, random_state=0)
sc.tl.umap(adata_combined)
sc.pl.umap(adata_combined, color=['batch'], ncols=1, palette='Set1')
