import scanpy as sc

# apply global scaling normalization
sc.pp.normalize_total(adata_combined, target_sum=1e4)  # scaling factor is 10,000 (1e4)
sc.pp.log1p(adata_combined)  # log transformation

# find the 2000 most highly variable genes
sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000, subset=True)

print(adata_combined) 
