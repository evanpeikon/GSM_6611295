import scanpy as sc

# Step 1: Filter out cells with too few or too many expressed genes
adata_combined = adata_combined[adata_combined.obs['n_genes_by_counts'] > 200, :]
adata_combined = adata_combined[adata_combined.obs['n_genes_by_counts'] < 3000, :]

# Step 2: Filter out genes expressed in fewer than or equal to 20 cells
sc.pp.filter_genes(adata_combined, min_cells=21)

# Step 3: Calculate the percentage of mitochondrial genes
adata_combined.var['mt'] = adata_combined.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_combined, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Step 4: Filter out cells with high mitochondrial gene content
adata_combined = adata_combined[adata_combined.obs['pct_counts_mt'] < 14, :]

print(adata_combined)
