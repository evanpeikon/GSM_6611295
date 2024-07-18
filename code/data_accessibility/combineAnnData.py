# concatenate the 8 AnnData objects along the observations axis (cells) and get summary stats
adata_combined = ad.concat([adata_sample1_pre, adata_sample1_post,adata_sample2_pre, adata_sample2_post, adata_sample3_pre, adata_sample3_post, adata_sample3_DEEP_pre, adata_sample3_DEEP_post], axis=0, join='outer', label='batch', keys=['subject1_pre', 'subject1_post','subject2_pre', 'subject2_post', 'subject3_pre', 'subject3_post', 'subject3_DEEP_pre', 'subject3_DEEP_post'])

num_cells = adata_combined.n_obs
print(f"Number of cells: {num_cells}") # print number of cells

num_genes = adata_combined.n_vars
print(f"Number of genes: {num_genes}") # print number of genes

num_batches = adata_combined.obs['batch'].nunique()
print(f"Number of libraries (batches): {num_batches}") # print of libraries (batches)

batch_distribution = adata_combined.obs['batch'].value_counts()
print("Batch distribution:")
print(batch_distribution) # insect the batch (library) distribution
