import scanpy as sc
import pandas as pd
import anndata as ad

# create anndata for subject 1 (pre)
adata_sample1_pre = sc.read_mtx('GSM6611295_P15306_5001_matrix.mtx').T
adata_sample1_pre.obs_names = pd.read_csv('GSM6611295_P15306_5001_barcodes.tsv', header=None)[0]
adata_sample1_pre.var_names = pd.read_csv('GSM6611295_P15306_5001_features.tsv', sep='\t', header=None)[0]
print(adata_sample1_pre) #should read AnnData object with n_obs × n_vars = 2462 × 33538

# create anndata for subject 1 (post)
adata_sample1_post = sc.read_mtx('GSM6611296_P15306_5002_matrix.mtx').T
adata_sample1_post.obs_names = pd.read_csv('GSM6611296_P15306_5002_barcodes.tsv', header=None)[0]
adata_sample1_post.var_names = pd.read_csv('GSM6611296_P15306_5002_features.tsv', sep='\t', header=None)[0]
print(adata_sample1_post) #should read AnnData object with n_obs × n_vars = 3070 × 33538

# create anndata for subject 2 (pre)
adata_sample2_pre = sc.read_mtx('GSM6611297_P14601_4004_matrix.mtx').T
adata_sample2_pre.obs_names = pd.read_csv('GSM6611297_P14601_4004_barcodes.tsv', header=None)[0]
adata_sample2_pre.var_names = pd.read_csv('GSM6611297_P14601_4004_features.tsv', sep='\t', header=None)[0]
print(adata_sample2_pre) #output = AnnData object with n_obs × n_vars = 4777 × 33538

# create anndata for subject 2 (post)
adata_sample2_post = sc.read_mtx('GSM6611298_P14601_4005_matrix.mtx').T
adata_sample2_post.obs_names = pd.read_csv('GSM6611298_P14601_4005_barcodes.tsv', header=None)[0]
adata_sample2_post.var_names = pd.read_csv('GSM6611298_P14601_4005_features.tsv', sep='\t', header=None)[0]
print(adata_sample2_post) # output = AnnData object with n_obs × n_vars = 1443 × 33538

# create anndata for subject 3 (pre)
adata_sample3_pre = sc.read_mtx('GSM6611299_P15306_5003_matrix.mtx').T
adata_sample3_pre.obs_names = pd.read_csv('GSM6611299_P15306_5003_barcodes.tsv', header=None)[0]
adata_sample3_pre.var_names = pd.read_csv('GSM6611299_P15306_5003_features.tsv', sep='\t', header=None)[0]
print(adata_sample3_pre) # output = AnnData object with n_obs × n_vars = 15668 × 33538

# create anndata for subject 3 (post)
adata_sample3_post = sc.read_mtx('GSM6611300_P15306_5004_matrix.mtx').T
adata_sample3_post.obs_names = pd.read_csv('GSM6611300_P15306_5004_barcodes.tsv', header=None)[0]
adata_sample3_post.var_names = pd.read_csv('GSM6611300_P15306_5004_features.tsv', sep='\t', header=None)[0]
print(adata_sample3_post) # output AnnData object with n_obs × n_vars = 9588 × 33538

# create anndata for subject 3 DEEP (pre)
adata_sample3_DEEP_pre = sc.read_mtx('GSM6611301_10X_20_003_matrix.mtx').T
adata_sample3_DEEP_pre.obs_names = pd.read_csv('GSM6611301_10X_20_003_barcodes.tsv', header=None)[0]
adata_sample3_DEEP_pre.var_names = pd.read_csv('GSM6611301_10X_20_003_features.tsv', sep='\t', header=None)[0]
print(adata_sample3_DEEP_pre) # AnnData object with n_obs × n_vars = 28474 × 33538

# create anndata for subject 3 DEEP (post)
adata_sample3_DEEP_post = sc.read_mtx('GSM6611302_10X_20_004_matrix.mtx').T
adata_sample3_DEEP_post.obs_names = pd.read_csv('GSM6611302_10X_20_004_barcodes.tsv', header=None)[0]
adata_sample3_DEEP_post.var_names = pd.read_csv('GSM6611302_10X_20_004_features.tsv', sep='\t', header=None)[0]
print(adata_sample3_DEEP_post) # AnnData object with n_obs × n_vars = 15564 × 33538
