import scanpy as sc
import pandas as pd
import anndata as ad

def create_anndata(matrix_file, barcodes_file, features_file):
    adata = sc.read_mtx(matrix_file).T
    adata.obs_names = pd.read_csv(barcodes_file, header=None)[0]
    adata.var_names = pd.read_csv(features_file, sep='\t', header=None)[0]
    return adata

files = [
    ('GSM6611295_P15306_5001_matrix.mtx', 'GSM6611295_P15306_5001_barcodes.tsv', 'GSM6611295_P15306_5001_features.tsv'),
    ('GSM6611296_P15306_5002_matrix.mtx', 'GSM6611296_P15306_5002_barcodes.tsv', 'GSM6611296_P15306_5002_features.tsv'),
    ('GSM6611297_P14601_4004_matrix.mtx', 'GSM6611297_P14601_4004_barcodes.tsv', 'GSM6611297_P14601_4004_features.tsv'),
    ('GSM6611298_P14601_4005_matrix.mtx', 'GSM6611298_P14601_4005_barcodes.tsv', 'GSM6611298_P14601_4005_features.tsv'),
    ('GSM6611299_P15306_5003_matrix.mtx', 'GSM6611299_P15306_5003_barcodes.tsv', 'GSM6611299_P15306_5003_features.tsv'),
    ('GSM6611300_P15306_5004_matrix.mtx', 'GSM6611300_P15306_5004_barcodes.tsv', 'GSM6611300_P15306_5004_features.tsv'),
    ('GSM6611301_10X_20_003_matrix.mtx', 'GSM6611301_10X_20_003_barcodes.tsv', 'GSM6611301_10X_20_003_features.tsv'),
    ('GSM6611302_10X_20_004_matrix.mtx', 'GSM6611302_10X_20_004_barcodes.tsv', 'GSM6611302_10X_20_004_features.tsv')]

adatas = [create_anndata(*file_pair) for file_pair in files]

adata_combined = ad.concat(adatas, axis=0, join='outer', label='batch', keys=['subject1_pre', 'subject1_post', 'subject2_pre', 'subject2_post', 'subject3_pre', 'subject3_post', 'subject3_DEEP_pre', 'subject3_DEEP_post'])

print(adata_combined)
