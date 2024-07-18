# Make a new project directory, then enter new directory
mkdir lovric_2022,
cd mkdir lovric_2022

# get sample 1 data (pre) then decompress files
wget -O GSM6611295_P15306_5001_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611295&format=file&file=GSM6611295%5FP15306%5F5001%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611295_P15306_5001_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611295&format=file&file=GSM6611295%5FP15306%5F5001%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611295_P15306_5001_matrix.mtx.gz	 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611295&format=file&file=GSM6611295%5FP15306%5F5001%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611295_P15306_5001_barcodes.tsv.gz
gunzip GSM6611295_P15306_5001_features.tsv.gz
gunzip GSM6611295_P15306_5001_matrix.mtx.gz

# get sample 1 data (post) then decompress files
wget -O GSM6611296_P15306_5002_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611296&format=file&file=GSM6611296%5FP15306%5F5002%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611296_P15306_5002_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611296&format=file&file=GSM6611296%5FP15306%5F5002%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611296_P15306_5002_matrix.mtx.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611296&format=file&file=GSM6611296%5FP15306%5F5002%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611296_P15306_5002_barcodes.tsv.gz
gunzip GSM6611296_P15306_5002_features.tsv.gz
gunzip GSM6611296_P15306_5002_matrix.mtx.gz

# get sample 2 data (pre) then decompress files
wget -O GSM6611297_P14601_4004_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611297&format=file&file=GSM6611297%5FP14601%5F4004%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611297_P14601_4004_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611297&format=file&file=GSM6611297%5FP14601%5F4004%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611297_P14601_4004_matrix.mtx.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611297&format=file&file=GSM6611297%5FP14601%5F4004%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611297_P14601_4004_barcodes.tsv.gz
gunzip GSM6611297_P14601_4004_features.tsv.gz
gunzip GSM6611297_P14601_4004_matrix.mtx.gz

# get sample 2 data (post) then decompress files
wget -O GSM6611298_P14601_4005_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611298&format=file&file=GSM6611298%5FP14601%5F4005%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611298_P14601_4005_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611298&format=file&file=GSM6611298%5FP14601%5F4005%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611298_P14601_4005_matrix.mtx.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611298&format=file&file=GSM6611298%5FP14601%5F4005%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611298_P14601_4005_barcodes.tsv.gz
gunzip GSM6611298_P14601_4005_features.tsv.gz
gunzip GSM6611298_P14601_4005_matrix.mtx.gz

# get sample 3 data (pre) then decompress files
wget -O GSM6611299_P15306_5003_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611299&format=file&file=GSM6611299%5FP15306%5F5003%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611299_P15306_5003_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611299&format=file&file=GSM6611299%5FP15306%5F5003%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611299_P15306_5003_matrix.mtx.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611299&format=file&file=GSM6611299%5FP15306%5F5003%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611299_P15306_5003_barcodes.tsv.gz
gunzip GSM6611299_P15306_5003_features.tsv.gz
gunzip GSM6611299_P15306_5003_matrix.mtx.gz

# get sample 3 data (post) then decompress files
wget -O GSM6611300_P15306_5004_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611300&format=file&file=GSM6611300%5FP15306%5F5004%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611300_P15306_5004_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611300&format=file&file=GSM6611300%5FP15306%5F5004%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611300_P15306_5004_matrix.mtx.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611300&format=file&file=GSM6611300%5FP15306%5F5004%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611300_P15306_5004_barcodes.tsv.gz
gunzip GSM6611300_P15306_5004_features.tsv.gz
gunzip GSM6611300_P15306_5004_matrix.mtx.gz

# get sample 3 data Deep Seq (pre) then decompress files
wget -O GSM6611301_10X_20_003_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611301&format=file&file=GSM6611301%5F10X%5F20%5F003%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611301_10X_20_003_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611301&format=file&file=GSM6611301%5F10X%5F20%5F003%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611301_10X_20_003_matrix.mtx.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611301&format=file&file=GSM6611301%5F10X%5F20%5F003%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611301_10X_20_003_barcodes.tsv.gz
gunzip GSM6611301_10X_20_003_features.tsv.gz
gunzip GSM6611301_10X_20_003_matrix.mtx.gz

# get sample 3 data Deep Seq (post) then decompress files
wget -O GSM6611302_10X_20_004_barcodes.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611302&format=file&file=GSM6611302%5F10X%5F20%5F004%5Fbarcodes%2Etsv%2Egz'
wget -O GSM6611302_10X_20_004_features.tsv.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611302&format=file&file=GSM6611302%5F10X%5F20%5F004%5Ffeatures%2Etsv%2Egz'
wget -O GSM6611302_10X_20_004_matrix.mtx.gz 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6611302&format=file&file=GSM6611302%5F10X%5F20%5F004%5Fmatrix%2Emtx%2Egz'
gunzip GSM6611302_10X_20_004_barcodes.tsv.gz
gunzip GSM6611302_10X_20_004_features.tsv.gz 
gunzip GSM6611302_10X_20_004_matrix.mtx.gz
