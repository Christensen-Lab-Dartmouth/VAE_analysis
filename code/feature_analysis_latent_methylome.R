###########################
# Analysis of new features generated from a latent methylome
#
# Author: Alexander J. Titus
# Date: September 12, 2017
###########################

### Setup the environment we need ###
# install.packages('glmnet')
# install.packages('data.table')
require(glmnet)
require(data.table)


### Load our data to work with ###
vae.file = '../data/encoded_methyl_onehidden_warmup_batchnorm_100K-100.tsv'
BRCA.covFile = '../data/BRCAtarget_covariates.csv'
vae_features = data.frame(fread(vae.file), row.names=1)
colnames(vae_features) = vae_features[1, ]
vae_features = vae_features[2:nrow(vae_features), ]
vae_features$Basename = rownames(vae_features)

BRCA.covs = data.frame(fread(BRCA.covFile), row.names=1)

full.data = merge(vae_features, BRCA.covs, by='Basename')

### Modeling ### 
