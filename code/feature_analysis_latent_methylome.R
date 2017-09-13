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


### Heatmap ###
require(heatmap3)
heatmap3(vae_features)


### LASSO Regression ### 
temp.data = full.data[,2:14 ]
temp.data = temp.data[, -which(colnames(temp.data) %in% c("X_SAMPLE_ID", "X"))]
temp.data = temp.data[!is.na(temp.data$PAM50.RNAseq), ]
temp.data = temp.data[temp.data$PAM50.RNAseq != '', ]
colnames(temp.data) = c( "one", "two", "three", "four", "five",  
                         "six", "seven", "eight", "nine", "ten",
                         "PAM50.RNAseq", "randu")

x <- as.matrix(model.matrix(PAM50.RNAseq~., temp.data)[,-1])
y <- temp.data$PAM50.RNAseq
lambda <- 10^seq(10, -2, length = 100)

set.seed(489)
train = sample(1:nrow(x), nrow(x)/2)
test = (-train)
ytest = y[test]

lasso.mod <- cv.glmnet(x[train,], y[train], family='multinomial', alpha=1, standardize=TRUE)
plot(lasso.mod)
plot(lasso.mod$glmnet.fit, xvar="lambda", label=TRUE)
lasso.mod$lambda.min
lasso.mod$lambda.1se
coeffic = coef(lasso.mod, s=lasso.mod$lambda.min)

lasso.pred <- predict(lasso.mod, newx = x[test,])

clusters = kmeans(lasso.pred, 5, nstart = 20)
clusters

clusters$cluster <- as.factor(clusters$cluster)
