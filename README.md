# VAE analysis using DNA methylation

## How to use this repository
In order to use this repository, you need to complete a few steps to ensure we are keeping the code clean. 

1. Create a fork of this repository for your own GitHub account. DO NOT edit the Christensen-Lab-Dartmouth repository directly.
1. Once you have created your own fork, edit the code and documentation that you would like to change. 
1. When you want to commit your changes back to the main repository, create a pull request and then assign the review to AlexanderTitus initially. 

## Analysis Steps
1. Select top N most variable CpGs by Median Absolute Deviation (MAD) - "1-data_preparation.ipynb"
1. Train the VAE model - "2-Model_training.ipynb"
1. Conducte t-SNE dimensionality reduction - "3-tSNE.ipynb"
1. Visualize the t-SNE results - "4-data_visualization.ipynb"
