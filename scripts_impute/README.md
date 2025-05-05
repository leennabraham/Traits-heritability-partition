DNA Methylation and Trait Heritability Partitioning Using Bayesian Regression

Project Overview

This project aims to explore the role of DNA methylation in the heritability of complex traits in Arabidopsis lyrata by integrating DNA methylation data with phenotypic traits, we utilize Bayesian regression models to partition variance into different components: additive, dominance, and residual variances. This helps to understand the genetic and epigenetic contributions to trait variation in populations.

Objectives
	•	Partitioning trait heritability: Estimating the additive and non-additive components of genetic variation (Va, Vd) using Bayesian regression.
	•	Modeling genetic and epigenetic interactions: Quantifying the effects of DNA methylation on trait variation and how it interacts with the genetic background.
	•	Identifying key genomic regions: Identifying DNA methylation sites that contribute significantly to phenotypic variation in complex traits using ML.


Data Collection

The project uses high-throughput Whole Genome Bisulfite Sequencing (WGBS) data for DNA methylation and RNA-seq data for gene expression. These datasets are complemented by phenotypic data for complex traits such as growth, yield, or disease resistance in Arabidopsis lyrata (or your model species). The dataset consists of [Number of samples] individuals with detailed phenotype measurements.

Bayesian Regression Models

Bayesian regression is employed to model the variance components in complex traits, using the brms package in R. The models are used to estimate the following:
	•	Additive genetic variance (Va)
	•	Dominance variance (Vd)
	•	Residual variance (Vr)