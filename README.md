
# POLARIS- Polygenic LD-adjusted Risk Score #


# SET UP:

To successfully run POLARIS, the following prerequisites are required:
1) MPI
2) ANACONDA (+ Mpi4py)
3) PLINK

Below are generic install instructions for open source software, however, other compilers can be used.

1) MPI:

OpenMPI installation:

This can be downloaded from www.open-mpi.org

Unzip the downloaded folder and save on your machine.

For MacOSX:

In Terminal, go into the downloaded folder and run the following commands:

- ./configure - -prefix=/opt/openmpi 2>&1 | tee config.out

- make -j 4 2>&1 | tee make.out

- sudo make install 2>&1 | tee install.out

- export PATH=/opt/openmpi/bin:$PATH (Add permanently to path for easy use)

- ompi_info


2) ANACONDA:

I used the Anaconda distribution of python which can be downloaded from:
https://www.continuum.io/downloads
Anaconda contains the majority of packages used: sys, numpy, pandas, os, time, math

The mpi4py package must be additionally installed:

pip install mpi4py

3) PLINK:

The POLARIS code uses PLINK2, this must also be installed and available to run from the command line. Source your executable so plink2 runs by typing the command plink2.
PLINK2 can be downloaded from https://www.cog-genomics.org/plink2

# SYNTAX:


FILES:

-- INFILE <filename> specify file name for genotype data in plink binary format (Data can be either all chromosomes combined or split by chromosome, if split, please only use naming structure filename_chr1 etc. for each chromosome and simply put filename as INFILE option)

-- OUTFILE <filename> specify output filename for all files, if none specified, the infile name will be used

-- SUMM <filename> summary data containing columns CHR, SNP, A1, A2, BETA, P, MAF; tab separated. Note A1 and A2 should be uppercase. No missing values.

-- COVAR <filename> file containing population covariates to adjust for in logistic regression model (must include “FID” and “IID” columns, from PLINK .fam, all other columns are to be included in the logistic regression model); tab separated



ANNOTATE SNPS TO GENES:

--ANNOT <filename>  if annotation is already carried out file must contain columns CHR, SNP, and GENE, if --RUN-ANNOT option is specified this file should contain columns CHR, BP_START, BP_END and GENE, Gencode is used by default if no ANNOT option is specified (including all known, protein coding genes) https://www.gencodegenes.org/releases/19.html

-- RUN-ANNOT  <lower bound, upper bound> performs annotation specifying border around the gene in kb i.e. -- RUN-ANNOT 35,10 is 35 kb downstream and 10 kb upstream of the gene 



THRESHOLD SNPS:

-- THR <p threshold> only include SNPs with a p-value less than p-threshold (default p=1 if the option is not specified)


# EXAMPLES:

To run the code on a specified number of processors (nproc) use the following:

mpiexec -np nproc python POLARIS_master.py


To run POLARIS using the gencode annotation file provided use the following syntax: (A different annotation file can be specified using formats described above)

mpiexec -np nproc python POLARIS_master.py --INFILE input_filename  -- OUTFILE output_filename  -- THR pval_thr  -- RUN-ANNOT lower bound,upper bound -- ANNOT gencode_annot -- SUMM summary_data_filename  -- COVAR covariate_filename


To run POLARIS with a pre-defined annotation file formatted as described above use the following syntax:

mpiexec -np nproc python POLARIS_master.py --INFILE input_filename  -- OUTFILE output_filename  -- THR pval_thr  -- ANNOT annotation_filename -- SUMM summary_data_filename  -- COVAR covariate_filename


Using test data and running annotation:

mpiexec -np nproc python POLARIS_master.py --INFILE POLARIS_test_data -- OUTFILE POLARIS_test_data  -- THR pval_thr  -- RUN-ANNOT lower bound,upper bound -- ANNOT gencode_annot -- SUMM POLARIS_test_summ_data -- COVAR covariate_filename

Using test data and an annotation file:

mpiexec -np nproc python POLARIS_master.py --INFILE POLARIS_test_data  -- OUTFILE POLARIS_test_data  -- THR pval_thr  -- ANNOT POLARIS_test_data.genes.annot -- SUMM POLARIS_test_summ_data   -- COVAR covariate_filename

