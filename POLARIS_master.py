#!/usr/bin/env python

#############################################################################
#                                                                           #
#               POLARIS: POlygenic Ld- Adjusted RIsk Score                  #
#                                                                           #
#############################################################################
#                                                                           #
#   Copyright (C) 2017 Emily Baker and Cardiff University                   #
#                                                                           #
#   This program is free software: you can redistribute it and/or modify    #
#   it under the terms of the GNU General Public License as published by    #
#   the Free Software Foundation, either version 3 of the License, or       #
#   (at your option) any later version.                                     #
#                                                                           #
#   This program is distributed in the hope that it will be useful,         #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#   GNU General Public License for more details.                            #
#                                                                           #
#   You should have received a copy of the GNU General Public License       #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.   #
#                                                                           #
#   BakerEA@Cardiff.ac.uk                                                   #
#   Emily Baker, MRC Centre for Neuropsychiatric Genetics and Genomics,     #
#   Hadyn Ellis Building, Maindy Road, Cardiff, Wales, UK, CF24 4HQ         #
#                                                                           #
#############################################################################

###########################
#     Import packages     #
###########################

import sys
import POLARIS_function as f
import pandas as pd
import numpy as np
from numpy import linalg as LA
import os
import os.path
import time
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE
import math

###########################
#   Find MPI parameters   #
###########################

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()

timer_start=time.time()

#################################
# Define command line arguments #
#################################

total_arg = len(sys.argv)

#print (sys.argv)

#print (total_arg)

options = ["--INFILE"]
options.append("--OUTFILE")
options.append("--SUMM")
options.append("--ANNOT")
options.append("--RUN-ANNOT")
options.append("--THR")
options.append("--COVAR")

#print (options)

# Check for even number of input parameters
if ((total_arg-1) % 2 != 0):
    print ("Incorrect number of input arguments")
    exit()

# Scan input arguments
input_options={}
input_options['--INFILE']=str("input")
input_options['--OUTFILE']=str("output")
input_options['--SUMM']=str("")
input_options['--ANNOT']=str("gencode_annot")
input_options['--THR']=str("1")

for i in range(1,total_arg):
    if sys.argv[i] in options:
        x=str(sys.argv[i])
        input_options[x]=str(sys.argv[i+1])

#print (input_options)

if (input_options['--INFILE'] == "input"):
    print ("No input file specified")
    exit()

if (input_options['--OUTFILE'] == "output"):
    print ("No output file specified")
    input_options['--OUTFILE'] = input_options['--INFILE']

if (input_options['--SUMM'] == ""):
    print ("No summary statistic file specified")
    exit()

if ('--ANNOT' not in input_options) & ('--RUN-ANNOT' not in input_options):
    print ("No annotation options specified")
    exit()

if (float(input_options['--THR'])>1) | (float(input_options['--THR'])<0):
    print ("Invalid p-value threshold given")
    exit()


input_filename=str(input_options['--INFILE'])
output_filename=str(input_options['--OUTFILE'])
summ_filename=str(input_options['--SUMM'])
annot=str(input_options['--ANNOT'])
thr=float(input_options['--THR'])
covar_filename=str(input_options['--COVAR'])

#print (annot)
#print (input_filename)
#print (output_filename)
#print (summ_filename)
#print (thr)
#print (covar_filename)

###########################
#    Check Data Format    #
###########################


filename= str(input_filename) + "_chr1.fam"
if os.path.isfile(filename):
    by_chr=1
else:
    by_chr=0

if rank==0:
    if by_chr==1:
        
        #Combine bim files
        command="for chr in {1..22} X Y; do cat " + str(input_filename)+ "_chr$chr.bim; done > " + str(input_filename) + ".bim"
        os.system(command)
        
        #Set overall fam file
        command="cat " + str(input_filename) + "_chr1.fam > " + str(input_filename) + ".fam"

MPI.COMM_WORLD.Barrier()

log_filename= "log_" + str(output_filename)

###########################
# Run Annotation Function #
###########################

if '--RUN-ANNOT' in input_options:
    l_border=float(str(input_options['--RUN-ANNOT']).split(',')[0])*1000
    u_border=float(str(input_options['--RUN-ANNOT']).split(',')[1])*1000
    
    if rank==0:
        log=open(log_filename, 'a')
        log.write("Running Annotation")
        log.write("\n")
        log.close()
        print ("Running Annotation")

    f.annotation(input_filename, output_filename, annot, l_border, u_border)

if rank==0:
    log=open(log_filename, 'a')
    log.write("Annotation Complete")
    log.write("\n")
    log.close()
    print ("Annotation Complete")

MPI.COMM_WORLD.Barrier()

###########################
# Create Unique Gene List #
###########################

if ('--RUN-ANNOT' not in input_options) & ('--ANNOT' in input_options):
    annot_filename=annot
else:
    annot_filename= str(output_filename) + ".annot"

annot_data=pd.read_table(annot_filename)

unique_filename= "unique_genes_" + str(output_filename)

uni_genes=annot_data.drop_duplicates(subset=['GENE'], keep='last')
uni_genes['GENE'].to_csv(unique_filename, header=None, index=None, sep='\t')

MPI.COMM_WORLD.Barrier()

###########################
# Perform Allele Matching #
###########################

if rank==0:
    log=open(log_filename, 'a')
    log.write("Perform Allele Matching")
    log.write("\n")
    log.close()
    print ("Perform Allele Matching")
    f.allele_match(summ_filename, input_filename, output_filename, thr)

    log=open(log_filename, 'a')
    log.write("Allele Matching Complete")
    log.write("\n")
    log.close()

MPI.COMM_WORLD.Barrier()

#############################
# Split summary file by chr #
#############################

if rank==0:
    command= "mkdir summ"
    os.system(command)
    
    filename= str(output_filename) + ".summ"
    
    data=pd.read_table(filename)
    
    filename=str(output_filename) + ".summ.snps"
    
    data['SNP'].to_csv(filename, header=None, index=None, sep='\t')
    
    for chr in range(1,23):
        
        chr_data=data[data.CHR==chr]
        filename="./summ/"+ str(output_filename) + "_chr" + str(chr) + ".summ"
        chr_data.to_csv(filename, header=True, index=None, sep='\t', na_rep="NA")

MPI.COMM_WORLD.Barrier()

###########################
#    Recode PLINK data    #
###########################

if rank==0:
    command= "mkdir data_" + str(output_filename)
    os.system(command)

MPI.COMM_WORLD.Barrier()

a= int(math.floor(22/size))
b=int(22-(a*size))

if ((rank+1)<=b):
    start_chr=((rank)*a)+ (rank+1)
    end_chr= start_chr + a
elif ((rank+1)==(b+1)):
    start_chr= ((rank)*a)+(rank+1)
    end_chr= start_chr + a - 1
else:
    start_chr= (b*a) + (b+1) + ((rank-b)*(a-1)) + (rank-b)
    end_chr= start_chr + a - 1

MPI.COMM_WORLD.Barrier()

for chr in range(start_chr,end_chr+1):
    
    if by_chr==0:
        
        command="plink2 --bfile " + str(input_filename) + " --extract " + str(output_filename) + ".summ.snps --recodeA --chr " + str(chr) + " --out ./data_"+ str(output_filename)+"/" + str(input_filename) + "_chr" +str(chr)
        os.system(command)
    
    elif by_chr==1:
        
        command5="plink2 --bfile " + str(input_filename) + "_chr" + str(chr) + " --extract " + str(output_filename) + ".summ.snps --recodeA --out ./data_"+ str(output_filename)+"/" + str(input_filename) + "_chr" + str(chr)
        os.system(command5)

    command2="sed '1d' ./data_"+ str(output_filename)+"/" + str(input_filename)+ "_chr" +str(chr)+ ".raw > ./data_"+ str(output_filename)+"/" + str(input_filename)+ "_chr" +str(chr) + "_nohead.raw"
    os.system(command2)

    command3= "head -n 1 ./data_"+ str(output_filename)+"/" + str(input_filename)+ "_chr" +str(chr)+ ".raw > ./data_"+ str(output_filename)+"/" + str(input_filename)+ "_head_chr" +str(chr)
    os.system(command3)
    
    command4= "rm ./data_"+ str(output_filename)+"/"+ str(input_filename)+ "_chr" +str(chr)+ ".raw"
    os.system(command4)

MPI.COMM_WORLD.Barrier()


for chr in range(start_chr, end_chr+1):
    filename1="data_"+ str(output_filename)+"/" + str(input_filename)+ "_head_chr" + str(chr)
    
    infile= open(filename1, 'r')
    firstline=infile.readline()
    head=firstline.split()
    infile.close()
    
    id=head[0:6]
    head=head[6:len(head)]
    head=[x[:-2] for x in head]
    header= np.array(id + head)
    header=header.tolist()
    
    outfile=open(filename1, 'w')
    outfile.write(" ".join(header))
    outfile.close()

MPI.COMM_WORLD.Barrier()

###########################
#       Run POLARIS       #
###########################

if rank==0:
    command= "mkdir results"
    os.system(command)
    log=open(log_filename, 'a')
    log.write("Running POLARIS")
    log.write("\n")
    log.close()
    print ("Running POLARIS")

f.polaris(input_filename, output_filename, annot_filename)

MPI.COMM_WORLD.Barrier()

if rank==0:
    log=open(log_filename, 'a')
    log.write("POLARIS Complete")
    log.write("\n")
    log.close()

###########################
# Run Logistic Regression #
###########################

if rank==0:
    log=open(log_filename, 'a')
    log.write("Running Logistic Regression")
    log.write("\n")
    log.close()
    print ("Running Logistic Regression on POLARIS")
f.logit(input_filename, output_filename, covar_filename)

MPI.COMM_WORLD.Barrier()

if rank==0:
    log=open(log_filename, 'a')
    log.write("Logit Complete")
    log.write("\n")
    log.close()

###########################
#     Delete old files    #
###########################

if rank==0:
    #command= "rm ./summ/*"
    #os.system(command)
            
    command2= "rm ./data_" + str(output_filename)+ "/*"
    os.system(command2)
        
    #command3= "rmdir summ"
    #os.system(command3)
        
    command4= "rmdir ./data_" + str(output_filename)+ "/*"
    os.system(command4)

MPI.COMM_WORLD.Barrier()

timer_end= time.time()-timer_start

if rank==0:
    print (timer_end)
