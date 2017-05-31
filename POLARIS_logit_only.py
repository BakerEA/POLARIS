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

MPI.COMM_WORLD.Barrier()

###########################
# Run Logistic Regression #
###########################

if rank==0:
    print ("Running Logistic Regression on POLARIS")
f.logit(input_filename, output_filename, covar_filename)

MPI.COMM_WORLD.Barrier()


timer_end= time.time()-timer_start

if rank==0:
    print (timer_end)
