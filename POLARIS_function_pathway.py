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

import pandas as pd
import numpy as np
from numpy import linalg as LA
from scipy import stats
import statsmodels.api as sm
import math
from mpi4py import MPI
from mpi4py.MPI import ANY_SOURCE

import time

np.set_printoptions(suppress=True)

#######################
# Annotation function #
#######################

def annotation(input_filename, output_filename, annot, l_border, u_border):

    annot=pd.read_csv(annot, sep='\t')
    
    filename= str(input_filename) + ".bim"
    
    snps=pd.read_csv(filename, header=None, sep='\t')

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    
    a= int(math.floor(len(annot.index)/size))
    b=int(len(annot.index)-(a*size))
    
    if ((rank+1)<=b):
        start_gene=((rank)*a)+ (rank+1)
        end_gene= start_gene + a
    elif ((rank+1)==(b+1)):
        start_gene= ((rank)*a)+(rank+1)
        end_gene= start_gene + a - 1
    else:
        start_gene= (b*a) + (b+1) + ((rank-b)*(a-1)) + (rank-b)
        end_gene= start_gene + a - 1

    MPI.COMM_WORLD.Barrier()

    results=pd.DataFrame()
    #results_list=[]
  
    for gene in range(start_gene-1, end_gene):
        
        gene_name=annot.iloc[gene]['GENE']
        gene_chr=annot.iloc[gene]['CHR']
        gene_start=annot.iloc[gene]['BP_START']
        gene_end=annot.iloc[gene]['BP_END']

        gene_info= str(gene_name)+"_"+str(gene_chr)+"_"+str(gene_start)+"_"+str(gene_end)

        gene_snps=snps.loc[(snps[0]==gene_chr) & (snps[3]>=(gene_start-float(l_border))) & (snps[3]<=(gene_end+float(u_border)))]

        #print(gene_snps)

        if len(gene_snps.index)!=0:
            gene_snps.loc[:,6]= str(gene_info)

        del gene_snps[2]

        if len(gene_snps.index)!=0:
            results=pd.concat([results,gene_snps], axis=0)

        #print(results)

        results_list=results.values.tolist()

    del results

    MPI.COMM_WORLD.Barrier()
            
    if rank==0:
    
        for proc in range(1,size):
            
            local_a=MPI.COMM_WORLD.recv(source=proc)
            
            results_list.extend(local_a)

            del local_a
    else:
        
        MPI.COMM_WORLD.send(results_list, dest=0)

    MPI.COMM_WORLD.Barrier()

    if rank==0:
        results=pd.DataFrame(results_list)

        results.columns=['CHR', 'SNP', 'BP', 'A1', 'A2', 'GENE']

        filename= str(output_filename)+ ".annot"
 
        results.to_csv(filename, header=True, index=None, sep='\t')

############################
# Allele Matching Function #
############################

def allele_match(summ_filename, input_filename, output_filename, thr):
    
    # Read in Unique Gene file #
    summ=pd.read_csv(summ_filename, sep='\t')
    #print (summ)
    #print (len(summ.index))
    
    summ=summ.drop_duplicates(subset=['SNP'], keep='last')
    
    #summ['A1'] = map(lambda x: x.upper(), summ['A1'])
    #summ['A2'] = map(lambda x: x.upper(), summ['A2'])
    
    summ=summ[summ.P<=thr]
    
    filename= input_filename + ".bim"
    
    repl=pd.read_csv(filename, header=None, sep='\t')
    
    #print (repl)
    
    data=repl.merge(summ, left_on=[1], right_on=['SNP'])
    
    #print (data)
    #print (len(data.index))
    
    for i in range(len(data.index)):
        if ((data.iloc[i][4]==data.iloc[i]['A2']) & (data.iloc[i][5]==data.iloc[i]['A1'])):
            data.set_value(i, 'BETA', -(data.iloc[i]['BETA']))
            data.set_value(i, 'A1', str((data.iloc[i][4])))
            data.set_value(i, 'A2', str((data.iloc[i][5])))

    for i in range(len(data.index)):
        if ((data.iloc[i][4]!=data.iloc[i]['A1']) | (data.iloc[i][5]!=data.iloc[i]['A2'])):
            if data.iloc[i]['A1']=="A":
                data.set_value(i, 'A1', "T")
            elif data.iloc[i]['A1']=="C":
                data.set_value(i, 'A1', "G")
            elif data.iloc[i]['A1']=="G":
                data.set_value(i, 'A1', "C")
            elif data.iloc[i]['A1']=="T":
                data.set_value(i, 'A1', "A")
            if data.iloc[i]['A2']=="A":
                data.set_value(i, 'A2', "T")
            elif data.iloc[i]['A2']=="C":
                data.set_value(i, 'A2', "G")
            elif data.iloc[i]['A2']=="G":
                data.set_value(i, 'A2', "C")
            elif data.iloc[i]['A2']=="T":
                data.set_value(i, 'A2', "A")

    for i in range(len(data.index)):
        if ((data.iloc[i][4]==data.iloc[i]['A2']) & (data.iloc[i][5]==data.iloc[i]['A1'])):
            data.set_value(i, 'BETA', -(data.iloc[i]['BETA']))
            data.set_value(i, 'A1', str((data.iloc[i][4])))
            data.set_value(i, 'A2', str((data.iloc[i][5])))

    data=data.loc[(data[4]==data['A1']) & (data[5]==data['A2'])]
    
    data=data[data.columns[6:len(data.columns)]]
    
    filename= str(output_filename) + ".summ"
    
    #print (len(data.index))
    
    data.to_csv(filename, header=True, index=None, sep='\t', na_rep="NA")


####################
# POLARIS Function #
####################

def polaris(input_filename, output_filename, annot_filename):
    
    # Read in Unique Gene file #
    unigene_names=["GENE"]
    unique_filename="unique_genes_" + str(output_filename)
    unigene=pd.read_csv(unique_filename, names=unigene_names, sep='\t')
    
    
    # Read in Annot Data #
    annot=pd.read_csv(annot_filename, sep='\t')
    
    filename= str(input_filename) + ".fam"
    
    fam=pd.read_csv(filename, names=['FID', 'IID', 'PID', 'MID', 'SEX', 'PHENOTYPE'], sep=' ')
    fam=fam[['FID', 'IID', 'PHENOTYPE']]
    
    Nind=len(fam.index)
    
    #Split genes by processor
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()

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


    for chr in range(start_chr, end_chr):
        
        filename1="data_" + str(output_filename) + "/" + str(input_filename)+ "_head_chr" + str(chr)
        filename2="data_" + str(output_filename) + "/" + str(input_filename)+ "_chr" + str(chr) + "_nohead.raw"
        filename3="summ/" + str(output_filename)+ "_chr" + str(chr) + ".summ"
            
        #Read in header to determine SNP positions
        infile= open(filename1, 'r')
        firstline=infile.readline()
        head=firstline.split()
        infile.close()
        header=np.array(head)
        
        infile= open(filename3, 'r')
        firstline=infile.readline()
        summ_head=firstline.split()
            
        #Find beta and maf locations in summ data
        for i in range(len(summ_head)):
            if str(summ_head[i])=="BETA":
                beta_loc=i
            if str(summ_head[i])=="MAF":
                maf_loc=i
        
        fs=open(filename3, 'r')
        beta=[]
        maf=[]
        row=0
        for line in fs:
            if row>0:
                beta.append(line.split()[beta_loc])
                maf.append(line.split()[maf_loc])
            row=row+1
        fs.close()
        
        beta=np.asmatrix(beta, dtype=float)
        maf=np.array(maf, dtype=float)
            
        f=open(filename2, 'r')
        data=[]
        for line in f:
            datalist=[]
            testlist=line.split()
            for j in range(len(testlist)):
                datalist.append(testlist[j])
            loc=[i for i, x in enumerate(datalist) if x == "NA"]
            for i in range(len(loc)):
                datalist[loc[i]]=2*float(maf[loc[i]-6])
            data.append(datalist)
        f.close()

        data=pd.DataFrame(data)
    
        filename4="data_" + str(output_filename) + "/" + str(input_filename)+ "_chr" + str(chr) + "_nomiss_nohead.raw"
        data.to_csv(filename4, header=False, index=None, sep=' ')

        MPI.COMM_WORLD.Barrier()



    a= int(math.floor(len(unigene)/size))
    b=int(len(unigene)-(a*size))
        
    if ((rank+1)<=b):
        start_gene=((rank)*a)+ (rank+1)
        end_gene= start_gene + a
    elif ((rank+1)==(b+1)):
        start_gene= ((rank)*a)+(rank+1)
        end_gene= start_gene + a - 1
    else:
        start_gene= (b*a) + (b+1) + ((rank-b)*(a-1)) + (rank-b)
        end_gene= start_gene + a - 1

    MPI.COMM_WORLD.Barrier()
        
    if rank==0:
        results=np.empty([Nind,len(unigene)],dtype=np.float64)

    else:
        results=np.empty([Nind,end_gene-start_gene+1],dtype=np.float64)


    genename=[]

    for gene in range(start_gene-1, end_gene):

        pwscore= np.zeros([Nind,1], dtype=np.float64)

        gene_snp = annot[annot.GENE==unigene.iloc[gene]['GENE']]
        gene_snp=gene_snp.sort_values(['CHR'])
        
        #Output gene chromosome
        gene_chr=gene_snp.drop_duplicates(subset=['CHR'], keep='last')
        gene_chr = gene_chr['CHR']


        for chr in range(len(gene_chr)):

            gene_chr2 = gene_chr.iloc[chr]

            filename1="data_" + str(output_filename) + "/" + str(input_filename)+ "_head_chr" + str(gene_chr2)
            filename2="data_" + str(output_filename) + "/" + str(input_filename)+ "_chr" + str(gene_chr2) + "_nomiss_nohead.raw"
            filename3="summ/" + str(output_filename)+ "_chr" + str(gene_chr2) + ".summ"
            
            #Read in header to determine SNP positions
            infile= open(filename1, 'r')
            firstline=infile.readline()
            head=firstline.split()
            infile.close()
            header=np.array(head)

            snp_loc=[]
            #Find SNP locations in .raw data
            for i in range(len(gene_snp.index)):
                for j in range(6,len(header)):
                    if str(gene_snp.iloc[i]['SNP'])==str(header[j]):
                        snp_loc.append(j)
        
            snp_loc= np.unique(snp_loc)

            gene_name=str(unigene.iloc[gene]['GENE'])

            if len(snp_loc) != 0:
                
                infile= open(filename3, 'r')
                firstline=infile.readline()
                summ_head=firstline.split()
                
                #Find beta and maf locations in summ data
                for i in range(len(summ_head)):
                    if str(summ_head[i])=="BETA":
                        beta_loc=i
                    if str(summ_head[i])=="MAF":
                        maf_loc=i
            
                fs=open(filename3, 'r')
                beta=[]
                maf=[]
                row=-1
                for line in fs:
                    if (row+6) in snp_loc:
                        beta.append(line.split()[beta_loc])
                        maf.append(line.split()[maf_loc])
                    row=row+1
                fs.close()
                
                beta=np.asmatrix(beta, dtype=float)
                maf=np.array(maf, dtype=float)

                f=open(filename2, 'r')
                data=[]
                for line in f:
                    datalist=[]
                    testlist=line.split()
                    for j in range(len(snp_loc)):
                        datalist.append(testlist[snp_loc[j]])
                    data.append(datalist)
                f.close()
                
                data=np.array(data)
                data=np.asmatrix(data, dtype=float)

                if len(snp_loc)==1:
                    
                    score= data * beta
                
                else:
                    
                    corr_mat=np.corrcoef(data, rowvar=0)
                    
                    eval, evect=LA.eig(corr_mat)
                    
                    idx = eval.argsort()[::-1]
                    eval = eval[idx]
                    evect = evect[:,idx]
                    
                    pc_snp=beta * evect
                    
                    pc_wgt=evect
                    for i in range(len(eval)):
                        pc_wgt[:,i] *= math.sqrt((1+(1/math.sqrt(Nind)))/(eval[i]+(1/math.sqrt(Nind))))
                
                    Bi= pc_wgt * pc_snp.transpose()

                    score= data * Bi
            
                #print(score)
                
                pwscore= pwscore+score
                
                #print(pwscore)

        gene_score=pd.DataFrame(pwscore, columns=['GeneName'])
            
        #print(gene_score)

        if rank==0:
            results[:,gene:gene+1]=np.array(gene_score)
        else:
            results[:,gene-start_gene+1:gene-start_gene+2]=np.array(gene_score)

        #print(results)

        genename.append(gene_name)
    
        MPI.COMM_WORLD.Barrier()

        if rank==0:
            
            a= int(math.floor(len(unigene)/size))
            b=int(len(unigene)-(a*size))
            
            for proc in range(1,size):
                if ((proc+1)<=b):
                    start_gene=((proc)*a)+ (proc+1)
                    end_gene= start_gene + a
                elif ((proc+1)==(b+1)):
                    start_gene= ((proc)*a)+(proc+1)
                    end_gene= start_gene + a - 1
                else:
                    start_gene= (b*a) + (b+1) + ((proc-b)*(a-1)) + (proc-b)
                    end_gene= start_gene + a - 1
            
                ncol=(end_gene-start_gene+1)
                
                local_a = np.zeros([Nind,ncol], dtype=np.float64)
                
                MPI.COMM_WORLD.Recv(local_a, source=proc)
                
                results[:,start_gene-1:end_gene]=local_a
        
                del local_a
    
        else:
        
            MPI.COMM_WORLD.Send(results[:,:], dest=0)
        
        MPI.COMM_WORLD.Barrier()
    
        #print(results)

        if rank==0:
            
            for proc in range(1,size):
                
                local_b=MPI.COMM_WORLD.recv(source=proc)
                
                genename.extend(local_b)
                
                #print (genename)
                
                del local_b
    
        else:
            
            MPI.COMM_WORLD.send(genename, dest=0)
    
        MPI.COMM_WORLD.Barrier()
    
        if rank==0:
            
            results=pd.DataFrame(results)
            
            #print(results)
            
            results.columns=genename
            
            final_results=pd.concat([fam, results], axis=1)
            
            #print(final_results)
            
            filename= "results/" + str(output_filename) + "_pathway.polaris"
            
            final_results.to_csv(filename, header=True, index=None, sep='\t')

    MPI.COMM_WORLD.Barrier()



##################
# Logit Function #
##################

def logit(input_filename, output_filename, covar_filename):
    
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()

    log_filename= "log_" + str(output_filename)
    
    # Read in Unique Gene file #
    unigene_names=["GENE"]
    unique_filename="unique_genes_" + str(output_filename)
    unigene=pd.read_csv(unique_filename, names=unigene_names, sep='\t')

    #print(unigene)

    filename= "results/" + str(output_filename) + "_pathway.polaris"


    if rank==0:

        prs_data=pd.read_csv(filename, sep='\t')
        #print(prs_data)


        prs_norm_data = prs_data[['FID', 'IID', 'PHENOTYPE']]
        #print(prs_norm_data)

        for pw in range(len(unigene)):

            prs = prs_data.iloc[:,pw+3].astype(float)
            #print(prs)

            prs_norm = (prs-prs.mean())/(prs.std())

            prs_norm_data=pd.concat([prs_norm_data,prs_norm], axis=1)


        filename= "results/" + str(output_filename) + "_pathway_norm.polaris"

        prs_norm_data.to_csv(filename, header=True, index=None, sep='\t')


        filename= "results/" + str(output_filename) + "_pathway_norm.polaris"

        prs_norm_data=pd.read_csv(filename, sep='\t')

        covar=pd.read_csv(covar_filename, sep='\t')

        no_covar=len(covar.columns)-2
        
        data=prs_norm_data.merge(covar, how='left', on=['FID', 'IID'])

        data['Intercept']=1.0

        #print(data)

        #MPI.COMM_WORLD.Barrier()

        results=pd.DataFrame(index=range(len(unigene)), columns=range(4))


        for gene in range(len(unigene)):

            text=data.columns[gene+3]
            
            #Extract gene info from column header
            
            #print(unigene.iloc[gene])
            results.iat[(gene),0] = unigene.iloc[gene]['GENE']
            
            #print(results)
            
            #Logit regression model
            pos=[]
            pos.append(gene+3)
            
            for i in range(len(unigene)+3, len(data.columns)):
                pos.append(i)
        
            train_cols=data.columns[[pos]]
            
            logit=sm.Logit((data['PHENOTYPE']-1), data[train_cols], missing='drop')

            log_reg=logit.fit()
            #print(log_reg.summary2())
            
            results.iat[(gene),1] = log_reg.params[0]
            results.iat[(gene),2] = log_reg.bse[0]
            results.iat[(gene),3] = log_reg.pvalues[0]


        filename= "results/" + str(output_filename) + "_pathway.polaris.logit"
        
        results.columns=['GENESET','BETA', 'SE', 'P']
        
        results=results.sort_values('P')
        
        results.to_csv(filename, header=True, index=None, sep='\t')


############################
# If on the main processor #
############################

if __name__=='__main__':

    annotation()

    allele_match()

    polaris()

    logit()
