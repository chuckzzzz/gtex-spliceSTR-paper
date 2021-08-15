import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os
from tqdm import tqdm
plt.rcParams["figure.figsize"] = (20,12)
import argparse
import gzip
import math
# from subprocess import Popen, PIPE, DEVNULL
import sys
import tabix
import statsmodels.api as sm
import subprocess
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f",help='Input STR file', required=True)
args = parser.parse_args()
F=args.f

SNP_GNTP_TEMPLATE="/storage/ydong/data/snp/chr*SNP.recode.vcf.gz"
# SNP_GNTP_TEMPLATE="/storage/czhang/projects/STR_expression/data/QC_hg38_SNP_genotypes/chr*.gz"
#SNP_GNTP_TEMPLATE="/storage/czhang/projects/STR_expression/data/genotypes/GTExNormalizedSNPGenotypes_chr*.table.gz"
STR_GNTP_TEMPLATE="/storage/czhang/projects/STR_expression/data/processed_STRs/*.csv"
PSI_TEMPLATE="/storage/czhang/projects/STR_expression/data/PSI/tibial_nerve/chr*_PSI.csv"

def GetZ(xvals, yvals):
    X = np.array(xvals)
    Y = np.array(yvals)
    X = sm.add_constant(X)
    mod_ols = sm.OLS(Y, X, missing="drop")
    res_ols = mod_ols.fit()
    # not sure why this happens, but sometimes the params or base only has 1 value, which causes error
    if(len(res_ols.params)==1 or len(res_ols.bse)==1):
        return np.nan
    return res_ols.params[1]/res_ols.bse[1]

def generate_GNTP_files(chrom):
    chrom=str(chrom)
    SNP_f=SNP_GNTP_TEMPLATE.replace("*",chrom)
    STR_f=STR_GNTP_TEMPLATE.replace("*",chrom)
    PSI_f=PSI_TEMPLATE.replace("*",chrom)
    return SNP_f,STR_f,PSI_f

def get_snp_samples(snp_header_file):
    # p=subprocess.Popen(["zcat",SNP_gntp_f,"|","head", "-n", "1"],shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    with open(snp_header_file) as f:
        headers=f.readlines()[0]
    # headers=headers.split('\t')[2:]
    headers=headers.split('\t')[9:]
    return headers

def get_common_samples_for_PSI(psi_df,cur_gene,snp_samples):
    psi_samples=list(psi_df.index)
    samples_to_keep=[x for x in psi_samples if x in snp_samples]
    return samples_to_keep,psi_df

def search_SNPs_around(SNP_gntp_f, window_size,STR_id):
    chrom,pos=STR_id.split("_")
    # temporary
    if(chrom=="chr1"):
        SNP_gntp_f="/storage/ydong/data/snp/chr1SNP_NoSegDup.recode.vcf.gz"
    tb=tabix.open(SNP_gntp_f)
    pos=int(pos)
    records=tb.query(chrom,int(pos-window_size),int(pos+window_size))
    return records

def run(target_file, window=1e4, snp_header_file="/storage/czhang/projects/STR_expression/data/QC_hg38_SNP_genotypes/header_temp.tsv",output_dir="/storage/czhang/projects/STR_expression/pipeline/CAVIAR/test/"):
    df=pd.read_csv(target_file,sep='\t')
    cur_chr_ls=df['chrom'].values
    cur_str_ls=df['str.id'].values
    cur_str_ls=np.array(["_".join(["chr"+str(x),y.split("_")[1]]) for x,y in zip(cur_chr_ls,cur_str_ls)])
    cur_str_coord=np.array([int(x.split('_')[1]) for x in cur_str_ls])
    cur_gene_ls=df['gene'].values
    psi_df=None
    psi_df_chr=None
    for cur_chr,cur_str,cur_str_coordinate,cur_gene in tqdm(zip(cur_chr_ls,cur_str_ls,cur_str_coord,cur_gene_ls),total=len(cur_chr_ls)):
        # if(cur_chr!=1):
        #     continue
        SNP_gntp_f,STR_gntp_f,PSI_f=generate_GNTP_files(cur_chr)
        
        if(cur_chr!=psi_df_chr):
            print("Changing current psi chr from {} to {}!".format(str(psi_df_chr),str(cur_chr)))
            psi_df=pd.read_csv(PSI_f,index_col = "Unnamed: 0")
            psi_df_chr=cur_chr
        
        # some genes are not in the psi df, skip those
        if(cur_gene not in psi_df.columns):
            continue
        
        snp_samples=get_snp_samples(snp_header_file)

        samples_to_keep,psi_df=get_common_samples_for_PSI(psi_df,cur_gene,snp_samples)
        cur_exon_psi=psi_df.loc[samples_to_keep,cur_gene]
        # print(cur_exon_psi)

        records=search_SNPs_around(SNP_gntp_f,window,cur_str)

        ld_df = pd.read_csv(STR_gntp_f,dtype = {"start":int},sep='\t')
        ld_df['key']=ld_df.apply(lambda row:row['chrom']+'_'+str(row['start']),axis=1)
        ld_df=ld_df[ld_df['key']==cur_str].T

        # drop key row
        ld_df=ld_df[:-1]
        ld_df=ld_df.drop(['chrom','start'],axis=0)

        ld_df.columns=[cur_str]
        # print(ld_df)

        # calculate STR linear regression Zscore
        cur_str_row = pd.read_csv(STR_gntp_f,dtype = {"start":int},sep='\t')
        cur_str_row['key']=cur_str_row.apply(lambda row:row['chrom']+'_'+str(row['start']),axis=1)
        cur_str_row=cur_str_row[cur_str_row['key']==cur_str].T
        cur_str_row=cur_str_row.drop(['chrom','start','key'])
        cur_str_samples=cur_str_row.index
        cur_str_data=pd.concat([cur_str_row,psi_df[cur_gene]],axis=1,sort=True).dropna()
        cur_str_data.columns=[cur_str,cur_gene]
        cur_str_data=cur_str_data.astype(float)
        # store Z score
        z_df=[]

        PSI_range=cur_str_data[cur_gene].max()-cur_str_data[cur_gene].min()
        if(PSI_range<10):
            continue

        cur_str_z=[cur_str,GetZ(cur_str_data[cur_str],cur_str_data[cur_gene])]
        z_df.append(cur_str_z)

        z_df_SNP_dict={}

        for cur_record in records:
            # print(cur_record)
            cur_SNP='_'.join(['SNP',str(cur_record[0]),str(cur_record[1])])
            # print(snp_samples)
            cur_snp_df=pd.DataFrame([cur_record[9:]]).T # first two cols are not useful, change to 2 after recode SNP data is available
            cur_snp_df.index=snp_samples

            cur_snp_df=cur_snp_df.loc[samples_to_keep,:]
            cur_snp_df=pd.concat([cur_snp_df,cur_exon_psi],axis=1)
            cur_snp_df.columns=[cur_SNP,'PSI'] 

            ### Temperary, update after apply recode to the vcf file ###
            cur_snp_df[cur_SNP]=cur_snp_df.apply(lambda row:row[cur_SNP].split(':')[0],axis=1)
            cur_snp_df=cur_snp_df[~cur_snp_df[cur_SNP].isin(["NA"])]
            cur_snp_df=cur_snp_df.astype(float)

            # # keep only rows with valid genotypes
            # cur_snp_df=cur_snp_df[~cur_snp_df[cur_SNP].isin(["NA"])]

            # # cast to float
            # # print(cur_snp_df)
            # cur_snp_df=cur_snp_df.astype(float)
            
            cur_Z=GetZ(cur_snp_df[cur_SNP],cur_snp_df['PSI'])
            if(np.isnan(cur_Z)):
                continue
            
            # if after removing nan values, there are only 1 genotype, drop this SNP
            test_df=pd.concat([ld_df[cur_str],cur_snp_df[cur_SNP]],axis=1,sort=True)
            test_df=test_df[~(test_df[cur_str].isnull())|(test_df[cur_SNP]).isnull()]
            unique_genotypes=test_df[cur_SNP].nunique()
            # if(cur_SNP=="SNP_chr1_36152264"):
            #     print(unique_genotypes)
            #     print(cur_snp_df)
            if(unique_genotypes<2):
                # print(cur_SNP)
                continue
            ld_df=pd.concat([ld_df,cur_snp_df[cur_SNP]],axis=1,sort=True)
            z_df_SNP_dict[cur_SNP]=[cur_SNP,cur_Z]
            # z_df.append([cur_SNP,cur_Z])

        # further discard SNPs that has NaN correlation with other SNPs
        ld_df=ld_df.astype(float)
        ld_df.to_csv(output_dir+"t.tsv",sep="\t")
        ld_mtx=ld_df.corr()
        is_NaN=ld_mtx.isnull()
        row_has_NaN = is_NaN.any(axis=1)
        snp_to_discard = ld_mtx[row_has_NaN].index.values
        # print(len(snp_to_discard))
        # print(snp_to_discard)
        ld_mtx=ld_mtx.drop(snp_to_discard,axis=1)
        # print(ld_mtx.shape)
        ld_mtx=ld_mtx.drop(snp_to_discard,axis=0)
        # print(ld_mtx.shape)
        ld_f=output_dir+"_".join([cur_str,cur_gene,])+".LDFILE"
        ld_mtx.to_csv(ld_f,header=False,index=False,sep="\t")

        snps_to_keep=ld_mtx.index.values[1:] # first is str
        for snp in snps_to_keep:
            z_df.append(z_df_SNP_dict[snp])

        z_df=pd.DataFrame(z_df)
        z_f=output_dir+"_".join([cur_str,cur_gene,])+".ZFILE"
        z_df.to_csv(z_f,header=False,index=False,sep="\t")
        # print(ld_df)

        # print(ld_mtx.shape)
        # ld_mtx_columns=ld_mtx.columns
        # ld_df_columns=ld_df.columns
        # print(ld_df_test)
        # print(ld_df_test.corr().shape)
        # diff=[x for x in ld_df_columns if x not in ld_mtx_columns]
        # print(diff)


        # CAVIAR
        output_prefix=output_dir+"_".join([cur_str,cur_gene])
        msg=subprocess.check_output(["CAVIAR","-l",ld_f,"-z",z_f,"-o",output_prefix])




if __name__=="__main__":
    run(F)