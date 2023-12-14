import pandas as pd 
import numpy as np 
from tqdm import tqdm,trange
import argparse
import os
import warnings
from ast import literal_eval
from len2bin import change_format,generate_bin
warnings.filterwarnings("ignore")


def three_region(input,output,file_version):
    input.rename(columns={"mRNA":"ref_id"},inplace = True)
    input["5'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["5'UTR_length"]))
    input["3'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["3'UTR_length"]))
    input["ALL"] = list(map(lambda x,y,z:f'1-{x+y+z}',input["5'UTR_total"],input["3'UTR_total"],input["cds_length"]))
    input["UTR5"] = list(map(lambda x:'' if x == 0 else f'1-{x}',input["5'UTR_total"]))
    input["CDS"] = list(map(lambda x,y:f'{x+1}-{x+y}',input["5'UTR_total"],input["cds_length"]))
    input['UTR3'] = list(map(lambda x,y,z: '' if y == 0 else f'{x+z+1}-{x+y+z}',input["5'UTR_total"],input["3'UTR_total"],input["cds_length"]))
    input = input[["ref_id","ALL","UTR5",'CDS',"UTR3"]]
    input.sort_values(by="ref_id",inplace = True)
    input.to_csv(f'{output}/mRNA_{file_version}_3region.csv',index = False)
    return

def length(input,output,file_version):
    input.rename(columns={"mRNA":"ref_id"},inplace = True)
    input["5'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["5'UTR_length"]))
    input["3'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["3'UTR_length"]))
    input["length"] = list(map(lambda x,y,z:f'{x+y+z}',input["5'UTR_total"],input["3'UTR_total"],input["cds_length"]))
    input = input[["ref_id","length"]]
    input.sort_values(by="ref_id",inplace = True)
    input.to_csv(f'{output}/mRNA_{file_version}_length.csv',index = False)
    return

def id_to_name(input,output,file_version,mRNA_file):
    name_id = {}
    with open(f"{mRNA_file}") as f:
            for row in f:
                if row.startswith(">"):
                    rowsplit = row.split(" ")
                    name_id[rowsplit[0].strip(">")] = rowsplit[1].split('=')[1].strip("\n")
    input.rename(columns={"mRNA":"Gene name"},inplace = True)
    input["Gene ID"] = list(map(lambda x:name_id[x],input["Gene name"]))
    input=input[["Gene name","Gene ID"]]
    input.sort_values(by="Gene name",inplace = True)
    input.to_csv(f'{output}/mRNA_{file_version}_IDtoName.csv',index = False)
    return


def boundary(input,output,file_version):
    input.rename(columns={"mRNA":"ref_id"},inplace = True)
    input["5'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["5'UTR_length"]))
    input["3'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["3'UTR_length"]))
    input['head'] = np.ones(len(input))
    input['tail'] = list(map(lambda x,y,z:f'{x+y+z}',input["5'UTR_total"],input["3'UTR_total"],input["cds_length"]))
    input['start codon'] = list(map(lambda x:f'{x+1}',input["5'UTR_total"]))
    input['stop codon'] = list(map(lambda x,y:f'{x+y}',input["5'UTR_total"],input["cds_length"]))
    input = input[['ref_id','head','tail','start codon','stop codon']]
    input.sort_values(by='ref_id',inplace = True)
    input.to_csv(f'{output}/mRNA_{file_version}_boundary.csv',index = False)
    return

def metagene(input,output,file_version):
    input.rename(columns={"mRNA":"ref_id"},inplace = True)
    input["5'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["5'UTR_length"]))
    input["3'UTR_total"] = list(map(lambda x:sum(literal_eval(x)),input["3'UTR_length"]))
    input["length"] = list(map(lambda x,y,z:int(x+y+z),input["5'UTR_total"],input["3'UTR_total"],input["cds_length"]))
    input.sort_values(by='ref_id',inplace = True)
    metadf = change_format(input,100)
    metadf.sort_values(by='ref_id',inplace = True)
    metadf.to_csv(f'{output}/mRNA_{file_version}_metagene.csv',index = False)
    return


if __name__ =='__main__':
    parser = argparse.ArgumentParser(
        prog="mRNA_sRNAanalyzer_ref.py",
        description="This code is to make the reference file for analysis(plot)",
        epilog="This code can be seen as step4, use directly after worm_mRNA.py, three_region.py and intron_analysis.py"
    )
    parser.add_argument("-i","--input",help="This is argument for the output file of the intron_analysis.py")
    parser.add_argument("-f","--file_version",help="This is argument for the version of data")
    parser.add_argument("-o","--output",help="This is argument for the path that save output file")
    parser.add_argument("-t","--type",help="This is argument for the type that output we want to make, include 3region,boundary,IDtoName,length,metagene")
    parser.add_argument("-m","--mRNA_file",nargs="?",help="This is argument for the mRNA_file that trnasfer the id to Name(optional)")
    args = parser.parse_args()
    input = pd.read_csv(args.input)
    output = str(args.output)
    file_version = args.file_version
    file_type = args.type
    if os.path.exists(f'{output}') == False:
        os.mkdir(output)
    if file_type =="3region":
        three_region(input,output,file_version)
    elif file_type =="length":
        length(input,output,file_version)
    elif file_type =="IDtoName":
        mRNA_file = str(args.mRNA_file)
        id_to_name(input,output,file_version,mRNA_file)
    elif file_type=='boundary':
        boundary(input,output,file_version)
    elif file_type=='metagene':
        metagene(input,output,file_version)