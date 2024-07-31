import pandas as pd 
import json 
import numpy as np 
from tqdm import tqdm,trange
import time
import argparse
import os
import warnings
warnings.filterwarnings("ignore")
"""
step 3.
This .py file is to analysis the intron,exon,3'UTR,5'UTR data from given .json file name as transcript name,
This will work nicely with three_region.py and worm_mRNA.py file.

input
    -- {version}/{transcript name}.json
    -- .csv file with given cds transcript as column named as transcript 
    -- data version
output
    -- intron_analysis_output/intron_analysis_{data version}.csv
"""
def organize_transcript(data,cds_list):
    ThreePrimeNumber = 0
    FivePrimeNumber = 0
    with open (data) as json_file:
        sequence  = json.load(json_file)
    try:
        strand = sequence['fields']['unspliced_sequence_context']['data']['strand']
    except:
        strand = sequence['fields']['unspliced_sequence_context']['data']['strand']
    #辨別正反股
    if strand =='+':
        try:
            exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context']['data']['positive_strand']['features'])
            sequence_length =len(sequence['fields']['unspliced_sequence_context']['data']['positive_strand']['sequence'])
            spliced_sequence_length = len(sequence['fields']['spliced_sequence_context']['data']['positive_strand']['sequence'])
        except:
            exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context_with_padding']['data']['positive_strand']['features'])
            sequence_length =len(sequence['fields']['unspliced_sequence_context_with_padding']['data']['positive_strand']['sequence'])
            spliced_sequence_length = len(sequence['fields']['spliced_sequence_context_with_padding']['data']['positive_strand']['sequence'])
    elif strand =='-':
        try:
            exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context']['data']['negative_strand']['features'])
            sequence_length =len(sequence['fields']['unspliced_sequence_context']['data']['negative_strand']['sequence'])
            spliced_sequence_length = len(sequence['fields']['spliced_sequence_context']['data']['negative_strand']['sequence'])
        except:
            exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context_with_padding']['data']['negative_strand']['features'])
            sequence_length =len(sequence['fields']['unspliced_sequence_context_with_padding']['data']['negative_strand']['sequence'])
            spliced_sequence_length = len(sequence['fields']['spliced_sequence_context_with_padding']['data']['negative_strand']['sequence'])
    exon_intron['length'] = list(map(lambda x,y:x-y+1,exon_intron['stop'],exon_intron['start']))
    exon_intron['region'] = list(map(lambda x,y:f'{y}-{x}',exon_intron['stop'],exon_intron['start']))
    try:
        FivePrime_df = exon_intron[exon_intron['type']=='five_prime_UTR']
        for j in range(len(FivePrime_df)):
            FivePrimeNumber += FivePrime_df.iloc[j]['stop'] - FivePrime_df.iloc[j]['start'] + 1
    except:
        FivePrimeNumber = 0
    # section 2
    try :
        ThreePrime_df = exon_intron[exon_intron['type']=='three_prime_UTR']
        for j in range(len(ThreePrime_df)):
            ThreePrimeNumber += ThreePrime_df.iloc[j]['stop'] - ThreePrime_df.iloc[j]['start'] + 1
    except:
        ThreePrimeNumber = 0 
    if data.split('/')[1].strip('.json') in cds_list:
        CDSLengthNumber = spliced_sequence_length - FivePrimeNumber - ThreePrimeNumber
    else:
        CDSLengthNumber = 0
    return(exon_intron,sequence_length,FivePrimeNumber,ThreePrimeNumber,CDSLengthNumber)

def organize_cds(data,cds_list):
    ThreePrimeNumber = 0
    FivePrimeNumber = 0
    with open (data) as json_file:
        sequence  = json.load(json_file)    
    strand = sequence['fields']['cds_sequence']['data']['strand']
    if strand =='+':
        exon_intron = pd.DataFrame(sequence['fields']['cds_sequence']['data']['positive_strand']['features'])
        sequence_length = len(sequence['fields']['cds_sequence']['data']['positive_strand']['sequence'])
        spliced_sequence_length = len(sequence['fields']['cds_sequence']['data']['positive_strand']['sequence'])

    elif strand =='-':
        exon_intron = pd.DataFrame(sequence['fields']['cds_sequence']['data']['negative_strand']['features'])
        sequence_length = len(sequence['fields']['cds_sequence']['data']['negative_strand']['sequence'])
        spliced_sequence_length = len(sequence['fields']['cds_sequence']['data']['negative_strand']['sequence'])

    exon_intron['length'] = list(map(lambda x,y:x-y+1,exon_intron['stop'],exon_intron['start']))
    exon_intron['region'] = list(map(lambda x,y:f'{y}-{x}',exon_intron['stop'],exon_intron['start']))
    try:
        FivePrime_df = exon_intron[exon_intron['type']=='five_prime_UTR']
        for j in range(len(FivePrime_df)):
            FivePrimeNumber += FivePrime_df.iloc[j]['stop'] - FivePrime_df.iloc[j]['start'] + 1
    except:
        FivePrimeNumber = 0
    # section 2
    try :
        ThreePrime_df = exon_intron[exon_intron['type']=='three_prime_UTR']
        for j in range(len(ThreePrime_df)):
            ThreePrimeNumber += ThreePrime_df.iloc[j]['stop'] - ThreePrime_df.iloc[j]['start'] + 1
    except:
        ThreePrimeNumber = 0 
    if data.split('/')[1].strip('.json') in cds_list:
        CDSLengthNumber = spliced_sequence_length - FivePrimeNumber - ThreePrimeNumber
    else:
        CDSLengthNumber = 0
    return(exon_intron,sequence_length,FivePrimeNumber,ThreePrimeNumber,CDSLengthNumber)

def filter_length(data,type,length):
    data = data[data['type'] ==type][data['length']>=length]
    return(data)
    
if __name__ =="__main__":
    parser = argparse.ArgumentParser(
        prog="intron_analysis.py",
        description="This code is to analysis the intron,exon,3'UTR,5'UTR data from given .json file name as transcript name",
        epilog="This code can be seen as step3, use directly after worm_mRNA.py and three_region.py"
    )
    parser.add_argument("--input",help="This is argument for the path that save the mRNA.json files(path)")
    parser.add_argument("--input2",help="This is argument for the .csv file that include the cds transcript(.csv)")
    parser.add_argument("--output",help="This is argument for the path that save the output files(path)")
    parser.add_argument("--file_version",help="This is argument for the version of data")
    parser.add_argument("--length_filter",help="This is argument for the length limit for filter the region data")
    args = parser.parse_args()
    input_path = str(args.input)
    transcript_list = os.listdir(input_path)
    transcript_total = []
    cds_file = pd.read_csv(str(args.input2))
    cds_list = list(cds_file['transcript'])
    for i in trange(len(transcript_list)):
        try:
            exon_intron,sequence_length,FivePrimeNumber,ThreePrimeNumber,CDSLengthNumber = organize_transcript(f'{input_path}/{transcript_list[i]}',cds_list)
        except:
            exon_intron,sequence_length,FivePrimeNumber,ThreePrimeNumber,CDSLengthNumber= organize_cds(f'{input_path}/{transcript_list[i]}',cds_list)
        intron = filter_length(exon_intron,'intron',int(args.length_filter))
        exon = filter_length(exon_intron,'exon',int(args.length_filter))
        five_prime = filter_length(exon_intron,'five_prime_UTR',int(args.length_filter))
        three_prime = filter_length(exon_intron,'three_prime_UTR',int(args.length_filter))
        transcript_total.append({
                                'mRNA':transcript_list[i].strip('.json'),
                                "start-end(5'UTR)":list(five_prime['region']),
                                "5'UTR_length":list(five_prime['length']),
                                'start-end(intron)':list(intron['region']),
                                'intron_length':list(intron['length']),
                                'start-end(exon)':list(exon['region']),
                                'exon_length':list(exon['length']),
                                "start-end(3'UTR)":list(three_prime['region']),
                                "3'UTR_length":list(three_prime['length']),
                                'intron_count':len(intron),
                                "sequence_length":sequence_length,
                                "cds_length":CDSLengthNumber,
                                })
    
    transcript_total_table = pd.DataFrame.from_dict(transcript_total)
    transcript_total_table.to_csv(f'{args.output}/intron_analysis_{args.file_version}.csv')
