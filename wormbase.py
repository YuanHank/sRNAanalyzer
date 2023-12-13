import urllib.request as request
import bs4
import json
import pandas as pd
import time 
import numpy as np
import itertools
#this one is for intron only not the total one

def wormbase_crawler(transcript,version):
    '''
    input : transcript --str ，為輸入的transcript 名稱
    output: exon_intron --df，以dataframe形式儲存unspliced情況的各區域資料
            SequenceLength --int,儲存spliced情況的sequence長度
    '''
    url = 'https://wormbase.org/rest/widget/transcript/'+transcript+'/sequences'
    req = request.Request(url,headers ={"User-Agent":"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36"})
    with request.urlopen(req) as response:
        sequence = response.read().decode("utf-8")
    # sequence  = json.loads(sequence)
    with open(f'mRNA_json_{version}/{transcript}.json','w') as f:
        f.write(sequence)
    # try:
    #     strand = sequence['fields']['unspliced_sequence_context']['data']['strand']
    # except:
    #     strand = sequence['fields']['unspliced_sequence_context']['data']['strand']
    # #辨別正反股
    # if strand =='+':
    #     try:
    #         exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context']['data']['positive_strand']['features'])
    #         SequenceLength =len(sequence['fields']['spliced_sequence_context']['data']['positive_strand']['sequence'])
    #     except:
    #         exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context_with_padding']['data']['positive_strand']['features'])
    #         SequenceLength =len(sequence['fields']['spliced_sequence_context_with_padding']['data']['positive_strand']['sequence'])
    # elif strand =='-':
    #     try:
    #         exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context']['data']['negative_strand']['features'])
    #         SequenceLength =len(sequence['fields']['spliced_sequence_context']['data']['negative_strand']['sequence'])
    #     except:
    #         exon_intron = pd.DataFrame(sequence['fields']['unspliced_sequence_context_with_padding']['data']['negative_strand']['features'])
    #         SequenceLength =len(sequence['fields']['spliced_sequence_context_with_padding']['data']['negative_strand']['sequence'])
    return("", "")


## 用來爬歸類在cds的transcript （不確定分類方式，但有多此種類別）
def wormbase_crawler_cds(transcript,version):
    '''
    input : transcript --str ，為輸入的transcript 名稱
    output: exon_intron --df，以dataframe形式儲存unspliced情況的各區域資料
            SequenceLength --int,儲存spliced情況的sequence長度
    '''
    url = 'https://wormbase.org/rest/widget/cds/'+transcript+'/sequences'
    req = request.Request(url,headers ={"User-Agent":"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36"})
    with request.urlopen(req) as response:
        sequence = response.read().decode("utf-8")
    # sequence  = json.loads(sequence)
    # print("this")
    with open(f'mRNA_json_{version}/{transcript}.json','w') as f:
        f.write(sequence)
    # print("that")
    # #print(sequence)
    # strand = sequence['fields']['cds_sequence']['data']['strand']
    # if strand =='+':
    #     exon_intron = pd.DataFrame(sequence['fields']['cds_sequence']['data']['positive_strand']['features'])
    #     SequenceLength = len(sequence['fields']['cds_sequence']['data']['positive_strand']['sequence'])
    # elif strand =='-':
    #     exon_intron = pd.DataFrame(sequence['fields']['cds_sequence']['data']['negative_strand']['features'])
    #     SequenceLength = len(sequence['fields']['cds_sequence']['data']['negative_strand']['sequence'])
    return("", "")    

# def (transcript):
#     '''
#     input : transcript --str ，為輸入的transcript 名稱
#     output: exon_intron --df，以dataframe形式儲存unspliced情況的各區域資料
#             SequenceLength --int,儲存spliced情況的sequence長度
#     '''
#     url = 'https://wormbase.org/rest/widget/transcript/'+transcript+'/sequences'
#     req = request.Request(url,headers ={"User-Agent":"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36"})
#     with request.urlopen(req) as response:
#         sequence = response.read().decode("utf-8")
#     sequence  = json.loads(sequence)




if __name__ =='__main__':
    ExonIntron,SequenceLength = wormbase_crawler(transcript='C05C12.6.1',version = 285)
    # ExonIntron,SequenceLength= wormbase_crawler_cds(transcript= 'B0213.1')
    # print(ExonIntron)
    # print(SequenceLength)
    # IntronCount =(ExonIntron.type == 'exon').sum()-1
    # FivePrimeNumber = 0
    # try:
    #     FivePrime_df = ExonIntron[ExonIntron['type']=='five_prime_UTR']
    #     for j in range(len(FivePrime_df)):
    #         FivePrimeNumber += FivePrime_df.iloc[j]['stop'] - FivePrime_df.iloc[j]['start'] + 1
    # except:
    #     FivePrimeNumber = 0
