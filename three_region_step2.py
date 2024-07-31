from wormbase import wormbase_crawler,wormbase_crawler_cds
import pandas as pd 
import numpy as np 
from tqdm import tqdm, trange
import time
import argparse
import os 
'''
step 2.
This code is use for crawl over the wormbase by given .csv file according to different transcript type,
and will save each transcript data as .json format in local, under {version}/{transcript name}.json

input
    -- path of .csv file(with column of transcript)
    -- version of wormabase 
output
    -- {version}/{transcript name}.json
'''
parser = argparse.ArgumentParser(
        prog="three_region.py",
        description="This code is use for crawl over the wormbase by given .csv file according to different transcript type,and will save each transcript data as .json format in local, under {version}/{transcript name}.json",
        epilog="This code can be seen as step2, use directly after worm_mRNA.py"
)
parser.add_argument("arg1",help="This is first argument,for the path that include transcript id as column(.csv file)")
parser.add_argument("arg2",help="This is second argument,for the version of wormbase")
args = parser.parse_args()

#Part 1
print('input file:',args.arg1)
TargetRNA = pd.read_csv(str(args.arg1),encoding = "ISO-8859-1")
Introns = []
TranscriptID = []
Error = []
RNALength = []
FivePrime = []
ThreePrime = []
CDSLength = []
CrawlerCount = 0
# 紀錄已經爬到哪筆資料，避免段網路時需要重新爬
if os.path.exists('./log') == False:
    os.mkdir('./log')
    
with open(f'./log/RecordCrawl_{str(args.arg2)}.txt','a+') as f:
    f.write("\ntranscript")

#Part 2 
for i in trange(len(TargetRNA)):
    ThreePrimeNumber = 0
    FivePrimeNumber = 0
    try:
        ExonIntron,SequenceLength= wormbase_crawler(transcript= TargetRNA['transcript'][i],version=str(args.arg2)) #爬回來的資料會先儲存成df形式
    except:
        try:
            ExonIntron,SequenceLength= wormbase_crawler_cds(transcript= TargetRNA['transcript'][i],version=str(args.arg2))
        except:
            with open(f'./log/Error_new_{str(args.arg2)}.txt','a+') as f :
                f.write(f"\n{TargetRNA['transcript'][i]}")
            Error.append(TargetRNA['transcript'][i])
            continue
    # # section 1
    # try:
    #     FivePrime_df = ExonIntron[ExonIntron['type']=='five_prime_UTR']
    #     for j in range(len(FivePrime_df)):
    #         FivePrimeNumber += FivePrime_df.iloc[j]['stop'] - FivePrime_df.iloc[j]['start'] + 1
    # except:
    #     FivePrimeNumber = 0
    # # section 2
    # try :
    #     ThreePrime_df = ExonIntron[ExonIntron['type']=='three_prime_UTR']
    #     for j in range(len(ThreePrime_df)):
    #         ThreePrimeNumber += ThreePrime_df.iloc[j]['stop'] - ThreePrime_df.iloc[j]['start'] + 1
    # except:
    #     ThreePrimeNumber = 0 
    # # section 3
    # if TargetRNA["type"][i] == 'Coding transcript':
    #     CDSLengthNumber = SequenceLength - FivePrimeNumber - ThreePrimeNumber
    # else:
    #     CDSLengthNumber = 0

    # FivePrime.append(FivePrimeNumber)
    # ThreePrime.append(ThreePrimeNumber)
    # RNALength.append(SequenceLength)
    # CDSLength.append(CDSLengthNumber)
    # TranscriptID.append(TargetRNA['transcript'][i])
    # section 4
    with open(f'./log/RecordCrawl_{str(args.arg2)}.txt','a+') as f:
        f.write(f'\n{TargetRNA["transcript"][i]}')
    CrawlerCount += 1
    #單純輸出爬行到哪一筆的中繼黨(以千筆為單位)
    if CrawlerCount % 1000 == 0:
        with open(f'./log/AlreadyCrawl_{str(args.arg2)}.txt','w') as f :
            for item in TranscriptID:
                f.write('{}\n'.format(TargetRNA['transcript'][i]))
    else:
        pass
    time.sleep(0.3)
# #Part 3 
# ExonIntron_df = pd.DataFrame(columns=['transcript','RNA_length',"5'UTR_length","CDS_length","3'UTR_length"])
# ExonIntron_df['transcript'] = TranscriptID
# ExonIntron_df['RNA_length'] = RNALength
# ExonIntron_df["5'UTR_length"] = FivePrime
# ExonIntron_df["CDS_length"] = CDSLength
# ExonIntron_df["3'UTR_length"] = ThreePrime
# ExonIntron_df.to_csv(args.arg2,index=False)
# with open('Error_new.txt','w') as f:
#     for item in Error:
#         f.write('{}\n'.format(item))


