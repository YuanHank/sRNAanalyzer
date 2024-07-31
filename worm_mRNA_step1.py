#### make worm_mRNA.fasta
# mRNA is ATCG while downloaded from wormbase
# ncRNA is AUCG 
"""
step 1.
This code is to transfer the .fa file from wormbase to .csv file with transcript as column
"""
import argparse
import pandas as pd
parser = argparse.ArgumentParser()
parser.add_argument('--input',help='This is for the mRNA input file that downloaded from wormbase(e.g. c_elegans.PRJNA13758.WS285.mRNA_transcripts.fa)')
parser.add_argument('--version',help='This is for the wormbase version')
parser.add_argument('--type',help="This is for the RNA type, including CDS, mRNA (depend on the input file)")
args = parser.parse_args()
seq=""
read = False
namelist=[]
with open(f"{args.input}") as f:
    with open(f"worm_{args.type}_{args.version}.fasta","w") as of:
        for row in f:
            if row.startswith(">"):
                rowsplit = row.split(" ")
                if seq != "":
                    of.write(seq+'\n')
                namelist.append(rowsplit[0].strip(">"))
                of.write(rowsplit[0]+'\n')
                seq=""
            else:
                seq+=row.strip("\n")
        else:
            of.write(seq+'\n')

df = pd.DataFrame()
df["transcript"] = namelist
df.to_csv(f"WS{args.version}_{args.type}.csv", index=False)
