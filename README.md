## flow 
使用檔案:
wormbase下載指定版本的mRNA or cds檔案
e.g. 
  * c_elegans.PRJNA13758.WS290.mRNA_transcripts.fa
  * c_elegans.PRJNA13758.WS290.CDS_transcripts.fa

### step 1. 整理檔案成.fasta格式與.csv格式
python worm_mRNA_step1.py --input {輸入檔案, 範例:c_elegans.PRJNA13758.WS285.mRNA_transcripts.fa} --version {wormbase資料版本, 範例:285} --type {檔案為mRNA or CDS, 範例:mRNA}

output : WS285_mRNA.csv, worm_mRNA_285.fasta
### step 2. 爬蟲(整理transcript成為.json並存放於資料夾,使用到wormbase.py)
python WS285_mRNA.csv 285

output : mRNA_json_285/

### step 3. 整理成輸出reference前中繼檔
python intron_analysis_step3.py --input {存放json檔案的位置,範例:mRNA_json_285/} --input2 {含有提供CDS transcript的檔案,範例:WS285_cds.csv} --output {輸出位置,範例:./} --file_vesrion {使用檔案的版本,範例:285} --length_filter {長度篩選,範例:100}

output : intron_analysis_285.csv

### step 4. 由步驟3.的中繼檔輸出最終reference(mRNA_sRNAanalyzer_ref.py & len2bin.py)
#### 出3區域,長度,IDtoName,metagene
python mRNA_sRNAanalyzer_ref.py --input {步驟3.的中繼檔,範例intron_analysis_285.csv} --file_version {wormbase檔案使用的版本,範例:285} --output {輸出位置,範例:./} --type {要輸出的類型,範例:metagene} --mRNA_file {轉換mRNA與gene的檔案，非必要輸入,範例:c_elegans.PRJNA13758.WS290.mRNA_transcripts.fa}
#### 為了mRNA_WS285_length.csv
python len2bin.py --input {步驟3.的中繼檔,範例intron_analysis_285.csv} --output {輸出檔名} -bin {bin的數量,100}



