##  Study: Curcumin Treatment of AOM-DSS Model
### Scientist: Yue, Renyi Wu
### Data Analysis: Renyi Wu, Davit Sargsyan 
### Created: 11/06/2017 

---

## Daily Logs
### 03/24/2018
* Added script to redraw figures 3D and 5A for the paper(*aom_dss_cur_figures_pub_v1.R*)    
* Added new documents (not synced to GitHub, see CVI computer): Renyi's request and publicaion proof    

###  11/27/2017
* Organized the repository    

## Stored Files (4Tb internal drive)
### File Legend
[sample list NGS AOMDSS_Yue.xlsx](/home/administrator/Documents/aom.dss.cur/docs)

### Original Data (FastQ) Location
1. [First 8 samples (Pool1, C1 through C70)](/media/administrator/datastorage/2016_Isilon/Kong/MethylSeq/FASTQ_Kong_MethylSeq_Pool1)    
2. [Next 8 samples(Pool2, C14 through C79)](/media/administrator/datastorage/2016_Isilon/Kong/MethylSeq/FASTQ_Kong_MethylSeq_Pool2)    
3. [Next 6 samples("Christmas Samples", C15 thorugh C59). NOTE: note good.](/media/administrator/datastorage/2016_Isilon/Kong/MethylSeq/FastQ Files_repeated_yue)
4. [Rerun of the last batch of 6 samples (April 2017, C15 through C59). NOTE: not good.](/media/administrator/datastorage/FastQ_2017/Methyl_seq/April)    
5. [Rerun of the last batch of 6 samples (October 2017, C15 through C59)](/media/administrator/datastorage/FastQ_2017/Methyl_seq/October)    
NOTE1: 6 samples from #1 & #2 above were used for the publication: "Curcumin AOM-DSS Mouse at 18 Weeks" (not published yet as of 11/27/2017)    
NOTE2: additional to samples in NOTE1 above, 6 samples from #4 were used for RO1 grant in September 2017    

### Processed BAM files
1. BAM files from batches #1 and #2 processed by John are [here](/media/administrator/datastorage/2016_Isilon/john/MethylSeq).   
2. BAM files from batches #1, #2, #3 and #5 processed by Davit are [here](/media/administrator/datastorage/Processed_BAM_Files/Yue_AOM_DSS_Cur_MethylSeq_Processed/Davit)    
3. BAM files from batch #5 above processed by Renyi are [here](/media/administrator/datastorage/Processed_BAM_Files/Yue_AOM_DSS_Cur_MethylSeq_Processed/Renyi)

## ToDo 
### 11/27/2017
* Install **formatR**, **rprojroot** and **rmarkdown** packages.
 
 
M=\log _{2}(R/G)=\log _{2}(R)-\log _{2}(G)
{\displaystyle A={\frac {1}{2}}\log _{2}(RG)={\frac {1}{2}}(\log _{2}(R)+\log _{2}(G))} A={\frac  12}\log _{2}(RG)={\frac  12}(\log _{2}(R)+\log _{2}(G))

## Data Files Original Locations
1. comb2.csv: first batch
/media/administrator/datastorage/John/MethylSeq

2. combined_yue_2nd_dedup.csv: before Christmas 2017 batch
/home/administrator/Documents/Renyi/yue_cur_methyl_2nd

3. Methyl_rep_2017.samdup.csv: October 2017 rerun
/home/administrator/Documents/Methyl_Oct2017

## Reference
MA Plot:
[https://en.wikipedia.org/wiki/MA_plot]

Volcano Plot:
[https://en.wikipedia.org/wiki/Volcano_plot_(statistics)]

### 11/09/2017
* Added annotation and plotted TNFa methylation by DMR in the 2 batches    

## Samples Legend
### Batch 1 8-week
ID	Sample name
C1	Negative Control 8wk
C2	Negative Control 8wk
C26	AOM+DSS 8wk
C29	AOM+DSS 8wk
C42	AOM+DSS+Cur 8wk
C47	AOM+DSS+Cur 8wk
C65	DSS 8wk
C70	DSS 8wk
C75	DSS+Cur 8wk
C79	DSS+Cur 8wk

### Batch 1 18-week
ID	Sample name
C14	Negative Control 18 wk
C20	Negative Control 18 wk
C34	AOM+DSS 18wk
C40	AOM+DSS 18wk
C54	AOM+DSS+Cur 18wk
C60	AOM+DSS+Cur 18wk

### Batch 2 Dec 2016 and Oct 2017 Repeat
ID	Sample name
C15	Negative Control 18wk
C19	Negative Control 18wk
C33	AOM/DSS 18wk
C36	AOM/DSS 18wk
C55	AOM/DSS+Cur 18wk
C59	AOM/DSS+Cur 18wk

Install!
formatR, rprojroot, rmarkdown. 
 
 
M=\log _{2}(R/G)=\log _{2}(R)-\log _{2}(G)
{\displaystyle A={\frac {1}{2}}\log _{2}(RG)={\frac {1}{2}}(\log _{2}(R)+\log _{2}(G))} A={\frac  12}\log _{2}(RG)={\frac  12}(\log _{2}(R)+\log _{2}(G))

## Data Files Original Locations
1. comb2.csv: first batch
/media/administrator/datastorage/John/MethylSeq

2. combined_yue_2nd_dedup.csv: before Christmas 2017 batch
/home/administrator/Documents/Renyi/yue_cur_methyl_2nd

3. Methyl_rep_2017.samdup.csv: October 2017 rerun
/home/administrator/Documents/Methyl_Oct2017

## Reference
MA Plot:
[https://en.wikipedia.org/wiki/MA_plot]

Volcano Plot:
[https://en.wikipedia.org/wiki/Volcano_plot_(statistics)]