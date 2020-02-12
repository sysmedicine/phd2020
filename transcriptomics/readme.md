##  Transcriptomics Workshop

In this workshop, we will do several steps to analyze transcriptomics data. We will start with data quantification, followed by differential expression and functional analysis. As for the data, we will be using publicly available data [PRJNA412001](https://www.ebi.ac.uk/ena/data/view/PRJNA412001).

Related Publication: [HIF1 Mediates a Switch in Pyruvate Kinase Isoforms After Myocardial Infarction](https://pubmed.ncbi.nlm.nih.gov/29652636-hif1-mediates-a-switch-in-pyruvate-kinase-isoforms-after-myocardial-infarction/)

### Quantification of the data

To quantify the data, we will use a software called Kallisto.

#### Install Kallisto

Kallisto has been installed in your pre-course preparation step. If you want to install it, please refer to [their website](https://pachterlab.github.io/kallisto/about)

#### Kallisto for FASTQ Quantifications

Our goal in this step is to translate the FASTQ files from RNA-seq analysis to understandable quantified numbers (raw count and TPM).

1. Download several files:

  * the fastq files from the [data repository](https://www.ebi.ac.uk/ena/data/view/PRJNA412001). **Note:** For the sake of the time, just download a condition (2 fastq files --> paired).

  * cDNA file from Ensembl website for Mouse (Mus Musculus). Look for "Mus_musculus.GRCm38.cdna.all.fa.gz"
  
  ![cdna](img/cdna_dl.png)
  
  * Mapping of transcript to gene name in [Biomart](https://www.ensembl.org/biomart/martview). Select "Ensembl Gene 99" --> "Mouse genes (GRCm38.p6)". In the "Filter" tab (left), select "Gene", then select "protein_coding" for both "Gene Type" and "Transcript Type". Next, select "Attributes" tab (left), select "Gene" --> select only "Transcript Stable ID Version", "Gene Stable ID", and "Gene Name". You can select whatever you like, on top of those three features. Click on "Count" (make sure you have "22191/56289 Genes") and then "Result" (above), download the Result as TSV.
  
  ![biomart](img/biomart.png)

2. Create Index File. Index file only needs to be created once, unless you want to change the reference genome/species

```shell
$ kallisto index -i [AnyIndexName] /your_path_to/Mus_musculus.GRCm38.cdna.all.fa.gz
```

3. Quantify the FASTQ files. You can process a paired and replicated fastq files together, but make sure that it is belong to one sample. This might take a while

```shell
$ kallisto quant -i [AnyIndexName] -o /path_where_you_want/to_store_your_result/ /your_path_to/fastq1.tar.gz /your_path_to/fastq2.tar.gz
```

4. Once the process is completed, in your output directory, you will have several files: "abundance.h5", "abundance.tsv", "run_info.json". Open "abundance.tsv" in excel, and see what data that you get from Kallisto (hint: look at the column name)

#### Consolidating multiple Kallisto samples

**Goal:** Converting single sample tsv into two big tables with rownames of gene name and column name of the value (one for raw count, one for TPM).

1. Download the Kallisto result from all 8 samples in [this link](data/KallistoResults.zip). Unzip it into a folder.

2. Open Rstudio and run this script. Don't forget to enter the path to the unzipped folder. This script is combining multiple samples into a big dataframe of count and TPM values. Remember, this is still in transcript level.

```R
result_dir = 'data/' # The directory that includes the extracted folders

count_tr = data.frame()
tpm_tr = data.frame()

for(i in list.dirs(result_dir,recursive=F)){
  temp = read.csv(paste0(i,'/abundance.tsv'),sep='\t',stringsAsFactors = F)
  temp_count = data.frame(temp$est_counts)
  temp_tpm = data.frame(temp$tpm)
  colnames(temp_count) = gsub("data//", "", i)
  colnames(temp_tpm) = gsub("data//", "", i)
  if(ncol(count_tr) == 0){
    count_tr = temp_count
    rownames(count_tr) = temp$target_id
    tpm_tr = temp_tpm
    rownames(tpm_tr) = temp$target_id
  } else {
    count_tr = cbind(count_tr, temp_count)
    tpm_tr = cbind(tpm_tr,temp_tpm)
  }
}
```

3. Combining the transcripts into the gene names: Merge the transcripts into gene names. Remember that earlier we chose only protein_coding gene and transcripts mapping from Biomart.

```R
biomart_file = 'mart_export.txt' #adjust this based on your file location
mapping = read.csv(biomart_file,sep='\t',stringsAsFactors = F,row.names = 1)
count_gn = merge(count_tr,mapping['Gene.name'], by=0,all=F) # merge only for the shared row names
count_gn = count_gn[,2:ncol(count_gn)]
count_gn = aggregate(.~Gene.name, count_gn, sum)
rownames(count_gn) = count_gn$Gene.name 
tpm_gn = merge(tpm_tr,mapping['Gene.name'], by=0,all=F) # merge only for the shared row names
tpm_gn = tpm_gn[,2:ncol(tpm_gn)]
tpm_gn = aggregate(.~Gene.name, tpm_gn, sum)
rownames(tpm_gn) = tpm_gn$Gene.name
```

**Question: How many protein coding transcripts and genes that we have?**

4. Save the count and TPM table as tab-separated format.

```R
write.table(count_gn,file='/path_to_save/count_gene.txt',sep = '\t', na = '',row.names = F)
write.table(tpm_gn,file='/path_to_save/tpm_gene.txt',sep = '\t', na = '',row.names = F)
```
### Data Exploration

Once we have the consolidated table, we have to look at the data to understand it better. In this step, we will do a simple PCA analysis to see the data separation of each condition.

1. Download the [metadata file](data/metadata.txt) and load it together with TPM table that we generated to R. 

**Question: Why TPM and not raw count?**

```R
metadata_file = 'metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

tpm_file = 'tpm_gene.txt' #adjust this based on your file location
tpm = read.csv(tpm_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)] # make sure that sequence of metadata is the same with tpm
```

2. Perform the PCA with the script below.

```R
PCA=prcomp(t(tpm), scale=F)
plot(PCA$x,pch = 15,col=c('blue','blue','red','red','lightgreen','lightgreen','black','black'))
```

**Question: What do you think of the sample separations?**


### Differential Expression Analysis

To perform differential expression analysis, we will use a package in R called [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html).

1. Load the count and metadata to R
**Question: Why count and not TPM?**

```R
metadata_file = 'metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

count_file = 'count_gene.txt' #adjust this based on your file location
count = read.csv(count_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)] # make sure that sequence of metadata is the same with tpm
```

2. Perform the Differential Expression Analysis

```R
library('DESeq2')
metadata=as.factor(metadata$condition)
coldata <- data.frame(row.names=rownames(metadata),metadata)
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(count)),colData=coldata,design=~conds)
dds <- DESeq(dds)
```

3. Retrieve the specific comparison results. For now, retrive comparison between MI_1D vs SHAM_1D.

```R
cond1 = 'MI_1D' #First Condition
cond2 = 'SHAM_1D' #Reference Condition
res=results(dds,contrast=c('conds',cond1,cond2))
res=data.frame(res)
```

4. Save it as TSV file

```R
write.table(res,file='/path_to_save/deseq_1D.txt',sep = '\t', na = '',row.names = T,col.names=NA)
```

5. Retrive and save comparison between "MI_3D" vs "SHAM_3D" as well

6. Open the saved result for "MI_1D" vs "SHAM_1D".

**Question 1: With Adjusted P-value < 0.05, how many genes are significantly differentially up-regulated? down-regulated? (HINT: Look at Log2FoldChange) What is the most affected gene? (HINT: sort it based on Adjusted P-Value)**
**Question 2: Why do we have several genes with NA in their statistical result?**

![deseq](img/deseq.png)

### Functional Analysis

The goal of this analysis is to map the changes in the transcriptome level to functional level (Gene Ontology, KEGG Pathways, etc). We will use a tool in R called [PIANO](https://bioconductor.org/packages/release/bioc/html/piano.html)

1. Load the differential expression result for "MI_1D" vs "SHAM_1D"

```R
deseq_file='/path_to_saved/deseq_1D.txt'
deseq = read.csv(deseq_file,sep='\t',stringsAsFactors = F,row.names = 1)
```

2. Load the KEGG Pathway gene-set collection. It can be downloaded from [this link](data/KEGG.gmt)

```R
library('piano')
library('Biobase')
library('snow')
library('RColorBrewer')
library('gplots')
library('visNetwork')

GSC='/path_to_downloaded/KEGG.gmt'
y=loadGSC(GSC)
```

3. Perform the functional analysis with Piano

```R
input_file=deseq[ ,c('log2FoldChange','pvalue')]
logFC=as.matrix(input_file[,1])
pval=as.matrix(input_file[,2])
rownames(logFC)=toupper(rownames(input_file))
rownames(pval)=toupper(rownames(input_file))
logFC[is.na(logFC)] <- 0
pval[is.na(pval)] <- 1
gsaRes <- runGSA(pval,logFC,gsc=y, geneSetStat="reporter", signifMethod="nullDist", nPerm=1000)
```

4. Retrieve and save the result of KEGG functional analysis as TSV

```R
res_piano=GSAsummaryTable(gsaRes)
```

5. Visualize the result as heatmap and save it as PDF (you'll probably see the "Error in plot.new() : figure margins too large", ignore it for now)
```R
pdf("heatmap.pdf") 
hm = GSAheatmap(gsaRes, adjusted = T)
dev.off()
```

5. Visualize the result as network of terms and save it as PDF 
```R
pdf("/path_to_save/network_plot.pdf") 
nw = networkPlot(gsaRes, class="distinct", direction="both",significance=0.00005, label="names")
dev.off() 
```

6. Generate Interactive Network and save as HTML 
```R
nw_int = networkPlot2(gsaRes, class="distinct", direction="both", significance=0.00005)
#if you want to show it without saving in R, type "nw_int"
visSave(nw_int, file = "network_plot_interactive.html", background = "white")
```

**NOTE: Adjust the significance level during network plotting if necessary (you can increase or decrease it) if necessary**

7. Redo the same analysis for GO Biological Process. The gene-set collection can be downloaded from [this link](data/GO.gmt). 

8. Redo all the steps above for "MI_3D" vs "SHAM_3D"

**Question 1: How many KEGG pathways are up-regulated and down-regulated in 1D? (HINT: look at distinct adjusted p-value up and down)**

**Question 2: How many GO Biological Processes are up-regulated and down-regulated in 1D? (HINT: look at distinct adjusted p-value up and down)**

**Question 3: How many KEGG pathways are up-regulated and down-regulated in 3D? (HINT: look at distinct adjusted p-value up and down)**

**Question 4: How many GO Biological Processes are up-regulated and down-regulated in 1D? (HINT: look at distinct adjusted p-value up and down)**

**Question 5: Is it getting worse or better with time? Also which processes and pathways are specific for 1D and for 3D?**

**Question 6: What can you summarize from this?**

### Multiple Factor in Transcriptomic Analysis.

Until now, the comparison that we are doing is a simple comparison, as in A vs B. In real-life, you might encounter a condition where we have multiple factors. Using the same example, what if we want to compare overall MI vs SHAM, and correcting it based on the day? Another example, if you have several batches with a clear batch effect, so we have to find a way to correct for it.

The goal of this lab is how to consider multiple factor in transcriptomic analysis. The only difference with the previous lab is on how we approach the problem in DESeq2.

1. Load the count and metadata to R (download new metadata [here](data/metadata_multi.txt))

```R
metadata_file = 'metadata_multi.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

count_file = 'count_gene.txt' #adjust this based on your file location
count = read.csv(count_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)] # make sure that sequence of metadata is the same with tpm
```

2. Create Model Matrix

```R
library('DESeq2')
model=model.matrix(~MI+day,metadata)
```

3. Perform Differential Expression analysis

**Question: What is the difference between this and the previous deseq step?**

```R
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(count)),colData=coldata,design=~conds)
dds <- DESeq(dds)
```

3. Retrieve the specific comparison results. For now, retrive comparison between MI_1D vs SHAM_1D.

```R
cond = 'MI1' #First Condition
res=results(dds,contrast=list(c(cond)))
res=data.frame(res)
```

4. Save it as TSV file

```R
write.table(res,file='/path_to_save/deseq_multifactor.txt',sep = '\t', na = '',row.names = T,col.names=NA)
```

5. Perform the Functional Analysis on this Multifactor analysis and compare the result with single factor analysis.

