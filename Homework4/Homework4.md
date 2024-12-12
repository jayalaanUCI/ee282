# Homework 4
### Julio Ayala Angulo


Drosophila Melanogaster FASTA file data from Flybase [dmel_r6.60] (http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/dmel-all-chromosome-r6.60.fasta.gz)
```
r6.60="http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/dmel-all-chromosome-r6.60.fasta.gz"
```
>corresponding [md5] (http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/md5sum.txt)

## Summarize Partitions of Genome Assembly

I first need to download the code from the flybase into my raw data folder.

```
cd ~/myrepos/ee282/Homework4/data/raw/

r6.60= "http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.60_FB2024_05/fasta/dmel-all-chromosome-r6.60.fasta.gz"

wget -O r.r60

```

Now I should be able to analyze the data and begin processing. I'll move to my main Homework4 directory to start summarizing.
Usign faFilter, I'll make 2 sets of of data. <= 100kb and >than 100kb. Then I can use faSize to summarize my data for each one before concatenating the files
together.

```
cd ~/myrepos/ee282/Homework4

faFilter -maxSize=100000 ~/myrepos/ee282/Homework4/data/raw/dmel-all-chromosome-r6.60.fasta.gz ~/myrepos/ee282/Homework4/data/processed/r6_max100kb.txt

faFilter -minSize=100001 ~/myrepos/ee282/Homework4/data/raw/dmel-all-chromosome-r6.60.fasta.gz ~/myrepos/ee282/Homework4/data/processed/r6_min100kb.txt
```

These two command produce a list of sequences below 100kb and one above 100kb. I can then summarize each separately and cat their summaries together. Next is to sort them by sequence size.

```
#in the same directory as above

faSize -detailed ~/myrepos/ee282/data/processed/r6_max100kb.txt ~/myrepos/ee282/data/processed/r6_max_summ.text

faSize -detailed ~/myrepos/ee282/data/processed/r6_min100kb.txt ~/myrepos/ee282/data/processed/r6_min_summ.text

cat ~/myrepos/ee282/Homework4/data/processed/r6_max_sum.txt ~/myrepos/ee282/Homework4/data/processed/r6_min_sum.txt > ~/myrepos/ee282/Homework4/data/processed/r6SUM.text

```

Output text file contains separated summary reports of sizes above and below 100kb seperately. 

## Plots for Sequences <= or > than 100kb
### Producing GC content 
I will be using bioawk to process the fasta data seperated by faFilter. First I need to sort seqs by size in reverse numberical order so larger seqs are on top.

Then I generate a report of the sequence and gc content of each seq along with its name and length then sorted by reverse numberical order of the seq size

```
bioawk -c fastx '{print $name, length($seq), gc($seq)}' ~/myrepos/ee282/Homework4/data/processed/r6_max100kb.txt | sort -k2,2 | > ~/myrepos/ee282/Homework4/data/processed/r6MaxGC.txt

bioawk -c fastx '{print $name, length($seq), gc($seq)}' ~/myrepos/ee282/Homework4/data/processed/r6_min100kb.txt | sort -k2,2 | > ~/myrepos/ee282/Homework4/data/processed/r6MinGC.txt

```
Now, I just have to throw the text files into ggplot in R. 
### Sequence Length Distribution
```
#in R
>p <- ggplot(data=r6mingc)
>p+geom_histogram(mapping=aes(x=length), bins=10) + scale_x_log() + ggtitle("R6<100kbp Length Distribution Plot")

>p <- ggplot(data=r6maxgc)
> p+geom_histogram(mapping=aes(x=length), bins=50) + scale_x_log10() + ggtitle("R6>100kbp Length Distribution Plot")

#plots save in ./Homework4/outputs/figures
```
### GC plots
```
# in Rstudio

r6maxgc <- read.delim("C:/Users/julio/OneDrive/Desktop/class/r6maxgc.txt", header=FALSE)
>   View(r6maxgc)
> colnames(r6maxgc) <- c("seq","length","GC")
> r6maxgc$GCBin <- cut(x=r6maxgc$GC, breaks=20)
> p <- ggplot(data=r6maxgc)
> p+geom_bar(mapping=aes(x=GC))
> p+geom_histogram(mapping=aes(x=GC), bin=100)

r6mingc <- read.delim("C:/Users/julio/OneDrive/Desktop/class/r6mingc.txt", header=FALSE)
>   View(r6mingc)
> colnames(r6mingc) <- c("seq","length","GC")
> p <- ggplot(data=r6mingc)
> p+geom_histogram(mapping=aes(x=GC), bins=100)
> p+geom_histogram(mapping=aes(x=GC), bins=10)
> p+geom_histogram(mapping=aes(x=GC), bins=20)
> p+geom_histogram(mapping=aes(x=GC), bins=30
+ )

```
### Culmulative Sequence Plot
Returning to bash
```
sort -rnk2,2 ~/myrepos/ee282/Homework4/data/processed/r6_max_summ.txt > data/processed/maxsorted.txt
sort -rnk2,2 ~/myrepos/ee282/Homework4/data/processed/r6_min_summ.txt > data/processed/minsorted.txt

plotCDF data/processed/maxsorted.txt output/figures/r6maxCDF.png
plotCDF data/processed/minsorted.txt output/figures/r6minCDF.png
```
cd..
# Genome Assembly

In Homework4 directory, I asked for more cores to run. I then converted the hifiasm file to FASTA
### Hifiasm assembly and N50
```
srun -A ee_282 -c 16 --pty /bash/bin -i

hifiasm -o data/processed/ISO1.hifiasm.asm -t 16 data/raw/ISO1_Hifi_AdaptorRem.40X.fasta.gz

awk '/^S/{print ">"$2;print $3}' data/processed/ISO1.hifiasm.asm.bp.hap1.p_ctg.gfa > data/proceessed/ISO1_ctg.fasta

faSize -detailed data/processed/ISO_ctg.fasta | sort -rnk 2,2 | awk '{tot=$2+tot ; print $2 "\t" $1 "\t" tot} END {print tot}' | sort -k1,1rn | awk 'NR==1{tot=$1} NR>1 {print $0 "\t" $3/tot}' > data/processed/sortedContigs.text

cut -f 1 data/processed/sortedContigs.text > data.processed/r6.onlysortedContigs
```

This outputs a file ending showing N50= 23014838 bp or 23Mbp as the the sum of the contigs up to 23014838bp contig covers ~55.5% of the genome. The genome 6 release reads [25.3Mbp](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001215.4/) as the N50.

### Flybase

This code is taken fromt the Pipelines Module as is sorts out the contigs and scaffolds from the Flybase dataset.

```
r6url= "http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz"

wget -O data/raw/r6.scaff.fa.gz -q $r6url

faSize -detailed data/raw/r6.scaff.fa.gz > data/processed/r6.scaff.unsorted.namesizes.txt

sort -rnk 2,2 data/processed/r6.scaff.unsorted.namesizes.txt > data/processed/r6.scaff.sorted.namesizes.txt

cut -f 2 data/processed/r6.scaff.sorted.namesizes.txt > data/processed/r6.scaff.sorted.sizes.txt

faSplitByN data/raw/r6.scaff.fa.gz data/raw/r6.ctg.fa 10

gzip data/raw/r6.ctg.fa

faSize -detailed data/raw/r6.ctg.fa.gz > data/processed/r6.ctg.unsorted.namesizes.txt

sort -rnk 2,2 data/processed/r6.ctg.unsorted.namesizes.txt > data/processed/r6.ctg.sorted.namesizes.txt

cut -f 2 data/processed/r6.ctg.sorted.namesizes.txt > data/processed/r6.ctg.sorted.sizes.txt
```

Now just use PlotCDF to produce the plots.
```
plotCDF data/processed/r6.onlysortedContigs.txt output/figures/hifiasmCDF.png
plotCDF data/processed/r6.ctg.sorted.sizes.txt output/figures/flybaseCDF.png
plotCDF data/processed/r6.scaff.sorted.sizes.txt  output/figures/scaffCDF.png
```

From this we see that the scaff data has a tighter distribution than the hifiasm or the flybase contigs. Hifiasm has a a bit of a peak as it begins to plateu.

### Busco Scores
Downloaded BUSCO on seperate env and will need to compare datasets to the diptera_odb10 dataset.
```
gunzip -c  data/raw/r6.scaff.fa.gz > data/raw/r6.scaff.fa
busco -i data/raw/r6.scaff.fa -m genome -c 16 -l diptera_odb10 -o r6.ctg

gunzip -c  data/raw/r6.ctg.fa.gz > data/raw/r6.ctg.fa
busco -i data/raw/r6.ctg.fa -m genome -c 16 -l diptera_odb10 -o r6.ctg


busco -i data/processed/ISO_ctg.fasta -m genome -c 16 -l diptera_odb10 -o r6.hifiasm

```
| Dataset    |  C    |  S  |  D |  F | M |  Total |
| --- | ----- | ---- | --- | --- | --- | ------ |
| r6.contig | 3276 | 3267 | 9 | 5 | 4 | 3285 |
|r6.scaff   | 3276| 3267| 9 | 5 | 4 | 3285 |
|r6.Hifiasm | 2777 | 2767 | 10 | 7 | 501 | 3285 |

#Extra Credit

Used Browser option for DGenie and downloaded fasta.gz files for the r6.contig and the r6.hifiasm datasets and ran them throught the aligner.
 Saved Map png and summary data.