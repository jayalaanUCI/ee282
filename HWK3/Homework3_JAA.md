# Summary and Data Processing
Made a classrepo directory for Homework3 usign the Createproject script so it mirrors what we have done in classrepo.

```
mamba activate ee282
cd ~myrepos/HWK3


wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/dmel-all-chromosome-r6.48.fasta.gz

```

## Integrity
To check integrity I can download the md5sum text file <br> and simply run the md5sum check command.
```

wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/fasta/md5sum.txt

md5sum -c md5sum.txt
```

This command outputs all the files within the checksum. <br> Only the chromosome file is reads ok which suggest the integrity is good.

## Calculate Summaries of the Genome

```
faSize  dmel-all-chromosome-r6.48.fasta.gz  > data/processed/r6chrom.summary.txt
```

This data outputs a summary (143726002 bases, 1870 sequences, 1152978 N's)

# Summarize Annotation Files
First I need to downloadd the GFF files into its own directory within the raw directory

```
cd ~/myrepos/HWK3/

mkdir GTF

wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/dmel-all-r6.48.gtf.gz
```
## Integrity

Download md5sum file and chemsum integrity

```
wget http://ftp.flybase.net/releases/FB2022_05/dmel_r6.48/gtf/md5sum.txt

md5sum -c md5sum.txt
```

## GFF Report

### Counts of Unique Features 
To output a count of each feature, I used $3 to print the feature column, sort it by reverse numerical order and then count using uniq

```
bioawk -c gff '{print $3}' dmel-all-r6.48.gtf.gz | sort -rn | uniq -c | sort -rn > r6gtf.count.txt
```

This Outputs a txt file with counts for each feature.

### Counts by chromosome

Using a similar ppraoch as before I can ourput the chromosome region, sort it, look for the chromosome patters specifically an duse uniq -c <b/> to output a text file with the per chromosome count.

```
bioawk -c gff '{print $1}' dmel-all-r6.48.gtf.gz | sort -rn | grep -e 'X' -e 'Y' -e '2L' -e '2R' -e '3L' -e '3R' |  uniq -c  > r6gtf.chrom.count.txt
```



