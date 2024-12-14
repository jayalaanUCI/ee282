# Final Project: Determining site of insertion of hCD1a transgene in CD1a mouse model
### Julio Ayala Angulo

## Introduction
Type 2 Diabetes (T2D) is the 8th leading cause of death, with cases expected to reach 1.3 billion globally by 2050 (Ong et al. 2023). Current therapies manage symptoms but fail to cure the disease, with efficacy declining after four years (GRADE Study, NEJM, 2022). The lack of a curative treatment stems from an incomplete understanding of T2D's root causes and complete mechanistic model. Research has largely focused on glycemic and insulinemic dysregulation as core drivers. Adiposity, dyslipidemia, and immune inflammation have seen a resurgence of interest but remain understudied partially due to poor predictive and mechanistic mouse models and the difficulty of elucidating pathophysiological mechanisms in humans. Increased Visceral Adipose Tissue (VAT) mass and dyslipidemia in the form of elevated free fatty acids (FFA), low and very low density lipoproteins (LDL/VLDL) and triglycerides (TG) in T2D patient plasma has been well established. The current literature largely focuses on the impact hyperinsulinemia has on adipose mass and dyslipidemia, therefore assuming that the change in lipids is the outcome and that insulin is the cause(Galicia-Garcia et al. 2020). However, there is evidence that the converse is also true, though it has been less studied. For example, plasma FFAs are predictive of hyperglycemia and strongly correlated with insulin resistance (IR)(Fryk et al. 2021)(Mahendran et al. 2013). Higher FFA and glycerol levels in obese patients without T2D compared to obese with T2D suggests a crucial role in lipid levels in the transition from Insulin Resistant/Prediabetic to T2D [Fryk et al. 2021). Onset of T2D starts with increasing obesity and FFA as hyperinsulinemia develops to maintain blood sugar and lipid levels. Eventually, a mixture of insulin resistance and pancreatic beta cell atrophy drive hyperglycemia (Johnson et al. 2021). However, the cause of the atrophy or why adiposity and FFA don't drive hyperglycemia earlier remains unknown. 
CD1a is one of 4 CD1 molecules in humans that allow the immune system to respond to lipid antigens. Mice only have 1 member of the CD1a family: CD1d, which is limited in its antigen presentation to microbial metabolites (Moody and Cotton, 2017). Consequently, mice are largely unable to adaptively respond to lipids. In a disease context where dyslipidemia is not only characteristic, but potentially causative, anti-lipid T cell responses can bridge the gap between mouse models and human disease. Of all CD1 molecules, most is known about CD1a. Interestingly, mice not only do not model T2D and oebsity well, but also do not model the low-chronic inflammation or tissue specif inflammation in humans (Elevated IL17a+ and TNFa+) where as mice only express tissue-specfic inflammation of a different type (IFNy+ and TNFa+) (Wiener et al. 2009). We hypothesized that human CD1a lipid antigen presentation to T cells drives progression of T2D through inflammation. To this extent, we acquired CD1a expressing B6 mice. Under a high fat diet, we saw changes in metabolism and inflammation in CD1a mice. However, it is still unknown if CD1a expression in these mice proproperly mirror humans. The original lab that developed these mice only validated its expression in the skin and thymus as described in humans (Kobayashi et al. 2011). The insertion included regulatory elements for CD1a in humans, but it's insertion is likely not targeted to chr1 as in humans. Initially, we used flow cytometry to characterized CD1a+ cells proportions in mice (Figure 1). This data validates the high expression in tissues previsly seen in humans (Thymus and skin) while also basal expression in liver and pancreas (not described) while notably absent in lymphatic tissues. Additonally, there is a clear lower expression of CD1a in female mice. Depending on the insertion site(s), this could be an artifact of the insertions location as well as its expression in pathologically relevant tissues. 

![CD1a Expression Across tissues](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/CD1a%20Expression.png?raw=true)

### Goals
1. Chromosomes Harboring CD1a transgene
    -Primarily, I wish to determine which chromosomes have copies of CD1a. This lets me know which trangenic mice, like RAG1KO mice,  I am able to breeed with the CD1a mice. It mght also shed light on differeces across sexes if some copies lie on the X or Y chromosomes.
    -I aim to output a table per sequence with chromosome and location afte raligbnemnt that I can graph with ggplot2.
2. Vizualize spread of of insertion across the mouse genome using [IGV](https://igv.org/). 
## Methods and Results
### Sample Preparation and Targeted Locus Amplification
DNA samples for sequencing where generated from CD1a mouse splenocytes, lymph nodes, and blood. Tissues were macerated and red blood cells were lysed with RBC lysis buffer. Mononuclear cells were then put through this [TLA protocol](https://link.springer.com/protocol/10.1007/978-1-4939-6442-0_13) for isolation of DNA. In short, cells were lysed in the presensence of paraformaldehyde to fix DNA to proximal transcription factors, enahncers, or another regulatory proteins. Samples were then digested and subsequebtly ligated to link DNA sequences attached to the same regualtory elements. Proteins were digeted before ligating DNA once more to form circular DNA. CD1a-enriched amplicons were generated through PCR usign inverted genotying primers. DNA amplicons were purified via AMPure XP beads before beign sent for Nanopore Sequencing at [Plasmidsaur](https://plasmidsaurus.com/).
### Required Packages and alnguages
I will be using Mamba as a package manager through which most of these packages were downloaded.
| Packacage/Utility |      Purpose      |
| ----------------- | ----------------- | 
| Samtools       |    managing sam files               |   
| NanoQC            |   Quality checks  |  
| Minimap2          |   Alignement      |   
| Bioawk            | Data handling     |   
| R and ggplot2     |  Graphing         |   
| IGV (online)      | Alignment         |

Data and output will be posted to [github](https://github.com/jayalaanUCI/ee282). 

## Code 
```
# !/bin/bash
#Uploaded raw FastQ 7PCS4R_1_CD1a_TLA_amplicon_primer_3.fastq to ./FinalProject/data/raw/ before chaning name before and compressing
cd ~/myrepos/FinalProject/data/raw
mv 7PCS4R_1_CD1a_TLA_amplicon_primer_3.fastq CD1aTLA.fastq | gzip 
```
Next, my goal is to run NanoQC to out out quality check graphs specific to the nanopore platform.

```
nanoqc -o output/figures/ data/raw/CD1aTLA.fasta.gz

```
![NanoQC plots](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/newplot.png?raw=true)

> Figure 2. NanoQC plots showing size and quality of reads.

Based on the QC data, Average read length is ~400bp  with several long reads and extremly long reads. 

Next is the alignment to a mouse reference genome using Minimap2. Mouse reference genome used is the GCRm39 genome from [Ensembl](https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/).
I also used bioawk to generate as a text file of chromosomes that the sequence aligned by extractting the rname field.

I also made sorted and indexed BAM files for IGV visuzalization.

```
cd ./FinalProject/data/raw/
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

cd ../../

minimap -a data/raw/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz data/raw/CD1aTLA.fasta.gz > CD1aTLA.sam

### Chromosome list
samtools sort data/processed/CD1aTLA.sam -o data/processed/sortedCD1aTLA.sam
bioawk -c sam '{print $rname}' data/processed/sortedCD1aTLA.sam > data/processed/TLAchr.text
### Vizualize files
samtools view -Sb CD1aTLA.sam > CD1aTLA.bam
samtools sort data/processed/CD1aTLA.bam -o data/processed/CD1aTLAsort.bam
samtools index CD1aTLAsort.bam
```
## Results

### Table of chromosomes
In R, usign ggplot2 I took the text file TLAchrt.txt and used the chromosome list to graph the chromosome that the reads aligned to.
```
### in R
>x <- read.table(TLAchr.txt)
>q+geom_density(mapping=aes(x=V1), fill="#69b3a2", color="#e9ecef", alpha=0.8)+xlab("Chromosome"+ylab("Reads Mapped") ### where V1 is the chr column
```
### Alignment chart from R
![Rplot](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/Rplot.png?raw=true)
>Figure 2. Density plots displaying alignment of Reads across Chromosomes # or unmapped * .

Files then were uploaded into IGV with a focus on chromosomes 1 and 10 based on the previsoly generated table.
### For Chromosome 1 alignment.
![fig3.1](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshotchr1.png?raw=true)
![Figure 3](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshot.png?raw=true)
![Figure 4](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshot2.png?raw=true)
![Fgure 5](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshot3.png?raw=true)
![Figure 6](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshot4.png?raw=true)
>Figure 3. IGV Visualization of TLA sequences on Chromosome 1.
### For Chromosome 10
![Figure 7](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshot6.png?raw=true)
![Figure 8](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshot7.png?raw=true)
![Figure 9](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/igv_snapshot8.png?raw=true)
>Figure 4. IGV Visualization of TLA sequences on Chromosome 10.



## Discussion
The data and visuzalizations shows alignement of CD1a TLA sequences to chromosomes 1 and 10 specifically. However, there are several reads that remained unmapped. This is likely due to the amplicons being a mixture of mouse and human elements as the TLA protocol would amplify any proximal DNA to the CD1a primer sites. Its likely at least some of these unaligned regions are due to the presensence of human regulatory sequences inserted with the CD1a transgene. 
Immediate future steps would be to filter out the unmapped sequences by size. If the reads are long, in silico restriction enzyme digestions could be a good options to produce smaller reads that might better align. Additionally, aligning to a human reference genome for Chr1, the original site of CD1a, should show the proportion of reads derived from these regulatory regions.
Furthermore, alignment to enhancer region using the [Enhancer ATLAS](http://www.enhanceratlas.org/downloadv2.php) database for immune cells references might even help in pinpointing specfic enhancers included in the transgene.
Regardless, as it stands, CD1a mice are breedable with RAG1KO mice as the RAG1 gene is on Chr 2. 

## Citations
1. [Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, et al. 2021\. Twelve years of SAMtools and BCFtools. *GigaScience* **10**: giab008.](https://www.zotero.org/google-docs/?Oy4FWS)  
2. [De Coster W, D’Hert S, Schultz DT, Cruts M, Van Broeckhoven C. 2018\. NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics* **34**: 2666–2669.](https://www.zotero.org/google-docs/?Oy4FWS)  
3. [Fryk E, Olausson J, Mossberg K, Strindberg L, Schmelz M, Brogren H, Gan L-M, Piazza S, Provenzani A, Becattini B, et al. 2021\. Hyperinsulinemia and insulin resistance in the obese may develop as part of a homeostatic response to elevated free fatty acids: A mechanistic case-control and a population-based cohort study. *EBioMedicine* **65**: 103264\.](https://www.zotero.org/google-docs/?Oy4FWS)  
4. [Galicia-Garcia U, Benito-Vicente A, Jebari S, Larrea-Sebal A, Siddiqi H, Uribe KB, Ostolaza H, Martín C. 2020\. Pathophysiology of Type 2 Diabetes Mellitus. *Int J Mol Sci* **21**: 6275\.](https://www.zotero.org/google-docs/?Oy4FWS)  
5. [Hottentot QP, van Min M, Splinter E, White SJ. 2017\. Targeted Locus Amplification and Next-Generation Sequencing. In *Genotyping: Methods and Protocols* (eds. S.J. White and S. Cantsilieris), pp. 185–196, Springer, New York, NY https://doi.org/10.1007/978-1-4939-6442-0\_13 (Accessed October 16, 2024).](https://www.zotero.org/google-docs/?Oy4FWS)  
6. [Johnson JD. 2021\. On the causal relationships between hyperinsulinaemia, insulin resistance, obesity and dysglycaemia in type 2 diabetes. *Diabetologia* **64**: 2138–2146.](https://www.zotero.org/google-docs/?Oy4FWS)  
7. [Kobayashi C, Shiina T, Tokioka A, Hattori Y, Komori T, Kobayashi-Miura M, Takizawa T, Takahara K, Inaba K, Inoko H, et al. 2011\. GM-CSF-independent CD1a expression in epidermal Langerhans cells: evidence from human CD1A genome-transgenic mice. *J Invest Dermatol* **132**: 10.1038/jid.2011.280.](https://www.zotero.org/google-docs/?Oy4FWS)  
8. [Li H. 2024\. lh3/bioawk. https://github.com/lh3/bioawk (Accessed December 13, 2024).](https://www.zotero.org/google-docs/?Oy4FWS)  
9. [Li H. 2018\. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* **34**: 3094–3100.](https://www.zotero.org/google-docs/?Oy4FWS)  
10. [Mahendran Y, Cederberg H, Vangipurapu J, Kangas AJ, Soininen P, Kuusisto J, Uusitupa M, Ala-Korpela M, Laakso M. 2013\. Glycerol and fatty acids in serum predict the development of hyperglycemia and type 2 diabetes in Finnish men. *Diabetes Care* **36**: 3732–3738.](https://www.zotero.org/google-docs/?Oy4FWS)  
11. [Moody DB, Cotton RN. 2017\. Four Pathways of CD1 Antigen Presentation to T cells. *Curr Opin Immunol* **46**: 127–133.](https://www.zotero.org/google-docs/?Oy4FWS)  
12. [Ong KL, Stafford LK, McLaughlin SA, Boyko EJ, Vollset SE, Smith AE, Dalton BE, Duprey J, Cruz JA, Hagins H, et al. 2023\. Global, regional, and national burden of diabetes from 1990 to 2021, with projections of prevalence to 2050: a systematic analysis for the Global Burden of Disease Study 2021\. *The Lancet* **402**: 203–234.](https://www.zotero.org/google-docs/?Oy4FWS)  
13. [Robinson JT, Thorvaldsdóttir H, Winckler W, Guttman M, Lander ES, Getz G, Mesirov JP. 2011\. Integrative genomics viewer. *Nat Biotechnol* **29**: 24–26.](https://www.zotero.org/google-docs/?Oy4FWS)  
14. [Wickham H. 2016\. *ggplot2*. Springer International Publishing, Cham http://link.springer.com/10.1007/978-3-319-24277-4 (Accessed December 13, 2024).](https://www.zotero.org/google-docs/?Oy4FWS)  
15. [Winer S, Chan Y, Paltser G, Truong D, Tsui H, Bahrami J, Dorfman R, Wang Y, Zielenski J, Mastronardi F, et al. 2009\. Normalization of obesity-associated insulin resistance through immunotherapy. *Nat Med* **15**: 921–929.](https://www.zotero.org/google-docs/?Oy4FWS)  
16. [R: The R Project for Statistical Computing. https://www.r-project.org/ (Accessed December 13, 2024a).](https://www.zotero.org/google-docs/?Oy4FWS)  
17. [Role of Insulin Resistance in Human Disease | Diabetes | American Diabetes Association. https://diabetesjournals.org/diabetes/article/37/12/1595/8592/Role-of-Insulin-Resistance-in-Human-Disease (Accessed November 26, 2024b).](https://www.zotero.org/google-docs/?Oy4FWS)  
18. [Underlying Cause of Death, 1999-2020 Results Form. https://wonder.cdc.gov/controller/datarequest/D76;jsessionid=B5D9C6A25586B834E944DF1C9884 (Accessed September 19, 2024c).](https://www.zotero.org/google-docs/?Oy4FWS)






