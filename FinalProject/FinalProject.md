# Final Project: Determining site of insertion of hCD1a transgene in CD1a mouse model
### Julio Ayala Angulo

## Introduction
Type 2 Diabetes (T2D) is the 8th leading cause of death, with cases expected to reach 1.3 billion globally by 2050[18,21]. Current therapies manage symptoms but fail to cure the disease, with efficacy declining after four years (GRADE Study). The lack of a curative treatment stems from an incomplete understanding of T2D's root causes and complete mechanistic model. Research has largely focused on glycemic and insulinemic dysregulation as core drivers. Adiposity, dyslipidemia, and immune inflammation have seen a resurgence of interest but remain understudied partially due to poor predictive and mechanistic mouse models and the difficulty of elucidating pathophysiological mechanisms in humans. Increased Visceral Adipose Tissue (VAT) mass and dyslipidemia in the form of elevated free fatty acids (FFA), low and very low density lipoproteins (LDL/VLDL) and triglycerides (TG) in T2D patient plasma has been well established[22]). The current literature largely focuses on the impact hyperinsulinemia has on adipose mass and dyslipidemia, therefore assuming that the change in lipids is the outcome and that insulin is the cause[23]. However, there is evidence that the converse is also true, though it has been less studied. For example, plasma FFAs are predictive of hyperglycemia and strongly correlated with insulin resistance (IR)[22,24,25]. Higher FFA and glycerol levels in obese patients without T2D compared to obese with T2D suggests a crucial role in lipid levels in the transition from Insulin Resistant/Prediabetic to T2D [26]. Onset of T2D starts with increasing obesity and FFA as hyperinsulinemia develops to maintain blood sugar and lipid levels. Eventually, a mixture of insulin resistance and pancreatic beta cell atrophy drive hyperglycemia[27]. However, the cause of the atrophy or why adiposity and FFA don't drive hyperglycemia earlier remains unknown. 
CD1a is one of 4 CD1 molecules in humans that allow the immune system to respond to lipid antigens. Mice only have 1 member of the CD1a family: CD1d, which is limited in its antigen presentation to microbial metabolites[51]. Consequently, mice are largely unable to adaptively respond to lipids. In a disease context where dyslipidemia is not only characteristic, but potentially causative, anti-lipid T cell responses can bridge the gap between mouse models and human disease. Of all CD1 molecules, most is known about CD1a. Interestingly, mice not only do not model T2D and oebsity well, but also do not model the low-chronic inflammation or tissue specif inflammation in humans (Elevated IL17a+ and TNFa+) where as mice only express tissue-specfic inflammation of a different type (IFNy+ and TNFa+) [1]. We hypothesized that human CD1a lipid antigen presentation to T cells drives progression of T2D through inflammation. To this extent, we acquired CD1a expressing B6 mice. Under a high fat diet, we saw changes in metabolism and inflammation in CD1a mice. However, it is still unknown if CD1a expression in these mice proproperly mirror humans. The original lab that developed these mice only validated its expression in the skin and thymus as described in humans [Sugita, CD1a atlas]. The insertion included regulatory elements for CD1a in humans, but it's insertion is likely not targeted to chr1 as in humans. Initially, we used flow cytometry to characterized CD1a+ cells proportions in mice (Figure 1). This data validates the high expression in tissues previsly seen in humans (Thymus and skin) while also basal expression in liver and pancreas (not described) while notably absent in lymphatic tissues. Additonally, there is a clear lower expression of CD1a in female mice. Depending on the insertion site(s), this could be an artifact of the insertions location as well as its expression in pathologically relevant tissues. 

![CD1a Expression Across tissues](https://github.com/jayalaanUCI/ee282/blob/FinalProject/FinalProject/output/figures/CD1a%20Expression.png?raw=true)

### Goals
1. Chromosomes Harboring CD1a transgene
    -Primarily, I wish to determine which chromosomes have copies of CD1a. This lets me know which trangenic mice, like RAG1KO mice,  I am able to breeed with the CD1a mice. It mght also shed light on differeces across sexes if some copies lie on the X or Y chromosomes.
    -I aim to output a table per sequence with chromosome and location afte raligbnemnt that I can graph with ggplot2.
2. Vizualize spread of of insertion across the mouse genome using [IGV](https://igv.org/). 
## Methods and Results
### Sample Preparation and Targeted Locus Amplification
DNA samples for sequencing where generated from CD1a mouse splenocytes, lymph nodes, and blood. Tissues were macerated and red blood cells were lysed with RBC lysis buffer. Mononuclear cells were then put through this [TLA protocol](https://link.springer.com/protocol/10.1007/978-1-4939-6442-0_13) for isolation of DNA. In short, cells were lysed in the presensence of paraformaldehyde to fix DNA to proximal transcription factors, enahncers, or aother regulatory proteins. Samples were then digested and subsequebtly ligated to link DNA sequences attached to the same regualtory elements. Proteins were digeted before ligating DNA once more to form circular DNA. CD1a-enriched amplicons were generated through PCR usign inverted genotying primers. DNA amplicons were purified via AMPure XP beads before beign sent for Nanopore Sequencing at [Plasmidsaur](https://plasmidsaurus.com/).
### Required Packages and alnguages
I will be using Mamba as a package manager through which most of these packages were downloaded.
| Packacage/Utility |      Purpose      |
| ----------------- | ----------------- |
| Perl              |                   |
| NanoQC            |   Quality checks  |
| Minimap2          |   Alignement      |
| Bioawk            | Data handling     |
| R and ggplot2     |  Graphing         |
| Fasize/Fafilter   | Seq reports and filtering |
| IGV (online)      | Alignment         |

Data and output will be posted to [github](https://github.com/jayalaanUCI/ee282). 

## Code 
```
# !/bin/bash
#Uploaded raw FastQ 7PCS4R_1_CD1a_TLA_amplicon_primer_3.fastq to ./FinalProject/data/raw/ before chaning name before and compressing
cd ~/myrepos/FinalProject/data/raw
mv 7PCS4R_1_CD1a_TLA_amplicon_primer_3.fastq CD1aTLA.fastq | gzip 
```
Next, my goal is to run NanoQC to out out quality check graphs specific to the nanopore platform and use faSize to get summary info of the reads.

```
nanoqc -o output/figures/ data/raw/CD1aTLA.fasta.gzip


```






