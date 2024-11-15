As a tool to understand the mechanism through which the immune system can contribute to Type 2 Diabetes proposal, we decided to use a previously generated human CD1a expressing mouse (hCD1a). CD1a is an antigen presenting receptor that can present lipid antigens to T cells, which in diseases or conditions resulting in dyslipidemia, could produce autoimmune responses.  The hCD1a mouse was originally made by  the Sugita lab to study skin immunity [(Kobayashi et al. 2011\)]). The Sugita lab verified that the receptor is expressed in mouse langerhans cells, the thymus, and mature dendritic cells as seen in humans. They also showed it is inducible by GM-CSF as in humans. They generated the mice by inserting a 5.8 kbp CD1a gene consisting of 6 exons flanked by \~ 1.4 kbp without specifying location or method of insertion. Our own data suggest sex specific differences in expression as well as expression in liver, adipose, and intestines of both male and female mice. With limited information on what regulatory elements and location of the insert, we sought to determine the extent the site or sites of insertion as well as what other elements that might affect expression were included.  

## Method of Analysis and the Data
To that end, I proposed using Targeted Locus amplification PCR. The protocol generates amplicons of nearby DNA regions based on their respective Topological Associating Domain through fixation to regulating proteins [(Hottentot et al. 2017\)]. Using reverse direction primers, any complications with CD1a primer complementary sequences are amplified. I then sent these amplicons ranging from100 bp ro 1.5kb to adapter-free nanopore sequencing through a company called Plasmidsaur. While they produced a “consensus sequence” on the assumption I was sequencing a plasmid, it is unlikely this sequence is accurate.
I currently have the raw reads in a FASTA format. Using linux based workflows:

>Use [NanoQC](https://github.com/wdecoster/nanoQC) for generation of quality control plots and output summary data of read quality.
>
>Using [chopper](https://github.com/wdecoster/chopper) I can trim low quality reads, and short reads that unlike to selectively align
>
>Using [BWA](https://github.com/lh3/bwa) to produce aligned data to both a human reference genome [GRCh38](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14) and  mouse reference genome [GRCm39](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz)
## Visualization

>IGV would suffice for general chromosome location and extrapolation. ![TLA](https://media.springernature.com/full/springer-static/image/chp%3A10.1007%2F978-1-4939-6442-0_13/MediaObjects/332385_1_En_13_Fig3_HTML.gif?as=webp)
>
>But to look at specific sequences, introns, and generate figures, [ggMSA](http://yulab-smu.top/ggmsa/) might be more useful. It would also allow me to practice utilizing R.
  
## Feasability and End Goal

My primary goal is to generate an alignment map of the insertion site to inform future crosses with other mouse strains (make sure the CD1a and RAG1 are on separate chromosomes) in the form of a histogram or heat map highlighting the concentration of overlapping reads. Secondary, determining the specific human regulatory elements: promoter sites, enhancer/super-enhancer sites, transcription binding sites, and/or repressor binding sites that could better inform our understanding of the models limitations and perhaps differences between sexes in the form of a list. The fact that I already have a FASTA file of raw reads means that the main limiting factor is the implementation and expertise. Thus, once I develop the confidence and skills carrying out the analysis should be relatively fast.

### References:
[Hottentot QP, van Min M, Splinter E, White SJ. 2017\. Targeted Locus Amplification and Next-Generation Sequencing. In *Genotyping: Methods and Protocols* (eds. S.J. White and S. Cantsilieris), pp. 185–196, Springer, New York, NY https://doi.org/10.1007/978-1-4939-6442-0\_13 (Accessed October 16, 2024).]  
[Kobayashi C, Shiina T, Tokioka A, Hattori Y, Komori T, Kobayashi-Miura M, Takizawa T, Takahara K, Inaba K, Inoko H, et al. 2011\. GM-CSF-independent CD1a expression in epidermal Langerhans cells: evidence from human CD1A genome-transgenic mice. *J Invest Dermatol* **132**: 10.1038/jid.2011.280.]


