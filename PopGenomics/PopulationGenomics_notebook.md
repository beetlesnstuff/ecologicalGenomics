# Title  

## Author: Alexander Kissonergis
### Affiliation: University of Vermont, Dept. of Plant and Soil Science
### E-mail contact: alexander.kissonergis@uvm.edu


### Start Date: 9.11.2023
### End Date: TBD
### Project Descriptions:   This notebook will document my workflow on the bioinformatics of the Population Genomics section of the Fall 2023 Ecological Genomics course.





# Table of Contents:   
* [Entry 1: 2023-09-11](#id-section1)
* [Entry 2: 2023-09-13](#id-section2)
* [Entry 3: 2023-09-18](#id-section3)


------    
<div id='id-section1'/>   


### Entry 1: 2023-09-11.   

- We reviewed the red spruce study system and the exome capture data
- We discussed the structure of fastq files (DNA sequence, plus the Qscores)
- Using the program FastQC, we analyzed the quality of the sequencing runs for one file.

------    
<div id='id-section2'/>   


### Entry 2: 2023-09-13.  

- After discussing the initial FastQC results, we saw good quality sequence data for most of the read length
- The initial 5 bp or so had more variable base frequencies, and the very end of the reads had slightly lower Q-scores
- Based on this, we set up an analysis to trim the reads using the 'fastp' program
- We ran the bam script 'fastp.sh' for this
- We looked at the html files produced by 'fastp' and compared pre- and post-trimming -- things looked good!
- We ended the day setting up our read mapping of the trimmed and cleaned reads using 'bwa'

------    
<div id='id-section3'/>   


### Entry 3: 2023-09-18.

- Created lab notebook, to be used every day to document the processes taking place in class
- Processed our mapping files using samtools and sambamba
- Updated VIM scripts 'bam_stats.sh' and 'process_bam.sh' (did not finish process_bam.sh, will do Wednesday)
------    
<div id='id-section4'/>  

### Entry 4: 2023-09-20.

-Calculated mapping stats
-Created ANGSD

------    
<div id='id-section5'/>  

### Entry 5: 2023-09-25.

-Calculated SFS and diversity stats
-Edited and adjusted ANGSD, created a new vim script 'ANGSD_doTheta.sh'
-Applied new files to filezilla
-Created R Script to summarize diversity stats
-Added data to 'EcologicalGenomics_ANGSD_Fall2023'

------    
<div id='id-section6'/>  

### Entry 6: 2023-09-27.

-Created new VIM script 'ANGSD_Fst.sh', this script estimates Fst between the red spruce pop 2103 and black spruce:
-Calculated NeW and NeP in 'EcologicalGenomics_ANGSD_Fall2023'
-Created new VIM script 'ANGSD_allRS_poly.beagle.gz'

------    
<div id='id-section7'/>  

### Entry 7: 2023-10-02.

-Walked through results of r script 'Tutorial.Script.9.27.23'
-Adjusted mutation rates in document 'EcologicalGenomics_ANGSD_Fall2023' and added FST scores



