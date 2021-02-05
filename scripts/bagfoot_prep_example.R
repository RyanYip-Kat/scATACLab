
# bagfoot_prep_example.R
# This script shows how to prepare required data files to run "BagFoot"
# 
# By Songjoon Baek
# Last Update: 2017/06/19

library('bagfoot');

bamfile1= 'example_data/GH1047_PE_sorted.bam';  

cc1 = countReadsBAM(bamfile1);  # counts the number of cuts in the sequence file.

cutcountfile1 = makeCutCountBAM(bamfile1);   # generate a BedGraph file with DNase Cleavages counts

#################### MAPPABILITY FILES ###############################################
# Mappability files were generated using the Mappability code from PeakSeq (Gerstein Lab, http://info.gersteinlab.org/PeakSeq)
# To download, use http://archive.gersteinlab.org/proj/PeakSeq/Mappability_Map/Code/


MAPPABILITY_FILES_DIRECTORY_HG19<<-'~/Data/hg19/35mers/';   # directory of mappability files for the human hg19 ref. genome ( chr1b.out, chr2b.out, etc..)
MAPPABILITY_FILES_DIRECTORY_MM9<<-'~/Data/mm9/35mers/';     # directory of mappability files for the mouse mm9 ref. genome ( chr1b.out, chr2b.out, etc..)
MAPPABILITY_FILES_DIRECTORY_MM10<<-'~/Data/mm10/50mers/';    # directory of mappability files for the mouse mm10 ref. genome ( chr1b.out, chr2b.out, etc..)


#################### HOTSPOT FILES ###################################################
# Hotspots were called by the DNASE2Hotspots software (Ref:  https://www.ncbi.nlm.nih.gov/pubmed/22183609 ), Or 
# you can use MACS to call hotspots and convert narrow peak files to a compatible csv file (use  macs_output2csv.R  script)

hotspotfile1 = 'example_data/fed1_10000_hotspot.csv';   # A peak sites file is comma-separated and must has headers with chromosome, start, end (1-based)
hotspotfile2 = 'example_data/fasted1_10000_hotspot.csv';

combinedhotspotfile = combineTwoHotspots( hotspotfile1, hotspotfile2, 'fed1', 'fasted1');  # Hotspots in hotspotfile1 and hotspotfile2 are combined and saved to a new file name

tabNoMappability = MakeBiasCorrectionTableBAM(        # This function call generates a hexamer bias frequency table from the given bamfile without mappability assumption
	bamfile=bamfile1,
	outfile="Hexamer_fed_mm9_withoutMap.txt",
	refgenome="mm9", 
	np=6, 
	mapdir='', # if mappability file directory is set to empty, mappability is not used. 
	atac=F     # Set TRUE for ATAC data. The default value is FALSE
	);


tabNoMappability = MakeBiasCorrectionTableBAM(        # This function call generates a hexamer bias frequency table from the given bamfile without mappability assumption
	bamfile=bamfile1,
	outfile="decimer_fed_mm9_withoutMap.txt",
	refgenome="mm9", 
	np=10, 
	mapdir='', # if mappability file directory is set to empty, mappability is not used. 
	atac=F     # Set TRUE for ATAC data. The default value is FALSE
	);

tabNoMappability = MakeBiasCorrectionTableBAM(        # This function call generates a hexamer bias frequency table from the given bamfile without mappability assumption
	bamfile=bamfile1,
	outfile="octamer_fed_mm9_withoutMap.txt",
	refgenome="mm9", 
	np=8, 
	mapdir='', # if mappability file directory is set to empty, mappability is not used. 
	atac=F     # Set TRUE for ATAC data. The default value is FALSE
	);

if (dir.exists(MAPPABILITY_FILES_DIRECTORY_MM9)) {
   tabMappability = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_fed_mm9_withMap.txt", 
   refgenome="mm9", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM9);
}

if (dir.exists(MAPPABILITY_FILES_DIRECTORY_MM10)) {
   tabMappability10 = MakeBiasCorrectionTableBAM(bamfile= bamfile1,
   outfile="Hexamer_fed_mm10_withMap.txt", 
   refgenome="mm10", 
   np=6,
   mapdir = MAPPABILITY_FILES_DIRECTORY_MM10);
}

