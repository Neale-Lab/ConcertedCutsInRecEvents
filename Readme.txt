
Note: the scripts in this repository require input from the pipeline from Neale-Lab/OctadRecombinationMapping

Instructions: Download all files from this repository and store them in a directory. 
Unzip any archived files. The scripts can then be run.
Check each script for additional instructions written as comments.

List of relevant output (example files are provided in this repository for each major step):

1. File containing all original annotated events (1.5kb threshold) for every genotype
masterAEventTable_16_8_18
This table is output from the recombination event mapping pipeline - see other repository Neale-Lab/OctadRecombinationMapping.

2. Subset of table 1. File contains annotated events (1.5kb threshold) for every genotype.
Chromatids that were not properly seperated during octad dissection have been filtered out, as they contain many false positives.
List_of_all_clean_events.txt
Made by "Combined_Seg_and_SNPs_heat V2.R"

3. Subset of table 2. File containing annotated events with at least one 6:2 for every genotype. 
List_of_clean_DC_candidate_events_only.txt
Made by "Combined_Seg_and_SNPs_heat V2.R"

4. Table 3 split into segments ready for manual annotation. These files are individual per genotype.
Pre-annotation: e.g. OM_DC_annotation_table.txt
Post-annotation: e.g. OM_DC_annotated_table.txt
Made by "Combined_Seg_and_SNPs_heat V2.R", then annotated manually

5. Master table of annotated events. 
DC_master_table.txt
Made by "Double cut table maker.R"

6. Version of table 5 with segments collapsed back down into individual events.
Collapsed DC master table
Made by "Collapser.R"


Imager scripts 

1. Imager to assist with annotation 
DC imager - pre annotation.R

2. Imager post-annotation (suitable for publication)
DC imager - post annotation.R
