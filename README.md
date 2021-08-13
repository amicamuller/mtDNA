# MTDNA Repo

CONTENTS OF THIS FILE
---------------------

 * Introduction
 * Requirements
 * Running the analysis

INTRODUCTION
------------
This repo contains two .py files in the **"Scripts"** folder/directory.
These scripts use output files from [mtDNA-Server](https://mtdna-server.uibk.ac.at/index.html) to test the hypothesis that more cases than controls have 'out-of-place' mtDNA variants that are haplogroup-defining but do not define the mtDNA background they are found on. As such, the data to be analysed should contain mtDNA variants from both cases and controls. 

The first script (```making_a_mitomaster_upload_file.py```) creates an input (.txt) file for Mitomaster (_homoplasmies_mitomaster_upload_file.txt_)containing all homoplasmic mtDNA variants and those with a heteroplasmy level >50%, called by the mtDNA-Server pipline. [Mitomaster](https://www.mitomap.org/foswiki/bin/view/MITOMASTER/WebHome) annotates the called variants, and .csv files with these annotated variants can easily be downloaded by the user for further processing. For the user's convenience, the script additionally creates a _heteroplasmies_mitomaster_upload_file_ containing all variants with a heteroplasmy level ≤50% as well as an additional _het_level_cnfrm.xlsx_ file containing all the heteroplasmic variants in your dataset that have previously been confirmed to be pathogenic according to criteria outlined on [Mitopmap] (https://mitomap.org/MITOMAP/ConfirmedCriteria). The _heteroplasmies_mitomaster_upload_file_ can be uploaded to Mitomaster for viewing by the user, but will not be included in the subsequent analyses to test the 'out-of-place' hypothesis.

The second script (```cleaned_making_a_mitomaster_upload_file.py```) tests the 'out-of-place' hypothesis, requiring the mtDNA-Server haplogroups.txt files and the annoated Mitomaster.csv files as input. It generates several output files that can be used for further processing and evaluation by the user. These include:
       1. An annotated mitomaster output file containing all variants included in analysis (i.e. homoplasmic mtDNA variants and those with a heteroplasmy level >90%, called by the mtDNA-Server pipline and annotated using Mitomaster) 
       2. A file listing samples removed from the analysis (e.g. samples with low haplogroup quality scores (≤ 90%) or samples not belonging to the Ancestry selected by the user)
       3. A file with info to make graphs
       4. Input and output files for Fisher's exact tests (testing the differences between the number of cases with and without 'out-of-place' mtDNA varaiants across their mtDNA genomes)
       5. A file with any variants that were not reported on Genbank at the time of analysis
       6. Graphs showing the number of cases and controls with 'out-of-place' mtDNA varaiants across their mtDNA genomes
       

REQUIREMENTS
------------
Before running the scripts in this repository, please ensure you have set up an environment with the requirements specified in the 'mtDNA_sequencing_Python_requirements.txt'.

RUNNING THE ANALYSIS
--------------------

## making_a_mitomaster_upload_file.py
Running this script will make a text file that can be uploaded to Mitomaster for further processing and the assignment of Genbank and Haplogroup frequencies of the individual variants.
### Before running the above-mentioned script, please ensure the following:
       1. You have run your sequencing files (in BAM or fastq format) on the [mtDNA-Server] (https://mtdna-server.uibk.ac.at/index.html).
       2. You have downloaded and saved all individual mtDNA-Server output files in a single folder on your computer.
       3. The names of your mtDNA-Server heteroplasmy files contain the word "heteroplasmies" and the files are in a .txt (i.e. text) format.
       4. The names of your mtDNA-Server homoplasmy (i.e. variant) files contain the word "variants" and are in a .txt format.
       5. The names of your mtDNA-Server haplogroup files contain the word "haplogroups" and the files are in a .txt (i.e. text) format.
       6. Your file names contain no spaces.
       
Follow any on-screen prompts while running the script.

### After sucessfully running the script, please do the following:

 Follow the instructions below to upload your file to Mitomaster, and to subsequently download the required output files for the next step of the analysis:
 
       1. Go to the Mitomaster SNV query page at https://www.mitomap.org/mitomaster/index_snvs.cgi 
       2. Under Step 1: Select 'Compute haplotype (only recommended if supply a full set of SNVs)'
       3. Under Step 2, Option 2: upload your 'homoplasmies_mitomaster_upload_file.txt' file by clicking on the 'Select file' button and choosing your file
       4. Click on "Submit" and wait for the Results page to load 
       5. Scroll down to "Sequence alignment", then click on the word "here" in "click here to show variant details for all"
       6. A new page with "Alignment Details" should load. Click on the "CSV" button to download details of all variants detected (in a .csv format)
       7. Save the CSV file (in the same folder as the files generated by running this script), with a file name that:
               - Has NO SPACES between words (use an underscore instead of spaces e.g. mitomaster_output_2020)
               - Includes the word "mitomaster" (all lower case letters) 
       
       
## cleaned_making_a_mitomaster_upload_file.py

Before proceeding, please ensure the following:

       1. You have run  the 'making_a_mitomaster_upload_file' 
       2. You have downloaded and saved your mitomaster output files (as instructed) to a single folder on your computer
       3. You should have saved a .csv (i.e. comma delimited) file from Mitomaster in the said folder:
         - A .csv file that contains the word "mitomaster" (all lower case letters) 
       4. You have created a .csv file (that contains the words 'sample_info' in the file name) with the details of
       your sequenced samples and the following column headers:
       'sample' (i.e. name of sample) , 'age', 'sex', 'status' (i.e. case or control)
                     
After completeing the instructions above, run the ```cleaned_making_a_mitomaster_upload_file.py``` script and follow any on-screen prompts.
