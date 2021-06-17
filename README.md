# MTDNA Repo

CONTENTS OF THIS FILE
---------------------

 * Introduction
 * Installation
 * Requirements

INTRODUCTION
------------
This repo contains two .py files.
These scripts use output files from mtDNA-Server(https://mtdna-server.uibk.ac.at/index.html) to test the hypothesis that more cases than controls have 'out-of-place' mtDNA variants that are haplogroup-definining but not for the mtDNA background they are found in. As such, the data to be analysed should contain mtDNA variants from both cases and controls. 

The first script (**making_a_mitomaster_upload_file.py**) creates an input (.txt) file for MITOMAP/MITOMASTER, containing all homoplamic mtDNA variants called by the mtDNA-Server pipline. MITOMASTER annotates the called variants and .csv files with  these annotated variants can easily be downloaded by the user for further processing.

The second script 

REQUIREMENTS
------------
Before running the scripts in this repository, please ensure you have set up an environment with the requirements specified in the 'requirements.txt file.

## making_a_mitomaster_upload_file.py
Running this script will make a text file that can be uploaded to Mitomaster for further processing and the assignment of Genbank and Haplogroup frequencies of the individual variants.
### Before running the above-mentioned script, please ensure the following:
       1. You have run your sequencing files (in BAM or fastq format) on the mtDNA-Server (https://mtdna-server.uibk.ac.at/index.html).
       2. You have downloaded and saved all individual mtDNA-Server output files in a single folder on your computer.
       3. The names of your mtDNA-Server heteroplasmy files contain the word "heteroplasmies" and the files are in a .txt (i.e. text) format.
       4. The names of your mtDNA-Server homoplasmy (i.e. variant) files contain the word "variants" and are in a .txt format.
       5. Your file names contain no spaces
       
Follow any on-screen prompts while running the script.

### After sucessfully running the script, please do the following:

 Follow the 10 instructions below to upload your file to Mitomaster, and to subsequently download the required output files for the next step of the analysis:
 
       1. Go to the Mitomaster SNV query page at https://www.mitomap.org/mitomaster/index_snvs.cgi 
       2. Under Step 1: Select 'Compute haplotype (only recommended if supply a full set of SNVs)'
       3. Under Step 2, Option 2: upload your 'mitomaster_upload_file.txt' file by clicking on the 'Select file' button and choosing your file
       4. Click on "Submit" and wait for the Results page to load 
       5. Scroll down to "Sequence alignment", then click on the "CSV"  button to download a list of all detected variants (in a .csv format)
           Save the CSV file with a file name that:
            - Has NO SPACES between words (use an underscore instead of spaces e.g. var_file_2020)
            - Includes the word "var_list" (all lower case letters, words separated by underscores) 
            - EXCLUDES the word "mitomaster" (all lower case letters) 
       7. Then again under "Sequence alignment" click on the word "here" in "click here to show variant details for all
       8. A new page with "Alignment Details" should load. Click on the "CSV" button to download details of all variants detected (in a .csv format)
       9. Save the CSV file in the same folder as the CSV file you saved in Step 6, with a file name that:
            - Has NO SPACES between words (use an underscore instead of spaces e.g. mitomaster_output_2020)
            - Includes the word "mitomaster" (all lower case letters) 
       10. Run the script 'processing_mitomaster_output_files'
