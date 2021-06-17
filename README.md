# MTDNA Repo

CONTENTS OF THIS FILE
---------------------

 * Introduction
 * Installation
 * Requirements

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
