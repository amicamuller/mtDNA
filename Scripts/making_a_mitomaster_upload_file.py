print("Running 'making_a_mitomaster_upload_file' script...")

import easygui, os, pandas as pd, numpy as np, glob

# EASYGUI POP-UPS WITH REQUIREMENTS

welcome_msg = ("Running this script will make a text file that can be uploaded to Mitomaster for further processing "
               "and the assignment of Genbank and Haplogroup frequencies of the individual variants.\n\n"
               )
welcome_title = 'Welcome!\n\n'

easygui.msgbox(welcome_msg, welcome_title)



# use easygui to allow user to select the folder where they have saved the data for this script
source_dir_of_required_files = easygui.diropenbox(msg="Welcome! For this script to run, please select the "
                                                      "\"Required_files\" folder that was downloaded together with "
                                                      "this script", title="Select folder")


msg = ("Before proceeding, please ensure the following:\n\n"
       "1. You have run your sequencing files (in BAM or fastq format) on the mtDNA-Server ("
       "https://mtdna-server.uibk.ac.at/index.html)\n\n "
       "2. You have downloaded and saved all individual mtDNA-Server output files in a single folder on your "
       "computer\n\n"
       "3. The names of your mtDNA-Server heteroplasmy files contain the word \"heteroplasmies\" and the files are in "
       "a .txt "
       "(i.e. text) format)\n\n"
       "4. The names of your mtDNA-Server homoplasmy (i.e. variant) files contain the word \"variants\" and are in a "
       ".txt format\n\n "
       "              Click \"Continue\" if you have ensured all of the above\n\n")

title = 'Please Confirm'
if easygui.ccbox(msg, title):  # show a Continue/Cancel dialog
    pass  # user chose Continue
else:  # user chose Cancel
    exit()


# Preparing mtDNA heteroplasmy files from mtDNA Server

# use easygui to allow user to select the folder with their mtDNA server files
source_dir = easygui.diropenbox(msg='Please select the folder in which you have saved all your mtDNA-Server'
                                    ' files.', title='Select a folder')

os.chdir(source_dir)

# check if directory contains heteroplasmies.txt files

# if there are no heteroplasmy files print error message and move on
if len(glob.glob('*heteroplasmies*.txt')) == 0:
    # bring up error message asking the user to double check whether they have heteroplasmy files for their samples.
    if easygui.ccbox(msg="OOPS. Looks like there aren't any heteroplasmy files in your selected folder\n\n"
                   "Please check:\n\n"
                   "1. That your samples have heteroplasmy files\n\n"
                   "2. That you have selected the correct folder\n\n"
                   "3. That your heteroplasmy files (saved from mtDNA-Server) are .txt files and contain the word "
                         "'heteroplasmies'\n\n "
                   "If you are sure that your samples don't have heteroplasmy files, click 'Continue', otherwise "
                         "click 'Cancel' and ensure the above\n\n",
                   title="ERROR: We couldn't find your files!", choices=["Continue", "Cancel"]):
        pass
    else:
        exit() # exit if the user has heteroplasmy files but didn't choose the right folder

# else make one dataframe from individual heteroplasmies.txt files
else:
    # import all heteroplasmyy files and concatenate them to make one pandas dataframe
    all_heteroplasmies_files = glob.glob(source_dir + '\*heteroplasmies*.txt')
    heteroplasmies = pd.concat((pd.read_csv(f, sep='\t') for f in all_heteroplasmies_files), sort=False)
    print('Processing heteroplasmies...')

    # Filter for reliable type 1 heteroplasmies (heteroplasmies not in Low complexity regions)

    filter_type1_heteroplasmies = heteroplasmies.loc[:, 'TYPE'] == 1
    type1_heteroplasmies = heteroplasmies[filter_type1_heteroplasmies]




    # list of columns we want from the filtered heteroplasmy file
    wanted_heteroplasmy_columns = ['ID', 'POS', 'rCRS', 'TOP-BASE-FWD', 'MINOR-BASE-FWD', 'HET-LEVEL','NUMTs']
    filtered_heteroplasmies = type1_heteroplasmies[wanted_heteroplasmy_columns]
    filtered_heteroplasmies = filtered_heteroplasmies.astype({'POS': str})
    # New Column (Minor is not Ref):	 =if(ref=minor fwd; '' ; minor fwd)
    filtered_heteroplasmies.loc[:, 'minor_not_ref'] = pd.np.where \
        ((filtered_heteroplasmies.loc[:, 'rCRS'] == filtered_heteroplasmies.loc[:, 'MINOR-BASE-FWD']),
         '', filtered_heteroplasmies.loc[:, 'MINOR-BASE-FWD'])

    # Step 2: identify heteroplasmic TOP SNPs

    # make a df for top snps
    top_snp = filtered_heteroplasmies.copy()
    # filter for values = '' in 'minor_not_ref' column
    filter_blanks = top_snp.loc[:, 'minor_not_ref'] == ''
    # apply filter
    top_snp_df = top_snp[filter_blanks]
    # make a copy to avoid iloc warning
    top_snp_df = top_snp_df.copy()
    # add SNP column
    top_snp_df.loc[:, 'SNP'] = top_snp_df.loc[:, 'POS'] + top_snp_df.loc[:, 'TOP-BASE-FWD']

    # work out the heteroplasmy level of top SNPs
    # The heteroplasmy level in the file gives the heteroplasmy % of the minor allele (i.e. minor forward base)
    top_snp_df.loc[:, 'HET-LEVEL'] = (1 - top_snp_df.loc[:, 'HET-LEVEL']) * 100

    # rename TOP-Base-fwd column to var
    final_top_snps = top_snp_df.rename(columns={'TOP-BASE-FWD': 'var'})
    # reset the index to start at 0
    final_top_snps.index = np.arange(0, len(final_top_snps))

    # Step 3: identify heteroplasmic MINOR SNPs

    # make a df for minor snps
    minor_snp = filtered_heteroplasmies.copy()
    # filter for values not = '' in 'minor_not_ref' column
    filter_blanks = minor_snp['minor_not_ref'] != ''
    # apply filter
    minor_snp_df = minor_snp[filter_blanks]
    # make a copy to avoid iloc warning
    minor_snp_df = minor_snp_df.copy()
    # add SNP column
    minor_snp_df.loc[:, 'SNP'] = minor_snp_df.loc[:, 'POS'] + minor_snp_df.loc[:, 'MINOR-BASE-FWD']

    # work out the heteroplasmy level of minor SNPs
    # The heteroplasmy level in the file gives the heteroplasmy % of the minor allele (i.e. minor forward base)
    minor_snp_df.loc[:, 'HET-LEVEL'] = minor_snp_df.loc[:, 'HET-LEVEL'] * 100

    # rename minor-Base-fwd column to var
    final_minor_snps = minor_snp_df.rename(columns={'MINOR-BASE-FWD': 'var'})
    # reset the index to start at 0
    final_minor_snps.index = np.arange(0, len(final_minor_snps))



    # Step 4: Join the minor and top SNP data frames

    # ask the user which percentage of heteroplasmy  they want to include in their analysis
    user_het_level = easygui.enterbox(
        msg="Please enter the heteroplasmy level you'd like to include in your analysis (e.g. 50).\n\n"
            "If you are looking at inherited mtDNA variation, we recommend using a 50% heteroplasmy cut "
            "off or higher", title='Percentage heteroplasmy cut-off', default='50',
        strip=True)
    user_het_level_int = int(user_het_level)

    processed_heteroplasmies = pd.concat([final_top_snps, final_minor_snps], ignore_index=True, sort=False)


    # df with heteroplasmies >= user heteroplasmy cut-off percent
    processed_heteroplasmies_more_than_50_percent = processed_heteroplasmies[
        processed_heteroplasmies.loc[:, 'HET-LEVEL'] >= user_het_level_int]
    # rename rCRS column to ref
    processed_heteroplasmies_more_than_50_percent = processed_heteroplasmies_more_than_50_percent \
        .rename(columns={'rCRS': 'ref'})
    final_heteroplasmies = processed_heteroplasmies_more_than_50_percent[['ID', 'POS', 'ref', 'var', 'SNP']]

    # Make a dataframe for a mitomaster upload file with heteroplasmic SNPs < user selected cut-off
    processed_heteroplasmies_less_than_50_percent = processed_heteroplasmies[(processed_heteroplasmies.loc[:, 'HET-LEVEL'] > 10) & (processed_heteroplasmies.loc[:, 'HET-LEVEL'] < user_het_level_int)]
    # rename rCRS column to ref-
    processed_heteroplasmies_less_than_50_percent = processed_heteroplasmies_less_than_50_percent \
        .rename(columns={'rCRS': 'ref'})
    final_minor_heteroplasmies = processed_heteroplasmies_less_than_50_percent[['ID', 'POS', 'ref', 'var', 'SNP']]
    # save df in the format required for mitomap
    final_minor_heteroplasmies = final_minor_heteroplasmies.rename(columns={'ID': 'sample', 'POS': 'pos'})
    final_minor_heteroplasmies = final_minor_heteroplasmies[['sample', 'pos', 'ref', 'var']]
    # make sure pos are integers not floats
    final_minor_heteroplasmies.loc[:, 'pos'] = pd.to_numeric(final_minor_heteroplasmies.loc[:, 'pos']).astype(int)
    # remove _rCRS from sample name if input files were fastq files
    final_minor_heteroplasmies['sample'] = final_minor_heteroplasmies['sample'].str.replace('_rCRS', '')
    # change reference at position 3107, to the desired mitomap input of a deletion (ie. :) instead of an N for
    # mitomaster to read reference
    final_minor_heteroplasmies = final_minor_heteroplasmies.copy()
    final_minor_heteroplasmies['ref'].mask(final_minor_heteroplasmies['ref'] == 'N', ':', inplace=True)

    # make a dataframe of heteroplamic variants that can be added to the mitomaster output file
    heteroplasmic_merge = processed_heteroplasmies[['ID', 'SNP', 'HET-LEVEL']]
    # remove _rCRS from sample name if input files were fastq files
    heteroplasmic_merge.loc[:,'ID'] = heteroplasmic_merge.loc[:,'ID'].str.replace('_rCRS', '')
    print(heteroplasmic_merge[['ID']])
    # add new column that can be used for merging dataframes later on
    heteroplasmic_merge.loc[:, 'haplo_merge_ref'] = heteroplasmic_merge.loc[:,"ID"] + '_' + heteroplasmic_merge.loc[:,"SNP"]
    # rename columns
    heteroplasmic_merge = heteroplasmic_merge.rename(columns = {'HET-LEVEL':'heteroplasmy_%' })
    print(heteroplasmic_merge)

# Step 5: Preparing mtDNAServer variant files

# check if there are homoplasmic (i.e. variant) txt files in the selected folder
if len(glob.glob('*variants*.txt')) == 0:
    easygui.msgbox("OOPS. Looks like there aren't any variant files in your selected folder\n\n"
                   "We can't create a Mitomaster upload file for you without your variant files...\n\n"
                   "So please check:\n\n"
                   "1. That you have selected the correct folder\n\n"
                   "2. That your variant files (saved from mtDNA-Server) are .txt files and contain the word 'variants'\n\n",
                   "ERROR: We couldn't find your files!", "OK")

print('Processing homoplasmic variants...')

# make one dataframe from individual variant.txt (i.e. homoplasmy) files
all_variant_files = glob.glob(source_dir + '\*variants*.txt')
variants = pd.concat((pd.read_csv(f, sep='\t') for f in all_variant_files), sort = False)
# fill na values with 0
variants = variants.fillna(0)
# drop empty rows
variants= variants[variants['POS'] != 0]
# make sure pos are integers not floats
variants = variants.astype({'POS': int})
# split Mut column into ref and var
variants[['ref', 'var']] = variants.MUT.str.split(">", expand=True)
# make POS column into str and make a SNP column
variants = variants.copy()
variants.loc[:, 'POS'] = variants.loc[:, 'POS'].apply(str)
variants.loc[:, 'SNP'] = variants.loc[:, 'POS'] + variants.loc[:, 'var']
final_variants = variants[['ID', 'POS', 'ref', 'var', 'SNP', 'Locus']]

# make a dataframe of homoplasmic variants that can be added to the mitomaster output file
homoplasmic_merge = final_variants[['ID', 'SNP']]
# remove _rCRS from sample name if input files were fastq files
homoplasmic_merge.loc[:,'ID'] = homoplasmic_merge.loc[:,'ID'].str.replace('_rCRS', '')
# add new column that can be used for merging dataframes later on
homoplasmic_merge.loc[:,'haplo_merge_ref'] = homoplasmic_merge.loc[:,"ID"] + '_' + homoplasmic_merge.loc[:,"SNP"]
# add a heteroplasmy-level column where values = 100
homoplasmic_merge.loc[:,'heteroplasmy_%'] = 100

print(homoplasmic_merge[['heteroplasmy_%']])

# Step 6: Make a Mitomaster upload file with identified variants

print('Making a Mitomaster upload file...')

if len(glob.glob('*heteroplasmies*.txt')) != 0:
    # combine variants df with final heteroplasmies file
    homoplasmies = final_variants.append(final_heteroplasmies, sort=False)
    # save df in the format required for mitomap
    homoplasmies = homoplasmies.rename(columns={'ID': 'sample', 'POS': 'pos'})
    mitomaster_upload_file = homoplasmies[['sample', 'pos', 'ref', 'var']]
    # make sure pos are integers not floats
    mitomaster_upload_file.loc[:, 'pos'] = pd.to_numeric(mitomaster_upload_file.loc[:, 'pos']).astype(int)
    # remove _rCRS from sample name if input files were fastq files
    mitomaster_upload_file.loc[:,'sample'] = mitomaster_upload_file.loc[:,'sample'].str.replace('_rCRS', '')
    # change reference at position 3107, to the desired mitomap input of a deletion (ie. :) instead of an N for
    # mitomaster to read reference
    mitomaster_upload_file = mitomaster_upload_file.copy()
    mitomaster_upload_file['ref'].mask(mitomaster_upload_file['ref'] == 'N', ':', inplace=True)
else:
    homoplasmies = final_variants
    # save df in the format required for mitomap
    homoplasmies = homoplasmies.rename(columns={'ID': 'sample', 'POS': 'pos'})
    mitomaster_upload_file = homoplasmies[['sample', 'pos', 'ref', 'var']]
    # make sure pos are integers not floats
    mitomaster_upload_file = mitomaster_upload_file.astype({'pos': int})
    # remove _rCRS from sample name if input files were fastq files
    mitomaster_upload_file['sample'] = mitomaster_upload_file['sample'].str.replace('_rCRS', '')
    # change reference at position 3107, to the desired mitomap input of a deletion (ie. :) instead of an N for
    # mitomaster to read reference
    mitomaster_upload_file = mitomaster_upload_file.copy()
    mitomaster_upload_file['ref'].mask(mitomaster_upload_file['ref'] == 'N', ':', inplace=True)


destination_dir = easygui.diropenbox(
    msg='Please select the folder in which you want to save your Mitomaster upload file')  # easygui is used to get a path

# change directory to the directory chosen by the user
os.chdir(destination_dir)

# save mitomaster upload file/s as .txt file without index column into the directory selected by the user
mitomaster_upload_file.to_csv("homoplasmies_mitomaster_upload_file.txt", sep='\t', index=False)  # exporting the dataframe to txt
final_minor_heteroplasmies.to_csv("heteroplasmies_mitomaster_upload_file.txt", sep='\t', index=False)

# make a file with all variant's heteroplasmy levels that can be merged with the mitomaster output file
all_het_levels = pd.concat([heteroplasmic_merge, homoplasmic_merge], axis=0, sort = False)

het_levels = all_het_levels[['haplo_merge_ref', 'heteroplasmy_%']]
# save df as .csv file in mtDNA-server folder
het_levels.to_csv(destination_dir + '\het_levels.csv', index=False)

# make a file with all variant's heteroplasmy levels that flags previously confirmed pathogenic variants
# use path chosen by user to the "Required_files" folder to open the file
cnfrm_var_file_location = source_dir_of_required_files + "\mitomap_cnfrm_mutations.txt"
# import mitomap_cnfrm_mutations.csv file
cnfrm_var = pd.read_csv(cnfrm_var_file_location, sep='\t', engine='python')

# merge dfs
het_levels_cnfrm = pd.merge(all_het_levels, cnfrm_var, on ='SNP', how ='inner')
het_levels_cnfrm = het_levels_cnfrm[['ID', 'SNP','heteroplasmy_%', 'confirmed_disease','last_ status_update']]
het_levels_cnfrm.to_excel(destination_dir + '\het_levels_cnfrm.xlsx', index=False)

print('Mitomaster upload file/s saved successfully to '+ destination_dir)


easygui.msgbox("SUCCESS!! \n\n"
               "Your files (homoplasmies_mitomaster_upload_file.txt, heteroplasmies_mitomaster_upload_file.txt and "
               "het_levels.csv) are complete and have been saved to the folder you have "
               "selected!\n\n"
               "Please follow the 10 instructions below to upload your file to Mitomaster and to subsequently download "
               "the required output files:\n\n "
               "1. Go to the Mitomaster SNV query page at https://www.mitomap.org/mitomaster/index_snvs.cgi \n\n"
               "2. Under Step 1: Select 'Compute haplotype (only recommended if supply a full set of SNVs)' \n\n"
               "3. Under Step 2, Option 2: upload your 'homoplasmies_mitomaster_upload_file.txt' file by clicking on the "
               "'Select file' button and choosing your file\n\n "
               "4. Click on \"Submit\" and wait for the Results page to load \n\n"
               "5. Scroll down to \"Sequence alignment\", then click on the \"CSV\"  button to download a list of all "
               "detected variants (in a .csv format)\n\n "
               "6. Save the CSV file (in the same folder as the files created using this script) with a file name that:\n\n"
               "  - Has NO SPACES between words (use an underscore instead of spaces e.g. var_file_2020)\n\n"
               "  - Includes the word \"var_list\" (all lower case letters, words separated by underscores) \n\n"
               "  - EXCLUDES the word \"mitomaster\" (all lower case letters) \n\n"
               "7. Then again under \"Sequence alignment\" click on the word \"here\" in \"click here to "
               "show variant details for all\"\n\n "
               "8. A new page with \"Alignment Details\" should load. Click on the \"CSV\" button to download details "
               "of all variants detected (in a .csv format)\n\n "
               "9. Save the CSV file in the same folder as the CSV file you saved in Step 6, with a file name that:\n\n"
               "  - Has NO SPACES between words (use an underscore instead of spaces e.g. mitomaster_output_2020)\n\n"
               "  - Includes the word \"mitomaster\" (all lower case letters) \n\n"
               "10. Run the script 'processing_mitomaster_output_files'\n\n",
               title='File saved successfully', ok_button='OK. Got it!')


