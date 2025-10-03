#
import pandas
import pandas as pd
import requests
from requests.auth import HTTPBasicAuth

def combineMergeDiag():
    # DXSUM_PDXCONV - Download : - ADNI - Download - Study Data - Assessments - Diagnosis - Diagnostic Summary
    getDiag = pandas.read_csv('DXSUM_PDXCONV.csv')
    # ADNIMERGE - Download : - ADNI - Download - Study Data - Study Info - File: ADNIMERGE - Key ADNI tables merged into one table [ADNI1,GO,2,3]
    getMerge = pandas.read_csv('ADNIMERGE.csv')

    # Create a common column 'PHASE_ORIGPROT'
    getDiag['PHASE_ORIGPROT'] = getDiag['PHASE']
    getMerge['PHASE_ORIGPROT'] = getMerge['ORIGPROT']

    # # Extract the year from the 'EXAMDATE' column
    # getDiag['Year'] = pandas.to_datetime(getDiag['EXAMDATE']).dt.year
    #
    # # Extract the year from the 'EXAMDATE' column
    # getMerge['Year'] = pandas.to_datetime(getMerge['EXAMDATE']).dt.year

    # Merge the data on 'PTID', 'VISCODE', and the new 'PHASE_ORIGPROT' column
    merged_data = pd.merge(getDiag, getMerge, on=["PTID", "VISCODE", "PHASE_ORIGPROT", "EXAMDATE"], how="outer")
    # merged_data.drop(columns="Year")
    merged_data.to_csv("final_data.csv")

# combineMergeDiag()
def extractMetaData():
    getData = pandas.read_csv('CSVData/ADNI_Gene_Expression_Profile.csv', header=None)
    first_8_metadata = getData.iloc[:8]
    #convert rows to cols and cols to rows for metadata
    metaData = first_8_metadata.transpose()

    # save data
    metaData.to_csv('GENE_METADATA.csv', header=False)

# extractMetaData()

def extractData():
    getData = pandas.read_csv('GENE_METADATA.csv')
    geneData = pandas.read_csv('ADNI_Gene_Expression_Profile.csv', header=None)
    geneData = geneData.iloc[9:]
    print(geneData)
    subject_id_year = getData['SubjectID'].astype(str)
    print(subject_id_year)
    geneData.columns = ["ProbeSet", "LocusLink", "Symbol"] + [None] * (geneData.shape[1] - 3)
    for i in range(3, len(geneData.columns)):
        geneData.columns.values[i] = subject_id_year[i - 3]

    print(geneData.head())
    geneData.to_csv('GENE_DATA.csv')

# extractData()
# create common file
def getCommonTable():
    diagnosis_data = pandas.read_csv('CSVData/ADNIMERGE_30Jul2024.csv')
    gene_meta_data = pandas.read_csv('GENE_METADATA.csv')

    gene_meta_data['Gene_Flag'] = 1

    # Extract the year from the 'EXAMDATE' column
    diagnosis_data['Year'] = pandas.to_datetime(diagnosis_data['EXAMDATE']).dt.year

    # Merge the two dataframes on 'Phase', 'Visit', and 'Year'
    merged_data = pandas.merge(diagnosis_data, gene_meta_data,
                           left_on=['ORIGPROT', 'PTID', 'Year'], right_on=['Phase', 'SubjectID', 'YearofCollection'], how="left")

    # Fill NaN values in 'gene_flag' with 0 (indicating no gene expression data)
    merged_data['Gene_Flag'].fillna(0, inplace=True)

    # Drop the 'Year' and 'YearofCollection' columns as they are no longer needed
    merged_data.drop(columns=['Year', 'YearofCollection',"Phase","SubjectID"], inplace=True)

    # Save the result to a new CSV file
    merged_data.to_csv('merged_data.csv', index=False)

#getCommonTable()

# count number of same data
def countGeneticData():
    merged_data = pandas.read_csv('merged_data.csv')

    count = merged_data[(merged_data['Gene_Flag'] == 1)].shape[0]
    print(count)

    # Count by phase where both Gene_Flag and MRI_Flag are 1
    count_by_phase = merged_data[(merged_data['Gene_Flag'] == 1)].groupby('ORIGPROT').size().reset_index(name='Count')
    print("\nCount by phase where Gene_Flag are 1:")
    print(count_by_phase)
#countGeneticData()

# get average isoforms
def averageISOForms():
    gene_data = pandas.read_csv('GENE_DATA.csv')
    # Assuming the numerical columns start from the 4th column onwards
    numerical_columns = gene_data.columns[4:-1]

    # Group by 'Symbol' and calculate the mean for each group
    averaged_gene_data = gene_data.groupby('Symbol')[numerical_columns].mean().reset_index()

    # Save the averaged data to a new CSV file
    averaged_gene_data.to_csv('averaged_gene_data.csv', index=False)

# averageISOForms()

# set MRI Image Flag
def setMRIFlag():
    data = pd.read_csv('merged_data.csv')
    mri_data = pd.read_csv('CSVData/MRI.csv')

    mri_data['MRI_Flag'] = 1

    # Normalize the 'PHASE' column in merged_data to match the format in mri_data
    # mri_data['Phase'] = mri_data['Phase'].str.replace(' ', '')

    # Split the 'Imaging Protocol' column into separate components
    # split_columns = mri_data['Imaging Protocol'].str.split(';', expand=True)

    # Extract the values for each new column
    # mri_data['Acquisition Plane'] = split_columns[0].str.split('=', expand=True)[1]
    # mri_data['Acquisition Type'] = split_columns[1].str.split('=', expand=True)[1]
    # mri_data['Weighting'] = split_columns[2].str.split('=', expand=True)[1]

    # Drop the original 'Imaging Protocol' column if no longer needed
    # mri_data.drop(columns=['Imaging Protocol'], inplace=True)

    # Extract the year from the 'EXAMDATE' column
    data['Year'] = pandas.to_datetime(data['EXAMDATE']).dt.year
    mri_data['Year'] = pandas.to_datetime(mri_data['Acq Date']).dt.year

    # Merge the two dataframes on 'Phase', 'Visit', and 'Year'
    merged_data = pandas.merge(data, mri_data,
                               left_on=['PTID', 'Year'],
                               right_on=['Subject', 'Year'],
                               how='left')

    # Fill NaN values in 'gene_flag' with 0 (indicating no gene expression data)
    merged_data['MRI_Flag'].fillna(0, inplace=True)

    merged_data.to_csv('Merged_MRI_Data.csv', index=False)
# setMRIFlag()

# count number of same data
def countMRIData():
    merged_data = pandas.read_csv('CSVData/Merged_MRI_Data.csv')

    count = merged_data[(merged_data['Gene_Flag'] == 1) & (merged_data['MRI_Flag'] == 1)].shape[0]
    print(count)

    # Count by phase where both Gene_Flag and MRI_Flag are 1
    count_by_phase = merged_data[(merged_data['Gene_Flag'] == 1) & (merged_data['MRI_Flag'] == 1)].groupby(
        'ORIGPROT').size().reset_index(name='Count')
    print("\nCount by phase where both Gene_Flag and MRI_Flag are 1:")
    print(count_by_phase)


countMRIData()

# set pet flag
def setPETFlag():
    data = pandas.read_csv('Merged_Data.csv')
    pet_data = pandas.read_csv('idaSearch_pet.csv')

    # set flag
    pet_data['PET_Flag'] = 1

    # Extract the year from the 'EXAMDATE' column
    data['Year'] = pd.to_datetime(data['EXAMDATE']).dt.year

    # Extract the year from the 'Study Date' column in MRI Image
    pet_data['Year'] = pd.to_datetime(pet_data['Study Date']).dt.year

    # Normalize the 'PHASE' column in merged_data to match the format in mri_data
    pet_data['PHASE'] = pet_data['PHASE'].str.replace(' ', '')

    # Merge the two dataframes on 'Phase', 'Visit', and 'Year'
    merged_data = pd.merge(data,
                           pet_data,
                           left_on=['PHASE', 'PTID', 'Year'],
                           right_on=['PHASE', 'PTID', 'Year'],
                           how='left')

    # Save the merged data to a new CSV file
    merged_data.to_csv('merged_pet_data.csv', index=False)
# setPETFlag()

def countPET():
    merged_data = pandas.read_csv('merged_pet_data.csv')

    count = merged_data[(merged_data['Gene_Flag'] == 1) & (merged_data['PET_Flag'] == 1)].shape[0]
    print(count)

    # Count by phase where both Gene_Flag and MRI_Flag are 1
    count_by_phase = merged_data[(merged_data['Gene_Flag'] == 1)].groupby(
        'PHASE').size().reset_index(name='Count')
    print("\nCount by phase where both Gene_Flag and PET_Flag are 1:")
    print(count_by_phase)
# countPET()


def setFMRIFlag():
    data = pandas.read_csv('Merged_Data.csv')
    fmri_data = pandas.read_csv('idaSearch_FMRI.csv')

    # set flag
    fmri_data['FMRI_Flag'] = 1

    # Extract the year from the 'EXAMDATE' column
    data['Year'] = pd.to_datetime(data['EXAMDATE']).dt.year

    # Extract the year from the 'Study Date' column in MRI Image
    fmri_data['Year'] = pd.to_datetime(fmri_data['Study Date']).dt.year

    # Normalize the 'PHASE' column in merged_data to match the format in mri_data
    fmri_data['PHASE'] = fmri_data['PHASE'].str.replace(' ', '')

    # Merge the two dataframes on 'Phase', 'Visit', and 'Year'
    merged_data = pd.merge(data,
                           fmri_data[['PHASE', 'PTID', 'Year', 'FMRI_Flag']],
                           left_on=['PHASE', 'PTID', 'Year'],
                           right_on=['PHASE', 'PTID', 'Year'],
                           how='left')

    # Set FMRI_Flag to 0 if it's not 1
    merged_data['FMRI_Flag'] = merged_data['FMRI_Flag'].apply(lambda x: 1 if x == 1 else 0)

    # Save the merged data to a new CSV file
    merged_data.to_csv('merged_fmri_data.csv', index=False)
# setFMRIFlag()

def countFMRI():
    merged_data = pandas.read_csv('merged_fmri_data.csv')

    count = merged_data[(merged_data['Gene_Flag'] == 1) & (merged_data['FMRI_Flag'] == 1)].shape[0]
    print(count)

    # Count by phase where both Gene_Flag and MRI_Flag are 1
    count_by_phase = merged_data[(merged_data['Gene_Flag'] == 1) & (merged_data['FMRI_Flag'] == 1)].groupby(
        'PHASE').size().reset_index(name='Count')
    print("\nCount by phase where both Gene_Flag and FMRI_Flag are 1:")
    print(count_by_phase)

# countFMRI()