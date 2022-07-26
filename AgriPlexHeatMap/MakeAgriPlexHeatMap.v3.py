####   getHeatmap4Agriplex.py   Creted by Miyuraj H. Hikkaduwa Withanage    #######
####   contact: mwithanage@inari.com Jun/ 16/2022
####   This script will take two input files
import sys
#        1. A SNP count file (first 12 lines deleted)
#input_count_file = '220519_VH00354_25_AAAL5FTHV.2.csv'
input_count_file = sys.argv[1]
#input_count_file = '220719_VH00354_32_AAANH3MHV.csv'
#        2. A SNP pass/fail file (First 12 lines deleted)
#input_call_file = '220519_VH00354_25_AAAL5FTHV.2.csv.calls.csv'
input_call_file = sys.argv[2]
#input_call_file ='220719_VH00354_32_AAANH3MHV.csv.calls.csv'
# logic:
#       First it will open the SNP count file as a Pandas dataframe and select the first three columns
#       Next, the same dataframe is subsetted so that it only contaisn the counts.
#       We get the row sum for the file.
#       The first three coulmns and row sum are merged and saved in memory as "first3col"
#       Next, we open the SNP pass/fail file.
#       It is subsetted so that first two columns are discarded. The index column remains.
#       It is easier to count "Fail" strings when the dataframe is transposed.
#       we transpose the dataframe and use value count function for counting the occurrences of "Fail".
#       to get the pass count we subtract that from the total
#
#
####
####
####
####
import sys

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# first argument: is SNP counts
# second argument: is SNP pass/fail (call file)
# open the snp call count file. sample names are index as index_col = 0
# open the snp call count file. sample names are index as index_col = 0
#dataset_counts = pd.read_csv('./220519_VH00354_25_AAAL5FTHV.2.csv',index_col=0)
dataset_counts = pd.read_csv(input_count_file,index_col=0)
# save first three columns[0th, 1st, 2nd]. These are meta data.
first3col = dataset_counts.iloc[:,0:2]

# create a new datafframe without the first three columns and getting the row sum
dataset_counts_without_first_3 = dataset_counts.iloc[:,2:dataset_counts.shape[1]]
first3col['sum'] = dataset_counts_without_first_3.iloc[:,1:dataset_counts_without_first_3.shape[1]].sum(axis=1)


####### pass/fail ####
# openning the pass/fail file
#dataset_passfail = pd.read_csv('./220519_VH00354_25_AAAL5FTHV.2.csv.calls.csv',index_col=0)
dataset_passfail = pd.read_csv(input_call_file,index_col=0)

# save first three columns
#first3col_fails = dataset_passfail.iloc[:,0:2]

# creating a df without first three metadata columns
dataset_fails_without_first_3 = dataset_passfail.iloc[:,2:dataset_passfail.shape[1]]

# iterate over rows of dataset_fails_without_first_3 to count "Fails"
# and subtract that from the total number of cell in the row
# define the index list for pass/fail dataframe
index_list_pass_fail = []
# define the pass count list for pass/fail dataframe
passCount_list = []
newpasscount_summary_df = pd.DataFrame()

# a list to hold percentage of passed snps
SNPspassedAfterthreshold = []
SNPspassedAfterthreshold_logical = []
for index, row in dataset_fails_without_first_3.iterrows():
    # transposing makes it easier to count the occurance of "Fail"
    therowt = row.transpose()
    #print(therowt)
    # this make a frequency table
    therowValcount = therowt.value_counts()
    #print(therowValcount)
    # getting the pass count
    passed = dataset_fails_without_first_3.shape[1] - therowValcount
    #print(passed)
    #print(index)
    # passed is an pandas series object. It looks like a freq table


    # adding the index label into the list
    index_list_pass_fail.append(index)
    # Now I want to return passed pandas series keys as a list and check if "FAIL" in that list
    if "FAIL" in passed.keys().to_list():
        # adding the pass count into the list
        passCount_list.append(passed['FAIL'])
        #print(passed['FAIL'])
        passpercentage = 100*(passed['FAIL'])/(dataset_fails_without_first_3.shape[1])
        SNPspassedAfterthreshold.append(passpercentage)
        SNPspassedAfterthreshold_logical.append(passpercentage > 75)
    else:
        passCount_list.append(dataset_fails_without_first_3.shape[1])
        SNPspassedAfterthreshold.append(100)
        SNPspassedAfterthreshold_logical.append(True)

print(first3col)
first3col.index = first3col.index.map(str)
newpasscount_summary_df['Sample Name'] = index_list_pass_fail
newpasscount_summary_df['pass'] = passCount_list
newpasscount_summary_df['pass_percentage'] = SNPspassedAfterthreshold
newpasscount_summary_df['good.sample.or.bad.sample'] = SNPspassedAfterthreshold_logical
#print(newpasscount_summary_df.shape)

newpasscount_summary_df['Sample Name'] = newpasscount_summary_df['Sample Name'].astype(str)
#first3col['Sample Name'] = first3col['Sample Name'].astype(str)
# joining two summary tables on "Sample Name" key
joinedDf = first3col.join(newpasscount_summary_df.set_index('Sample Name'), on='Sample Name')
print(joinedDf)


#########   Creating the Heatmap Layout ###############

rowlabel_list = []
collabel_list = []
for each in joinedDf['Well Id']:
    rowlabel_list.append(each[0])
    collabel_list.append(each[1:3])

joinedDf['row'] = rowlabel_list
joinedDf['col'] = collabel_list




#print(joinedDf.columns)
joinedDf.to_csv("updated.summary."+input_count_file,index=False)

#######

####### creating supplimentary summary files


########
##########
# group by plate id
gb = joinedDf.groupby('Plate Id')
#print(gb.get_group('UDI15PL013'))

#exit()

#print(unique_plateID)
#exit()

def heatmapGenerator(plateID,groups):


    platedf = groups.get_group(plateID)


    # print([gb.get_group(x) for x in gb.groups][0])

    ##### defining the dataframe for heatmap

    heatmapDF_readDepth = pd.DataFrame(columns=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                                                '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                                                '21', '22', '23', '24'],
                                       index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                                              'P'])

    heatmapDF_pass = pd.DataFrame(columns=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                                           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                                           '21', '22', '23', '24'],
                                  index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                                         'P'])

    heatmapDF_pass_threshold = pd.DataFrame(columns=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10',
                                           '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                                           '21', '22', '23', '24'],
                                  index=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',
                                         'P'])



    # filling up the dataframe with read depth (total)
    for ind, val in platedf.iterrows():
        heatmapDF_readDepth.loc[val[6], val[7]] = val[2]
        heatmapDF_pass.loc[val[6], val[7]] = val[3]

    # annotating missing values with -999
    heatmapDF_readDepth = heatmapDF_readDepth.fillna(-999)
    #print(heatmapDF_readDepth)
    # masking
    # for that the dataframe needs to be converted to numpy array

    heatmapDF_readDepth_np = heatmapDF_readDepth.to_numpy()
    masked_readdepth_df = np.ma.masked_where(heatmapDF_readDepth_np == -999, heatmapDF_readDepth_np)
    # now this masked_readdepth_df contains several variables
    # we need masked_readdepth_df.mask
    f,ax = plt.subplots(figsize=(10,5))
    yticks = heatmapDF_readDepth.index
    xticks = heatmapDF_readDepth.columns
    #print(yticks)

    #ax.tick_params(axis='y', rotation=90)
    #sns.set(font_scale=2)
    #ax.set_yticklabels(yticks, rotation=0)
    sns.heatmap(heatmapDF_readDepth, yticklabels=yticks,xticklabels=xticks,annot=True, annot_kws={"size": 4},fmt='g',
                mask=masked_readdepth_df.mask).set(title='Read Depth-Based Heatmap (Plate: ' +plateID+')')

    plt.yticks(rotation=0)
    plt.xticks(rotation=90, size=7)
    ax.xaxis.tick_top()
    f.savefig(plateID+".readDepth.heatmap.pdf", format="pdf", bbox_inches="tight")

    #plt.xaxis.set_ticks_position("top")
    #plt.show()
    #exit()

    ###########################  pass figure ######
    # filling up the dataframe with pass
    #for ind, val in platedf.iterrows():
    #    heatmapDF_pass.loc[val[6], val[7]] = val[3]

    # annotating missing values with -999
    #print(heatmapDF_pass)
    heatmapDF_pass = heatmapDF_pass.fillna(-999)
    # print(heatmapDF_pass)
    # masking
    # for that the dataframe needs to be converted to numpy array

    heatmapDF_pass_np = heatmapDF_pass.to_numpy()
    #print(heatmapDF_pass_np)
    masked_pass_df = np.ma.masked_where(heatmapDF_pass_np == -999, heatmapDF_pass_np)

    # now this masked_readdepth_df contains several variables
    # we need masked_readdepth_df.mask
    #f1 = plt.figure(figsize=(10, 5))
    f1, ax = plt.subplots(figsize=(10, 5))
    yticks1 = heatmapDF_pass.index
    xticks1 = heatmapDF_pass.columns
    #print(yticks)
    # plt.show()
    #plt.yticks(rotation=0)
    #plt.xticks(rotation=90, size=7)
    #ax.xaxis.tick_top()
    sns.heatmap(heatmapDF_pass, mask=masked_pass_df.mask,yticklabels=yticks1,xticklabels=xticks1,annot=True, annot_kws={"size": 4},fmt='g',
                ).set(title='Pass-Based Heatmap (Plate: ' +plateID+')')
    plt.yticks(rotation=0)
    plt.xticks(rotation=90, size=7)
    ax.xaxis.tick_top()
    #plt.show()
    f1.savefig(plateID+".pass.heatmap.pdf", format="pdf", bbox_inches="tight")






#heatmapGenerator('UDI15PL013', gb)

### get unique plate IDS
unique_plateID = list(set(joinedDf['Plate Id']))
for f in unique_plateID:
    print(f)
    #print(gb)
    heatmapGenerator(f, gb)