import pandas as pd
# requires openpyxl

dataset = pd.read_excel('2022-06-27_Soy_PPO12_5samples_NanoporeRun.xlsx',index_col=0)

samplenames = dataset['Sample Name']
indexID = dataset['Index ID']
flowcelID = list(set(dataset['FlowCell ID']))[0]




# newcode list will save the index ids
newcode = []
for f in indexID:
    tempindex = f.split(sep='_')[1]
    # replace NB with barcode
    tempindex = tempindex.replace('NB','barcode')
    newcode.append(tempindex)



newdf = pd.DataFrame()
newdf['sample'] = samplenames
newdf['barcode'] = newcode

newdf.to_csv(flowcelID+'_names.txt',sep='\t',index=False,header=False)

