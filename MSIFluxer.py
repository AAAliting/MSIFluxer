import pandas as pd
import os
from rpy2.robjects import r
r('source("C:/Users/Administrator/A-MFA/MSIFluxer.R")')

# read Excel
file = 'Result_corrected.xlsx'
sheet_name = 'Normalized'
df = pd.read_excel(file, sheet_name=sheet_name)

# processing
df['Name'] = df.iloc[:, 0].str.replace(r'_[0-9]+$', '', regex=True)
df['Spectrum'] = 'm' + df['C_Label'].astype(str)
df['Select'] = 1
df['Std'] = 1

# save
output_folder = 'MDV_MFA'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for Glc_Label_num in range(1, 4):
    df['MDV'] = df[f'Glc_Label_{Glc_Label_num:02d}']

    output_file = os.path.join(output_folder, f'Glc_Label_{Glc_Label_num:02d}.txt')
    df[['Name', 'Spectrum', 'Select', 'MDV', 'Std']].to_csv(output_file, sep='\t', index=False, header=True)

print("Success!")

