Introduction
MSIFluxer is a python package which includes labeled metabolites identification and metabolic flux analysis. The package was originally designed and tested with the lab-built ambient airflow-assisted desorption electrospray ionization (AFADESI)–MSI system but can accommodate all types of MSI data as long as the input of correct format.
Please refer to the documents (https://github.com/fumiomatsuda/mfapy and https://github.com/zjgt/DietaryFluxomics.git)  for what is mfapy.
Install
1.Install 64 bit version of Anaconda3.
2.Create virtual environment such as named “mfapy” in Anaconda Prompt and install required packages.
> conda create -n mfapy python=3.9 numpy scipy matplotlib joblib rpy2
> conda activate mfapy
3.Download A-MFA.zip and unzip.
Data preparation
A folder named “A-MFA” is created on drive C of your computer containing the following files:
mydata.xlxs includes MSI data and target for labeled metabolites identification.
SC_MDV.txt and SC_status.csv for solver optimization.
Metabolic network model (MFA_Model.txt), constraint file (Status_Folder: MFA_status) and Carbon Source (MFA_Carbon_Source) for metabolic flux analysis.
MSITracer.R, MSITracer.py and MFA.py for data processing.
A Result_MFA folder for storing result files.
Data processing in Anaconda Prompt
1.Go to your data folder (e.g., A-MFA) and 
> cd C:\A-MFA
2.Run MSIFluxer.py using following code
> python MSITracer.py
3.Updata SC_MDV.txt and MFA_Model.txt
4.Run MFA.py using following code
> python MFA.py
