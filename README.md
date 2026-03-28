# MSIFluxer
## Introduction
MSIFluxer is a python package which includes labeled metabolites identification and metabolic flux analysis. The package was originally designed and tested with the lab-built ambient airflow-assisted desorption electrospray ionization (AFADESI)–MSI system but can accommodate all types of MSI data as long as the input of correct format. The docker image ting.tar contains entire envorienment for running MSIFluxer.
## Install
* Download A-MFA.zip and unzip.
+ Install Docker.
- load image
```docker pull ghcr.io/aaaliting/ting:latest```
## Data preparation
A folder named “A-MFA” is created containing the following files:
1. mydata.xlxs includes MSI data and target for labeled metabolites identification.
2. SC_MDV.txt and SC_status.csv for solver optimization.
3. Metabolic network model (MFA_Model.txt), constraint file (Status_Folder: MFA_status) and Carbon Source (MFA_Carbon_Source) for metabolic flux analysis.
4. MSITracer.R, MSITracer.py and MFA.py for data processing.
5. A MDV_MFA folder and a Result_MFA folder are created for storing result files.
## Data processing with MSIFluxer
* Open Powershell
- Run docker using following code
```docker run -it --rm -v D:/desktop/A-MFA:/MSIFluxer ghcr.io/aaaliting/ting```

  
  
  
  
  