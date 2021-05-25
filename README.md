# PDX-ML
Code and preprocessed pharmacogenomic data to predict PDX response to drug treatments

The scripts are in the folder ./scripts 
The processed data files are in ./data
Other input file required for model execution is in ./input, including:
 - List of PDXs belonging to each cancer type(pdxe_cancertype.csv)
 - List of treatments and the corresponding number of PDXs tested for that treatment (pdxe_treatment50.csv)

To run the scripts, please be sured that the three folders are put at the same level and set the working directory to the source file location.

There are two scripts available:
1. LOOCV_RF_allvsOMC: for executing all RF models: generating 2 excel-file results
  + 2 cancer types x 13 compounds x 2 models (all vs OMC) = 52 models
  + each excel-file contains results of all models relating 1 cancer type
 The resulted excel file includes the following columns:
 - Treatment_ID: The ID of the treatment
 - Treatment_Name: The name of the treatment
 - Number of PDXs: The number of PDXs available for that drug-cancer type pair
 - Median nfeats rep 1 - 10: The median of the optimized number of features resulted from each of all LOOCV RF-OMC models of the inner loop for that experiment, 
   for each replicate 1 to 10
 - MCC.OMC.Rep1 - 10: The MCC calculated when trained and tested RF-OMC by the outer loop of nested LOOCV of that experiment, for each replicate 1 to 10.
 - MCC.OMC.med: The median MCC value of RF-OMC model across 10 replicates 
 - MCC.all.Rep1 - 10:  The MCC calculated when trained and tested RF-all using LOOCV of that experiment, for each replicate 1 to 10.
 - MCC.all.med: The median MCC value of RF-all models across 10 replicates
 - PREC.OMC.med: The median Precision value of RF-OMC models across 10 replicates
 - RECALL.OMC.med: The median Recall value of RF-OMC models across 10 replicates
 - SPEC.OMC.med: The median Specificity value of RF-OMC models across 10 replicates
 - F1.OMC.med: The median F1 value of RF-OMC models across 10 replicates
 - nResp.OMC: The median number of PDXs of that data set that were predicted sensitive by RF-OMC model across 10 replicates.

2. LOOCV_RF_allvsOMC_1case: for executing one specific case: requires users' input.
  + Users need to input the information relating to the case that they want to execute
  under ## INPUT BY USERS
      _ cancertype : can be either "BRCA" or "CRC"
      _ feat.type: can be either "SNV", "CNA", "CN" or "GEX"
      _ treatmentName: name of the drugs/ treatment. 
  + For BRCA, models are available for "BGJ398", "BKM120", "BYL719", "BYL719 + LJM716",
    "CLR457", "HDM201", "INC424", "LEE011", "LLM871", "binimetinib", "paclitaxel"
  + For CRC, models are available for "BKM120", "BYL719", "BYL719 + LJM716",
    "CLR457", "HDM201", "LEE011", "LKA136", "binimetinib", "cetuximab", "CKX620", "encorafenib",
    "LFW527 + binimetinib"
 The resulted excel file includes the following columns:
 - Replicates: The number of replicate
 - Median.nfeats: The median of the optimized number of features resulted from each of all LOOCV RF-OMC models of the inner loop for that experiment, 
   for the current replicate
 - MCC.OMC: The MCC value of the RF-OMC model for that replicate
 - PREC.OMC: The Precision value of the RF-OMC model for that replicate
 - RECALL.OMC: The Recall value of the RF-OMC model for that replicate
 - F1.OMC: The F1 value of the RF-OMC model for that replicate
 - MCC.all: The MCC value of the RF-all model for that replicate
 - PREC.all: The Precision value of the RF-all model for that replicate
 - RECALL.all: The Recall value of the RF-all model for that replicate
 - F1.all: The F1 value of the RF-all model for that replicate
