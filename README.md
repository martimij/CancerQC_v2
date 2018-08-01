# CancerQC
**Shiny App for visualising whole genome sequencing QC metrics of tumour samples**

Includes:  
- Sample summaries  
- QC metric plots  

Plot descriptions:  

- **AT DROPOUT**: AT sequence dropout rate
- **UNEVENNESS**: Unevenness of coverage, calculated as local coverage RSMD 
- **MAPPING RATE**: Percentage of mapped reads
- **CHIMERIC READS**: Percentage of chimeric reads  
- **DEAMINATION**: Percentage of deamination mutations (C:G to T:A)
- **GENOME COVERAGE**: Median whole genome coverage
- **COSMIC COVERAGE**: Percentage of COSMIC Cancer Genes at <30X coverage
- **DUPLICATION**: Percentage of PCR duplicates

Sample types:

- FF = "Fresh frozen" tumour tissue  
- FFPE = "Formalin-fixed paraffin embeded" tumour tissue  

Test data includes a small set of annonymised samples.
Dummy input data table available here:  
https://github.com/martimij/CancerQC_v2/blob/master/QC_dummy_table.csv




