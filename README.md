# Mussel_feeding

The full analysis of the submitted paper: ***Invasive dreissenid mussels induce shifts in bacterioplankton diversity through selective feeding on high nucleic acid bacteria*** can be found in the **Analysis.Rmd** RMarkdown. 

Before starting the analysis please unzip **16S.zip** and download the FCM data for the mussel experiment from [here](https://flowrepository.org/experiments/1034) and store them in a directory named **data_mussel**. The reference FCM data (approx. 1GB) for the cooling water system and Lake Michigan/Muskegon Lake are available from [here](https://flowrepository.org/experiments/746) and [here](https://flowrepository.org/experiments/1047). These data should be stored in a directory called **data_reference** with subdirectories **FCM_CW** for the cooling water data and **FCM_MI** for the Lake Michigan and Muskegon Lake data.

File structure should be 

.
├── Analysis.Rmd
├── Analysis.html
├── Analysis.md
├── functions.R
├── /extra_scripts
├── /files
├── /data_reference
│   ├── **/FCM_CW**
│   └── **/FCM_MI**
└── **/16S**

The phenotypic and taxonomic diversities (`REF_diversity.csv`, `Lakes_diversity16S_F.csv` and `Lakes_diversityFCM_F.csv`) were calculated by means of the external scripts in `/extra_scripts` as the computation can take some time.

**Notice:** The bootstrap repeats are set rather low (`R = 3`) to facilitate faster run times. For exactly reproducing the output described in the manuscript please adjust the parameters to those specified in the methods section (`R = 100`). The output from this markdown file is provided in html format.
