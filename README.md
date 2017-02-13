# Mussel_feeding

The full analysis of the submitted paper: ***Invasive dreissenid mussels induce shifts in bacterioplankton diversity through selective feeding on high nucleic acid bacteria*** can be found in the **Analysis.Rmd** RMarkdown. 

Before starting the analysis please unzip **16S.zip** and download the FCM data for the mussel experiment from [here](https://flowrepository.org/experiments/1034) and store them in a directory named **data_mussel**. The reference FCM data (approx. 1GB) for the cooling water system and Lake Michigan/Muskegon Lake are available from [here](https://flowrepository.org/experiments/746) and [here](https://flowrepository.org/experiments/1047). These data should be stored in a directory called **data_reference** with subdirectories **FCM_CW** for the cooling water data and **FCM_MI** for the Lake Michigan and Muskegon Lake data.

**Notice:** Some bootstrap repeats are set rather low to facilitate faster run times. For exactly reproducing the output described in the manuscript please adjust the parameters to those specified in the methods section. The output from this markdown file is provided in html format.
