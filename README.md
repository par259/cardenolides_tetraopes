For the assessment of cardenolide content using a genin vector pipeline, LC-MS data requires the following preprocessing steps. Raw LC-MS files were first converted from .raw to .mzML using ProteoWizard MSconvert v3.0 and then processed in MZmine 4.2.0 (mzio GmbH) using the "mzwizard" tool (Box 2 in Heuckeroth et al., 2024)  without parameter optimization. The overall feature table generated was exported for molecular network analysis. The resulting quantification file ("_quantif.csv") was refined by removing non-essential columns (D-M), renaming headers of column A (row ID), B (row mz), C (row retention time) with id, mz, and rt respectively. The number of decimals was reduced in mz (column B) to optimize filtering in R Markdown. The final processed dataset was saved as a .csv file in the "/data" folder within the appropriate R project environment. In R Studio, we filtered the feature table (.csv) obtained from the LC-MS analysis using a predefined cardenolide genins vector (see a summary version of the R markdown below). The use of the genins vector in the R markdown is limited to samples from Asclepias related origin (plant or insect material). Retention times that contained mz associated with genins were further investigated and compared with MS² spectra of cardenolides isolated from Asclepias syriaca . 
resources.R: a script with the following packages/libraries "tidyverse", "janitor", "ggpubr", "ggfortify", GGally, "summarySE", "doBy","scales", sjmisc, "stringr", function: `%nin%` = Negate(`%in%`)

R markdown for cardenolide pipeline
source(file = "script/ resources.R ")
genins_3 <-  read_csv("data/genins_3.csv", show_col_types = FALSE) %>%
    clean_names()
genins_3_vector= genins_3[['mz']] 
print(genins_3_vector)
tetraopes<- read_csv("data/data_iimn_gnps_quant.csv",show_col_types = FALSE) %>%
    clean_names()
data_genins_tetraopes<-tetraopes%>%
filter(mz%in%genins_3_vector)
convert_rt_tetraopes = data_genins_tetraopes [['rt']] 
print(convert_rt_tetraopes)
tetraopes_ cards<-tetraopes %>%
filter(rt%in%convert_rt_tetraopes)
tetraopes _cards%>%
  filter(between(rt,0,8))%>%
    pivot_longer( cols = -c(id,mz,rt),
      names_to = "sample",
      values_to = "abundance"
    )%>%
    select(sample,abundance)%>%
  group_by(sample)%>%
  summarise_each(funs(sum))
