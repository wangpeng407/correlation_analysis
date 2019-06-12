- **Description**

  - This script is used for correlation analysis between taxa abundance and phenotypic data, and visulization.
  
  - There are two types of environmental, clinical or phenotypic data, numerical and binary variables.
  
  - For numerical variables, simply correlation such as spearman or pearson coefficient was caculated.
  
  - For binary variables, they can be as factors or groups for **wilcox test**, then the W-statistics and p-value can be acquired and visualized.
  
- **Dependencies**

  - "vegan","ggplot2","psych", "pheatmap", "plyr", "optparse"
  
- **Example**

  See ***example*** file

    ```
    Rscript corr.stat.plot.r --help
    ```
    
    ```
    Rscript corr.stat.plot.r --taxa_abund cag.abundance.matrix.xls --pheno clinical_data.xls --taxa_info map.txt --outdir output 1>1.log 2>2.log
    ```
