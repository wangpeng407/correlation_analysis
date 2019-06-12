Rscript ../corr.stat.plot.r --taxa_abund cag.abundance.matrix.xls --pheno clinical_data.xls --taxa_info map.txt --outdir output 1>1.log 2>2.log
Rscript ../corr.stat.plot.r --taxa_abund cag.abundance.matrix.xls --pheno clinical2.data.xls --taxa_info map.txt --outdir output > 11.log 2> 22.log
Rscript ../corr.stat.plot.r --taxa_abund cag.abundance.matrix.xls --pheno clinical3.data.xls --taxa_info map.txt --outdir output3 1>111.log 2>222.log &
