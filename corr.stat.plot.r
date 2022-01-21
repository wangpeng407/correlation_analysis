library(optparse)
option_list <- list(
  make_option(c("--taxa_abund"), action="store", type="character", default=NULL, help="Input the taxa abundance table, coloums are samples, row are taxa"),
  make_option(c("--clr"), action="store_true", default=FALSE, help="if transform data using center log ratio"),
  make_option(c("--pheno"), action="store", type="character", default=NULL, help="Input the pheno data or env data or clinical data"),
  make_option(c("--taxa_info"), action="store", type="character", default=NULL, help = "input the taxa_info, first colum is id and second column is map[Phylum]"),
  make_option(c("--na_percent"), action="store", default=0.1, type="double", help="if one pheno contains 10% NA, then remove it, default is 0.1"),
  make_option(c("--low_abun"), action="store", default=0.0001, type="double", help="if one taxa's average abundance < 0.0001, then remove it"),
  make_option(c("--low_freq"), action="store", default=0.1, type="double", help="if the frequency of one taxa < 0.1, then remove it"),
  make_option(c("--top_taxa"), action="store", default=100, type="integer", help="choose top 100 taxa for subsequent analysis"),
  make_option(c("--corr_method"), action="store", default='spearman', type="character", help="the method of correlation, default is spearman"),
  make_option(c("--prefix"), action="store", default=NULL, type="character", help="input the prefix of output files"),
  make_option(c("--outdir"), action="store", default="./", type="character", help="The output dirctory, default is ./"),
  make_option(c("--distance"), action="store", default="bray", type="character", help="The output dirctory, default is bray, other choices are manhattan, euclidean, canberra, clark, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup, binomial, chao, cao or mahalanobi")
)
opt <- parse_args(OptionParser(usage="%prog [options] file\n", option_list=option_list))

path_get <- function(){
    res <- list()
    allargs <- commandArgs(trailingOnly = FALSE)
    scr.name <- sub('--file=', '', allargs[grepl('--file=', allargs)])
    abs_path <- normalizePath(scr.name)                                                                                        
    res$scr.path <- dirname(abs_path)
    res$scr.name <- scr.name
    return(res)
}
subRsc <- paste(path_get()$scr.path, "src/cor.stat.R", sep="/")
if(file.exists(subRsc)){
    source(subRsc)
}else{   
    cat("NO ", subRsc, 'exists, check please!\n')
    quit()
}                  

package_list <- c("vegan","ggplot2","psych", "pheatmap", "plyr")
for(pk in package_list){
  if(!suppressWarnings(suppressMessages(require(pk, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    stop("WARNING: Please install ", pk, "!\n")
  }else{
#    cat("YES:", pk, "was succesfully installed and loaded!\n")
  } 
} 

for(fl in c(opt$taxa_abund, opt$pheno)){
  if(is.null(fl) | !file.exists(fl)){
    cat("WARNING:", fl, " not exist!\n")
    quit()
  }
}
if(is.na(file.info(opt$outdir)$isdir)) dir.create(opt$outdir, recursive = TRUE)

micro <- read.table(opt$taxa_abund, sep = '\t', header = T, row.names = 1, stringsAsFactors = F, comment.char = "", check.names = FALSE)

clr.micro <- as.data.frame(apply(micro+1, 2, clr))
#if(!is.null(opt$trans)){micro <- as.data.frame(t(micro))}

pheno <- read.table(opt$pheno, sep = '\t', header = T, row.names = 1, blank.lines.skip =  F, stringsAsFactors = F)

env.name <- colnames(pheno)
#remove vars contain 10% nas
clean_pheno <- pheno[, apply(pheno, 2, na.prop) < as.numeric(opt$na_percent)]
#numeric transform
clean_pheno <- as.data.frame(apply(clean_pheno, 2, as.numeric))
rownames(clean_pheno) <- rownames(pheno)
#fill NAs
clean_pheno <- as.data.frame(apply(clean_pheno, 2, na.replace))
#remove non-varied variables
clean_pheno <- clean_pheno[, apply(clean_pheno, 2, cv) != 0]
na10.vars <- env.name[apply(pheno, 2, na.prop) >= as.numeric(opt$na_percent)]
nonvaried.vars <- env.name[apply(clean_pheno, 2, cv) == 0]
cat("Vars contain >= 10% NA:", na10.vars, "\n\n")
cat("Vars with non-variatiton:",nonvaried.vars,"\n\n")

#seperate binary and numerical vars
btf <- find_binary_var_index(clean_pheno, 2)
all.vars <- colnames(clean_pheno)
if(any(btf)){
	binary_pheno <- clean_pheno[,btf, drop=FALSE]
	cat("Binary vars: ", all.vars[btf], "\n\n")
}else{
	binary_pheno <- NULL
}
if(any(!btf)){
	numerical_pheno <- clean_pheno[,!btf, drop=FALSE]
	cat("Numerical vars: ", all.vars[!btf], "\n\n")
}else{
	numerical_pheno <- NULL
}

#filter taxa table
clean_micro <- as.data.frame((filter_by_abun_freq(micro, low_abun = opt$low_abun, low_freq = opt$low_freq)))
if(is.integer(opt$top_taxa)){
  if(nrow(clean_micro) < opt$top_taxa){
    warning('Top taxa number must <= row number of filtered micro matrix')
    opt$top_taxa <- nrow(clean_micro)
  }
  top_micro <- as.data.frame(t(top_choose(clean_micro, margin = 1, top = opt$top_taxa)))
}else{
  top_micro <- as.data.frame(t(clean_micro))
}

if(opt$clr){
	top_micro  <- clr.micro[rownames(clr.micro) %in% colnames(top_micro), ]
	top_micro <- t(top_micro)
	cat("Perform CLR transformation...\n")
}
if(!is.null(opt$taxa_info)) {
	taxa_map <- read.table(opt$taxa_info, sep = '\t', header = F, stringsAsFactors = F)
	colnames(taxa_map) <- c('taxa', 'map')
	top_taxa_map <- taxa_map[match(colnames(top_micro), taxa_map$taxa), ]
	#get annotation colors
	annotation = data.frame(Note=factor(top_taxa_map$map))
	rownames(annotation) <- top_taxa_map$taxa
	ann_col <- getRandomColor(length(unique(annotation$Note)))
	names(ann_col) <- unique(annotation$Note)
}else{
	annotation <- ann_col <- NULL	
}
#caculate correlation matrix
#write.table("haha.xls",file=top_micro, sep="\t", quote=F, row.names = F)
if(!is.null(numerical_pheno)){
	res <- corr.test(top_micro, numerical_pheno, method = opt$corr_method, adjust = 'BH')
	rmat <- res$r
	pmat <- res$p
	opt$prefix <- ifelse(is.null(opt$prefix), '', paste0(opt$prefix, '.'))
	rout <- paste(opt$outdir, '/', opt$prefix, 'numerical.cor.r.csv', sep = '')
	pout <- paste(opt$outdir, '/', opt$prefix, 'numerical.cor.p.csv', sep = '')
	write.csv(rmat, rout, quote = F)
	write.csv(pmat, pout, quote = F)
	x.n <- ncol(rmat)
	y.n <- nrow(rmat)
	fold <- max(nchar(rownames(rmat))) / 5
	width1 <- (8 + 0.2 * x.n ) + fold/1.5
	height <- (10 + 0.2 * y.n) 
	outfig1 <- paste(opt$outdir, '/', opt$prefix, 'numerical.cor.heat.pdf', sep = '')
	pheatmap(rmat, 
			#color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         	display_numbers = matrix(ifelse(pmat>0.01&pmat<0.05, '*',ifelse(pmat<=0.01, '**', ' ')),nrow = nrow(pmat)),
	        annotation_row = annotation, border_color = NA, fontsize=15, fontsize_number=20,number_color='grey20',
	        annotation_colors = list(Note = ann_col),
		    filename = outfig1, width = width1, height = height, limitsize = FALSE)
}

if(!is.null(binary_pheno)){
	stat_mean <- apply(top_micro, 2, w.res.get, restype = 'statistics', binary_mat = binary_pheno, stat_mean = T)
	stat_mean.mat <- ldply(stat_mean, data.frame)
	names(stat_mean.mat) <- c('Taxa_id', 'Env_id', 'Type', 'Diff_mean', 'W-statistics', 'p-value', 'q-value')
	outfile1 <- paste(opt$outdir, '/', opt$prefix, 'wilcox.test.diff.mean.xls', sep = '')
	write.table(stat_mean.mat, file = outfile1, quote = F, sep = '\t', row.names = F)

	wres.list <- stat_mat_get(top_micro, binary_mat = binary_pheno)
	wmat <- as.data.frame(t(wres.list$WMAT))
	p.mat <- as.data.frame(t(wres.list$PMAT))
	x.n <- ncol(wmat)
	y.n <- nrow(wmat)
	fold <- max(nchar(rownames(rmat))) / 5
	width <- (8 + 0.2 * x.n ) + fold/1.5
	height <- (10 + 0.2 * y.n) 
	outfig3 <- paste(opt$outdir, '/', opt$prefix, 'binary.diff.heat.pdf', sep = '')
	pheatmap(wmat, scale = 'column', 
		 #color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         display_numbers = matrix(ifelse(p.mat>0.01&p.mat<0.05, '*',ifelse(p.mat<=0.01, '**', ' ')),nrow = nrow(p.mat)),
         number_color = 'black', annotation_row = annotation, border_color = 'NA',
         annotation_colors = list(Note = ann_col), fontsize=15, fontsize_number=20, number_color='grey20',
         filename = outfig3, width = width, height = height, limitsize = FALSE)
}
size = 15
if(!is.null(numerical_pheno)){
	pca <- tryCatch(rda(top_micro, scale = T))
	nmds <- tryCatch(metaMDS(top_micro, distance = opt$distance, k=2, try = 100, trace = 0))
	for(res in c('pca', 'nmds')){
		tmp <- get(res)
  		evnfit.res <- envfit(tmp, numerical_pheno)
	    sign = ifelse(evnfit.res$vectors$pvals <= 0.01, '**', ifelse(evnfit.res$vectors$pvals >0.01 & evnfit.res$vectors$pvals <=0.05, '*', ' '))
	    pd.ef <- data.frame(
		    env = names(evnfit.res$vectors$r), r2 = evnfit.res$vectors$r, 
		    pval = evnfit.res$vectors$pvals,
	    	pos = evnfit.res$vectors$r +  0.02*max(evnfit.res$vectors$r),
		    sign = sign)
	  levels(pd.ef$sign) <- c(" ", "*", "**")
	  p <- 
	    ggplot(pd.ef, aes(reorder(env, r2), r2, fill = sign, color = sign)) + 
    	geom_bar(stat = 'identity', color = 'white', width = 1, show.legend = FALSE) +
	    scale_fill_manual(values = c('grey30', 'lightseagreen', 'limegreen')) +
    	xlab('') + 
	    geom_text(aes(reorder(env, r2), pos), label = pd.ef$sign, angle = 90, size = 8, color = 'firebrick1', show.legend = FALSE)+
    	coord_flip() + 
		theme(
			 axis.text = element_text(size = size),
			 axis.title = element_text(size = size+3)
		)
	  outfile <- data.frame( env = pd.ef$env, r2 = pd.ef$r2, pval = pd.ef$pval, sign = sign)
	  outfile2 <- paste(opt$outdir, '/', opt$prefix, res, '.envfit.xls', sep = '')
	  write.table(outfile, file = outfile2, sep = '\t', quote = F, row.names = F)
	  outfig4 <- paste(opt$outdir, '/', opt$prefix, res, '.envfit.bar.pdf', sep = '') 
	  ggsave(plot = p, filename = outfig4, height = width1, width = 8)
	}
}
