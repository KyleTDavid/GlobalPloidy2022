#script used for all primary analyses 

#dependencies
library(readxl)
library(ape)
library(dplyr)
library(geiger)
library(effsize)

#create directory to hold results
dir.create("Amphibia_results")

#read environmental data table (for example: HERE)
geo <- read.csv("Amphibia_geodata.txt", sep='\t')

#load ploidy inferences (Table S3 in manuscript)
ploidy <- read_excel("Table_S3.xlsx", sheet = "Amphibia", col_names = c("species", "ploidy", "source", "comments"), na = "NA")

#alternative ploidy assignment using a more strict definition of polyploidy (see manuscript for details)
#ploidy[!is.na(ploidy$comments),]$ploidy <- NA

#remove unnecessary columns
ploidy <- ploidy[1:2]

#merge environmental data table with ploidy assignments
df <- merge(geo, ploidy, all = TRUE)
df[!(df$species %in% ploidy$species),]$ploidy <- 0
rownames(df) <- unlist(df$species)

#alternative ploidy assignment using only species with confirmed chromosome counts (see manuscript for details)
#count <- read.csv('Amphibia_count.txt', sep='\t', col.names = c('species', 'chr_count'), header=FALSE)
#diploids <- df[df$ploidy == 0,]$species
#suspect_diploids <- diploids[!(diploids %in% count$species)]
#df <- df[!(df$species %in% suspect_diploids),]

# read tree
tree <- read.tree("Amphibia_tree.nh")
tree$tip.label <- gsub("_", " ", tree$tip.label)

#for every environmental variable, calculate Cohen's D and perform phylogenetic ANOVA 
var_names <- colnames(df[c(2:15)])

cohens_d <- c()
anova_p <- c()

for (var in var_names) {
    
x <- cohen.d(na.omit(df[df$ploidy==1,][[var]]), na.omit(df[df$ploidy==0,][[var]]))
cohens_d <- c(cohens_d, x$estimate)

y <- df[[var]]
names(y) <- unlist(df$species)
group <- factor(df$ploidy)
names(group) <- unlist(df$species)
test <- aov.phylo(y ~ group, keep.tip(tree, df[!(is.na(df[[var]])),]$species), nsim = 1000)
anova_p <- c(anova_p, summary(test)$coefficients[2,4])
    
    }

#grab summary statistics and write output file
results_table <- data.frame(var = var_names, cohens.d = cohens_d, anova.p = anova_p)
write.table(results_table, 'Amphibia_results/Amphibia_table.txt', quote = FALSE, row.names = FALSE, sep = '\t')

#get every combination of variables found significant following Bonferroni correcion
sig_vars <- results_table[results_table$anova.p < 0.05/length(var_names),]$var
res <- Map(combn, list(sig_vars), seq_along(sig_vars), simplify = FALSE)
res <- unlist(res, recursive = FALSE)

#center data for MANOVA
df <- df %>%
  mutate_at(var_names, scale, scale = FALSE)

#perform MANOVA on every combination of significant variables and report summary AIC scores
for (vars in res) {
    manova.df <- na.omit(df)
    dat <- data.frame(manova.df[unlist(vars)])
    rownames(dat) <- unlist(manova.df$species)
    group <- factor(manova.df$ploidy)
    names(group) <- unlist(manova.df$species)
    phy.manova <- aov.phylo(dat~group, tree, nsim = 1)
    write(paste(c(paste(unlist(vars), collapse='+'), extractAIC(phy.manova)[2]), collapse='\t'),"Amphibia_results/Amphibia_models.txt", append=TRUE)}

#assign ploidy randomly by sampling without replacement from empirical inferences
random_polyploids <- sample(df$species, length(df[which(df$ploidy==1),]$species), FALSE)
random_NA <- sample(df$species, length(df[which(is.na(df$ploidy)),]$species))

#run same analysis as above on null dataset
df$dummy_ploidy <- 0
df[which(df$species %in% random_polyploids),]$dummy_ploidy <- 1
df[which(df$species %in% random_NA),]$dummy_ploidy <- NA

cohens_d <- c()
anova_p <- c()

for (var in var_names) {
    
x <- cohen.d(na.omit(df[df$dummy_ploidy==1,][[var]]), na.omit(df[df$dummy_ploidy==0,][[var]]))
cohens_d <- c(cohens_d, x$estimate)

y <- df[[var]]
names(y) <- unlist(df$species)
group <- factor(df$dummy_ploidy)
names(group) <- unlist(df$species)
test <- aov.phylo(y ~ group, keep.tip(tree, df[!(is.na(df[[var]])),]$species), nsim = 1000)
anova_p <- c(anova_p, summary(test)$coefficients[2,4])
    
    }

results_table <- data.frame(var = var_names, cohens.d = cohens_d, anova.p = anova_p)
write.table(results_table, 'Amphibia_results/Amphibia_null_table.txt', quote = FALSE, row.names = FALSE, sep = '\t')

