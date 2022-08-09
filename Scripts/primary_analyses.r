#script used for all primary analyses 
library(readxl)
library(ape)
library(dplyr)
library(geiger)
library(effsize)

# load data
dir.create("Amphibia_polyver_results")
geo <- read.csv("Amphibia_geodata.txt", sep='\t')

# load ploidy inferences
ploidy <- read_excel("canonical_polyploids.xlsx", sheet = "Amphibia", col_names = c("species", "ploidy", "source", "comments"), na = "NA")

#use more conservative polyploid inferences
ploidy[!is.na(ploidy$comments),]$ploidy <- NA

ploidy <- ploidy[1:2]
ploidy <- ploidy[ploidy$species %in% geo$species,]

# merge
df <- merge(geo, ploidy, all = TRUE)
df[!(df$species %in% ploidy$species),]$ploidy <- 0
rownames(df) <- unlist(df$species)

# remove diploids without chromosome count data
#count <- read.csv('Amphibia_count.txt', sep='\t', col.names = c('species', 'chr_count'), header=FALSE)
#diploids <- df[df$ploidy == 0,]$species
#suspect_diploids <- diploids[!(diploids %in% count$species)]
#df <- df[!(df$species %in% suspect_diploids),]

# read tree
tree <- read.tree("Amphibia_tree.nh")
tree <- multi2di(tree)
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree <- keep.tip(tree, df$species)
tree$tip.label[is.na(tree$tip.label)] <- "trim"
tree <- drop.tip(tree, c("trim"))

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

results_table <- data.frame(var = var_names, cohens.d = cohens_d, anova.p = anova_p)
write.table(results_table, 'Amphibia_polyver_results/Amphibia_table.txt', quote = FALSE, row.names = FALSE, sep = '\t')

sig_vars <- results_table[results_table$anova.p < 0.05/length(var_names),]$var
res <- Map(combn, list(sig_vars), seq_along(sig_vars), simplify = FALSE)
res <- unlist(res, recursive = FALSE)

datalist = list()

#center data for MANOVA
df <- df %>%
  mutate_at(var_names, scale, scale = FALSE)

for (vars in res) {
    manova.df <- na.omit(df)
    dat <- data.frame(manova.df[unlist(vars)])
    rownames(dat) <- unlist(manova.df$species)
    group <- factor(manova.df$ploidy)
    names(group) <- unlist(manova.df$species)
    phy.manova <- aov.phylo(dat~group, tree, nsim = 1)
    write(paste(c(paste(unlist(vars), collapse='+'), extractAIC(phy.manova)[2]), collapse='\t'),"Amphibia_polyver_results/Amphibia_models.txt", append=TRUE)}
    
random_polyploids <- sample(df$species, length(df[which(df$ploidy==1),]$species), FALSE)
random_NA <- sample(df$species, length(df[which(is.na(df$ploidy)),]$species))

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
write.table(results_table, 'Amphibia_polyver_results/Amphibia_null_table.txt', quote = FALSE, row.names = FALSE, sep = '\t')

