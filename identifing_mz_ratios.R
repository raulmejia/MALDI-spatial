#########
## This program 
## 
########
#####
# Data given by the user
# The user should enter 2 numbers and 1 character string by each desired phenoytpe
# And one input table of course:
#   Avoid whitespaces in the name of of the columns and file name
#   The recomended columns names (in that order): 
#   The output is joined table of all the subtables cut according the parameters given by the user 
#   
####
# It needs the column "
#####
# loading required libraries
#####
if (!require("BiocManager")) {
  install.packages("BiocManager", ask =FALSE)
  library("BiocManager")
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}
if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}
if (!require("limma")) {
  BiocManager::install("limma", dependencies = TRUE)
  library(limma)
}
rm(list=ls())

#########
## Functions defined by the user
#########
extract_unique_rownames_from_list_of_matrices <- function(list_of_matrices){
  vector_of_names <- vector()
  for(k in 1:length(list_of_matrices)){
    vector_of_names <- c(vector_of_names, rownames(list_of_matrices[[k]]) )
  }
  union_compund_names <- sort( unique( vector_of_names) , decreasing =  TRUE)  # here are the rownames
  return(union_compund_names)
}  

######################
## Module: Testing the consistency of the parameters given by the user  
######################
myargs <- commandArgs(trailingOnly = TRUE)

intable_path <- myargs[1] 
out_table_path <- myargs[2] 
results_path <- myargs[3] 
# Parameters in Spatial_MALDI.sh

x1_p1 <- 80
x2_p1 <- 109
y1_p1 <- 25
y2_p1 <- 56
p1_name <-"124_Ctrl"

x1_p2 <- 78
x2_p2 <- 109
y1_p2 <- 70
y2_p2 <- 101
p2_name <-"124_CQ_20_uM"

x1_p3 <- 126
x2_p3 <- 155
y1_p3 <- 25
y2_p3 <- 54
p3_name <-"127"

x1_p4 <- 124
x2_p4 <- 156
y1_p4 <- 73
y2_p4 <- 105
p4_name <-"127_CQ_20_uM"

x1_p5 <- 194
x2_p5 <- 220
y1_p5 <- 29
y2_p5 <- 57
p5_name <-"CAU"

x1_p6 <- 191
x2_p6 <- 216
y1_p6 <- 75
y2_p6 <- 104
p6_name <-"CAU_CQ_20_uM"

x1_p7 <- 236
x2_p7 <- 265
y1_p7 <- 32
y2_p7 <- 59
p7_name <-"124_2_Ctrl"

x1_p8 <- 236
x2_p8 <- 266
y1_p8 <- 76
y2_p8 <- 105
p8_name <-"124_2_CQ_20_uM"

##############
# Testing consistency of the input parameters
##############


##############
# Loading the data and building the extracted data
##############
results_path <- normalizePath(results_path); dir.create(results_path , recursive = TRUE)

intable <-read.table( intable_path , header = TRUE, sep=",") # reading the data
p1 <- intable %>% filter(Coordinate.X > x1_p1) %>% filter(Coordinate.X < x2_p1) %>% filter(Coordinate.Y > y1_p1) %>% filter(Coordinate.Y < y2_p1) # extracting subtable for phenotype 1
p2 <- intable %>% filter(Coordinate.X > x1_p2) %>% filter(Coordinate.X < x2_p2) %>% filter(Coordinate.Y > y1_p2) %>% filter(Coordinate.Y < y2_p2)
p3 <- intable %>% filter(Coordinate.X > x1_p3) %>% filter(Coordinate.X < x2_p3) %>% filter(Coordinate.Y > y1_p3) %>% filter(Coordinate.Y < y2_p3)
p4 <- intable %>% filter(Coordinate.X > x1_p4) %>% filter(Coordinate.X < x2_p4) %>% filter(Coordinate.Y > y1_p4) %>% filter(Coordinate.Y < y2_p4)
p5 <- intable %>% filter(Coordinate.X > x1_p5) %>% filter(Coordinate.X < x2_p5) %>% filter(Coordinate.Y > y1_p5) %>% filter(Coordinate.Y < y2_p5)
p6 <- intable %>% filter(Coordinate.X > x1_p6) %>% filter(Coordinate.X < x2_p6) %>% filter(Coordinate.Y > y1_p6) %>% filter(Coordinate.Y < y2_p6)
p7 <- intable %>% filter(Coordinate.X > x1_p7) %>% filter(Coordinate.X < x2_p7) %>% filter(Coordinate.Y > y1_p7) %>% filter(Coordinate.Y < y2_p7)
p8 <- intable %>% filter(Coordinate.X > x1_p8) %>% filter(Coordinate.X < x2_p8) %>% filter(Coordinate.Y > y1_p8) %>% filter(Coordinate.Y < y2_p8)

list_subtables <-list(p1,p2,p3,p4,p5,p6,p7,p8) # putting the sutables in a list
names(list_subtables) <-c(p1_name,p2_name,p3_name,p4_name,p5_name,p6_name,p7_name,p8_name) # adding the names to the previous list
table_with_phenotype_column_2 <- do.call(rbind.data.frame, list_subtables) # pasting the previous extracted subtables ( that only conintain the desired phenotype data)

table_with_phenotype_column_clean <- cbind(rownames(table_with_phenotype_column_2),table_with_phenotype_column_2) # adding a column at the beginning to mark the phenotypes
colnames(table_with_phenotype_column_clean)[1]<-"phenotype"
table_with_phenotype_column_clean$phenotype<- gsub("\\..*","",table_with_phenotype_column_clean$phenotype)

write.table(file=paste0(results_path,"/",basename(intable_path),"_cutted_and_annoted_according_phenotypes_coordinates.tsv")
            , table_with_phenotype_column_clean, sep=",", col.names = TRUE, row.names = FALSE ) # writing down the new generatad table

RI100 <- intable %>% filter( Relative.Intensity.... == 100 ) # Extracting the compounds with RelInt = 100 (the most abundant isomers)
length( unique(table_with_phenotype_column_clean$ID)) # number of pixels in the input table
length( unique( intable$Name) ) # number of different compounds
length( unique( RI100$ID)) #  pixels in the 100 relative intensity compounds 
length( unique(RI100$Name)) # Compound names in the the compounds with RelInt = 100 (the most abundant isomers)


####################
# Exploratory Analysis
####################
## Number of pixels per phenotype
pixels_perphenotyphe <- lapply(list_subtables,dim)
pixels_perphenotyphe_list <- lapply(pixels_perphenotyphe, function(x) x[[1]])
pixels_perphenotyphe_vec <- as.vector(unlist(pixels_perphenotyphe_list))
names(pixels_perphenotyphe_vec) <- names(pixels_perphenotyphe )
pixels_perphenotyphe_vec/max(pixels_perphenotyphe_vec)

ggplot(table_with_phenotype_column_clean, aes(x=as.factor(phenotype) )) +
  geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7) ) + ggtitle("Number of pixels per phenotype") + xlab("Phenotypes")

# Number of m/z compounds per phenotype
list_subtables_phenotype_columns_included <- split(table_with_phenotype_column_clean , table_with_phenotype_column_clean$phenotype ) # splitting by phenotype
list_subtables_unique_mz <- lapply(list_subtables_phenotype_columns_included, function(K){ K[ !duplicated( K$Experimental.m.z ) , ] })
unlist(lapply(list_subtables_unique_mz,function(x) dim(x)[1]))
list_subtables_unique_mz_tabled <- do.call(rbind.data.frame, list_subtables_unique_mz)
ggplot(list_subtables_unique_mz_tabled, aes(x= as.factor(phenotype) )) +
  geom_bar(color="blue", fill = rgb(0.1,0.4,0.5,0.7) ) + ggtitle("Number of m/z per phenotype") + xlab("Phenotypes") # compounds per phenotype

# compounds per pixel
listfrom_table_with_phenotype_column_2_splitted_by_id <- split(table_with_phenotype_column_clean, table_with_phenotype_column_2$ID)
listfrom_table_with_phenotype_column_2_splitted_by_id_unique_Name <- lapply(listfrom_table_with_phenotype_column_2_splitted_by_id, function(K){ K[ !duplicated( K$Name) , ] })
tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_Name <- do.call(rbind.data.frame, listfrom_table_with_phenotype_column_2_splitted_by_id_unique_Name )
ggplot( tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_Name , aes(x=as.factor(ID) )) +
  geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7) ) + ggtitle("Number of compounds per pixel") + xlab("Phenotypes") # compounds per pixel across all phenotypes

compunds_per_pixel <- unlist(lapply( listfrom_table_with_phenotype_column_2_splitted_by_id_unique_Name , function(x){dim(x)[1]}))

compunds_per_pixel_sorted <- sort(compunds_per_pixel, decreasing = TRUE)
df_compunds_per_pixel <- data.frame(  compunds_per_pixel_sorted , rep(NA , length(compunds_per_pixel_sorted)) )
rownames( df_compunds_per_pixel ) <- names( compunds_per_pixel_sorted ) ; colnames( df_compunds_per_pixel) <- c( "Comp_per_pix" , "Phenotype" )
df_compunds_per_pixel$pixelID <- rownames( df_compunds_per_pixel) 

ID_phenotype <- tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_Name[ !duplicated(tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_Name$ID), ]


for( k in df_compunds_per_pixel$pixelID ){
  df_compunds_per_pixel[k,"Phenotype"] <- ID_phenotype %>% filter(ID ==k) %>% select(phenotype)
}

theme_set(theme_classic())
# Denstity plot 
gk <- ggplot(data=df_compunds_per_pixel , aes(Comp_per_pix))
gk + geom_density(aes(fill=factor(Phenotype)), alpha=0.8) + 
  labs(title="Density plot", 
       subtitle="Different compounds per pixel per Phenotype",
       caption="Source: Slide 1",
       x="Number of compounds",
       fill="# Phenotypes")

# Boxplot 
g <- ggplot( df_compunds_per_pixel[,c(1,2)], aes(x=Phenotype, y=Comp_per_pix,color=Phenotype))
g + geom_boxplot(varwidth=T) + 
  labs(title="Different compounds per pixel per Phenotype", 
       subtitle="Box Plot",
       caption="Source: Slide 1",
       x="Phenotypes",
       y="Number of compounds")

# m/z per pixel
listfrom_table_with_phenotype_column_2_splitted_by_id <- split(table_with_phenotype_column_clean, table_with_phenotype_column_2$ID)
listfrom_table_with_phenotype_column_2_splitted_by_id_unique_mz <- lapply(listfrom_table_with_phenotype_column_2_splitted_by_id, function(K){ K[ !duplicated( K$Experimental.m.z) , ] })
tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_mz <- do.call(rbind.data.frame, listfrom_table_with_phenotype_column_2_splitted_by_id_unique_mz )
ggplot( tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_mz , aes(x=as.factor(ID) )) +
  geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7) ) + ggtitle("Number of m/z per pixel") + xlab("Phenotypes") # compounds per pixel across all phenotypes

mz_per_pixel <- unlist(lapply( listfrom_table_with_phenotype_column_2_splitted_by_id_unique_mz , function(x){dim(x)[1]}))

mz_per_pixel_sorted <- sort(mz_per_pixel, decreasing = TRUE)
df_mz_per_pixel <- data.frame(  mz_per_pixel_sorted , rep(NA , length(mz_per_pixel_sorted)) )
rownames( df_mz_per_pixel ) <- names( mz_per_pixel_sorted ) ; colnames( df_mz_per_pixel ) <- c( "mz_per_pix" , "Phenotype" )
df_mz_per_pixel$pixelID <- rownames( df_mz_per_pixel) 

ID_phenotype_mz <- tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_mz[ !duplicated(tabled_listfrom_table_with_phenotype_column_2_splitted_by_id_unique_mz$ID), ]


for( k in df_mz_per_pixel$pixelID ){
  df_mz_per_pixel[k,"Phenotype"] <- ID_phenotype_mz %>% filter(ID ==k) %>% select(phenotype)
}

theme_set(theme_classic())
# Denstity plot 
gk <- ggplot(data=df_mz_per_pixel , aes(mz_per_pix))
gk + geom_density(aes(fill=factor(Phenotype)), alpha=0.8) + 
  labs(title="Density plot", 
       subtitle="Different mz per pixel per Phenotype",
       caption="Source: Slide 1",
       x="Number of compounds",
       fill="# Phenotypes")

# Boxplot 
g <- ggplot( df_mz_per_pixel[,c(1,2)], aes(x=Phenotype, y=mz_per_pix,color=Phenotype))
g + geom_boxplot(varwidth=T) + 
  labs(title="Different mz per pixel per Phenotype", 
       subtitle="Box Plot",
       caption="Source: Slide 1",
       x="Phenotypes",
       y="Number of different m/z")











############
##  Normalization:
##  Extracting and normalizing relative intensities according to pixels
############
### Calculating the "average pixel per phenotype"
list_subtables_splitted_by_Name <- lapply( list_subtables, function(xx){split(xx,xx$Name)} ) # Split each phenotype into their "compound Names".
list_subtables_Unique_Names <- lapply( list_subtables, function(xx){ xx[ !duplicated(xx$Name), ]} ) # eliminating duplicated names inside a same phenotype

list_sums <- list()
for(uoo in 1:length(list_subtables_splitted_by_Name) ){
  for( w in names(list_subtables_splitted_by_Name[[uoo]] ) ){
    result <- lapply( list_subtables_splitted_by_Name[[uoo]] , function(ww){ sum(ww$Absolute.Intensity) })
  }
  list_sums[[uoo]] <-  result
}
names(list_sums) <- names(list_subtables_splitted_by_Name)

list_sums_phenotypes <- list() # the list of sums of Absolute concentration per phenotype
for(tax in names( list_sums) ){
  list_sums_matrix <- do.call(rbind.data.frame, list_sums[[tax]])
  list_sums_vec <- as.vector(list_sums_matrix)
  rownames(list_sums_vec) <- names(list_sums[[tax]])
  colnames(list_sums_vec) <- tax
  list_sums_phenotypes[[tax]] <- list_sums_vec
}

# Now divide each vector by the respective number of pixels
pixels_perphenotyphe_vec_MaxNorm <- pixels_perphenotyphe_vec/max(pixels_perphenotyphe_vec) # If I divide by the number of pixels the result will be very small, that´s why I multiply that by the max num of pixels reached

list_average_pix_per_pheno <- list() # list the absolute quantity in the average pixel
for( QQ in  names(list_sums_phenotypes) ){
  list_average_pix_per_pheno[[QQ]] <- as.matrix(list_sums_phenotypes[[QQ]])/pixels_perphenotyphe_vec_MaxNorm[QQ]
}

#------- Getting the proportions from the previous one
Sum_AbsInten_by_Pheno <- unlist(lapply(list_average_pix_per_pheno,sum)) 
plot( unlist(lapply(list_average_pix_per_pheno,sum))  , xlab=c("Phenotypes"), ylab="sum Absolute intensities")

list_proportions_of_absoluteQ_in_the_average_pix_per_pheno <- list() # Here we have in a list the absolute quantity in the average pixel
for( PP in  names(list_average_pix_per_pheno) ){
  mat_div_by_sum <- as.matrix( list_average_pix_per_pheno[[PP]] )/Sum_AbsInten_by_Pheno[[PP]]
  list_proportions_of_absoluteQ_in_the_average_pix_per_pheno[[ PP ]] <- mat_div_by_sum / max(mat_div_by_sum) # dividing by the max inside that vector to avoid very small numbers
} 

# Here we have in a list the absolute quantity in the average pixel using -Log to expand the differences
Minuslog_list_proportions_of_absoluteQ_in_the_average_pix_per_pheno <- list() 
for( PP in  names(list_average_pix_per_pheno) ){
  Minuslog_list_proportions_of_absoluteQ_in_the_average_pix_per_pheno[[ PP ]] <- -log(as.matrix( list_average_pix_per_pheno[[PP]] )/Sum_AbsInten_by_Pheno[[PP]])
} 

############
#### Formating results
#### Making the previous in matrices
############
# Creating the reference name vector for all compounds´ Names
vector_of_names <- vector()
for(k in 1:length(list_sums)){
  vector_of_names <- c(vector_of_names, names(list_sums[[k]]) )
}
union_compund_names <- sort( unique( vector_of_names) , decreasing =  TRUE)  # here are the rownames

#------ Getting the matrix of the "Absolute concentration" in the average pixel -------
vector_of_reference_names <- union_compund_names # this is the same for all the matrices
list_of_titrated_matrices <- list()

for( zz in 1:length(list_average_pix_per_pheno)  ){
  titrated_matrix <- matrix( rep( NA, length( vector_of_reference_names ) ) , nrow= length( vector_of_reference_names ) ) # creating an empty matrix to store the titrated vector
  colnames(titrated_matrix) <- names(list_average_pix_per_pheno)[zz] ; rownames(titrated_matrix) <- vector_of_reference_names # 
  for( k in rownames(list_average_pix_per_pheno[[ colnames(titrated_matrix) ]]) ){
    titrated_matrix[k, colnames(titrated_matrix)] <- list_average_pix_per_pheno[[ colnames(titrated_matrix) ]][k,] 
  }
  list_of_titrated_matrices[[zz]] <- titrated_matrix[order( rownames(titrated_matrix) , decreasing = TRUE),]
}

names(list_of_titrated_matrices) <- names( list_average_pix_per_pheno)
titrated_matrix <- matrix( unlist(list_of_titrated_matrices) , ncol= length(list_of_titrated_matrices), byrow = FALSE)
colnames(titrated_matrix) <- names(list_of_titrated_matrices)
rownames(titrated_matrix) <- union_compund_names

head(titrated_matrix)
heatmap( titrated_matrix[,c("124_Ctrl","124_CQ_20_uM","124_2_Ctrl","124_2_CQ_20_uM")])
write.table( file= paste0(results_path,"/",basename(intable_path),"_Relative_intesity_Normalized_by_pixel.tsv" ), titrated_matix, sep="\t", col.names = TRUE)

# --- Getting Proportions of compounds inside the relative intensity of the average pixel ----


extract_unique_rownames_from_list_of_matrices(list_proportions_of_absoluteQ_in_the_average_pix_per_pheno) # test


make_a_matrix_out_of_me <- function(list_with_matrices ){ 
  
  union_row_names <- extract_unique_rownames_from_list_of_matrices(list_with_matrices)
  a_list_of_titrated_matrices<- list()
  
  for( zz in 1:length(list_with_matrices)  ){
    a_titrated_matrix <- matrix( rep( NA, length( union_row_names ) ) , nrow= length( union_row_names ) ) # creating an empty matrix to store the titrated vector
    colnames(a_titrated_matrix) <- names(list_with_matrices)[zz] ; rownames(a_titrated_matrix) <- union_row_names # 
    for( k in rownames(list_with_matrices[[ colnames(a_titrated_matrix) ]]) ){
      a_titrated_matrix[k, colnames(a_titrated_matrix)] <- list_with_matrices[[ colnames(a_titrated_matrix) ]][k,] 
    }
    a_list_of_titrated_matrices[[zz]] <- a_titrated_matrix[order( rownames(a_titrated_matrix) , decreasing = TRUE),]
  }
  
  names(a_list_of_titrated_matrices) <- names( list_with_matrices)
  a_titrated_matrix <- matrix( unlist(a_list_of_titrated_matrices) , ncol= length(a_list_of_titrated_matrices), byrow = FALSE)
  colnames(a_titrated_matrix) <- names(a_list_of_titrated_matrices)
  rownames(a_titrated_matrix) <- union_row_names 
  return(a_titrated_matrix)
}  

matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno <- make_a_matrix_out_of_me(list_proportions_of_absoluteQ_in_the_average_pix_per_pheno)
head(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)
heatmap(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)

heatmap(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno[,c("124_Ctrl","124_CQ_20_uM","124_2_Ctrl","124_2_CQ_20_uM")])

colnames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)
rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)
grep("phosphoethanolamine",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno))

# Chafeo el clustering
#Minuslog_Matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno <- make_a_matrix_out_of_me(Minuslog_list_proportions_of_absoluteQ_in_the_average_pix_per_pheno)
#head(Minuslog_Matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)
#heatmap(Minuslog_Matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)


#write.table( file= paste0(results_path,"/",basename(intable_path),"_Relative_intesity_Normalized_by_pixel.tsv" ), titrated_matix, sep="\t", col.names = TRUE)



### Extracting relevant rows
PE_rows <- grep("phosphoethanolamine",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno))

rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)[grep("glycero-3-phosphoethanolamine",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno))]
rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)[grep("-sn-glycero-3-phosphoethanolamine",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno))]
<- rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)[grep("anoyl-sn-glycero-3-phosphoethanolamine",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno))]

grep("diac",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)[grep("-sn-glycero-3-phosphoethanolamine",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno))])

sn_PE_rows<-rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno)[grep("-sn-glycero-3-phosphoethanolamine",rownames(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno))]


"1,2-diacyl-sn-glycero-3-phosphoethanolamine"

head(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno[,c("124_Ctrl","124_CQ_20_uM","124_2_Ctrl","124_2_CQ_20_uM")][ PE_rows, ])
heatmap(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno[,c("124_Ctrl","124_CQ_20_uM","124_2_Ctrl","124_2_CQ_20_uM")][ PE_rows, ])
heatmap(matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno[,c("124_Ctrl","124_CQ_20_uM","124_2_Ctrl","124_2_CQ_20_uM")][ sn_PE_rows, ])
head(table_with_phenotype_column_clean)

matrix_proportions_of_absoluteQ_in_the_average_pix_per_pheno[,c("124_Ctrl","124_CQ_20_uM","124_2_Ctrl","124_2_CQ_20_uM")]["PE(19:0/19:0);1.2-dinonadecanoyl-sn-glycero-3-phosphoethanolamine", ]



