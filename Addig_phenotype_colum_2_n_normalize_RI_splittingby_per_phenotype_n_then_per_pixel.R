#########
## This program adds a column of phenotypes based on the X and Y coordinates given by the user
########
#####
# Data given by the user
# The user should enter 2 numbers and 1 character string by each desired phenoytpe
# And one input table of course:
#   Avoid whitespaces in the name of of the columns and file name
#   The recomended columns names (in that order): 
# The output is joined table of all the subtables cut according the parameters given by the user 
####
# It needs the column "Name" with the compound's name
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
rm(list=ls())
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
## Number of pixels
pixels_perphenotyphe <- lapply(list_subtables,dim)
pixels_perphenotyphe_list <- lapply(pixels_perphenotyphe, function(x) x[[1]])
pixels_perphenotyphe_vec <- as.vector(unlist(pixels_perphenotyphe_list))
names(pixels_perphenotyphe_vec) <- names(pixels_perphenotyphe )
pixels_perphenotyphe_vec/max(pixels_perphenotyphe_vec)

ggplot(table_with_phenotype_column_clean, aes(x=as.factor(phenotype) )) +
  geom_bar(color="blue", fill=rgb(0.1,0.4,0.5,0.7) ) + ggtitle("Number of pixels per phenotype") + xlab("Phenotypes")

head(table_with_phenotype_column_clean)

# Number of compounds 
lapply(list_subtables,dim)
list_subtables_compounds <- lapply(list_subtables,function(x) x$Name)
list_subtables_UNIQUEcompounds <- lapply(list_subtables_compounds,unique)
lapply(list_subtables_UNIQUEcompounds , length)
str(list_subtables_UNIQUEcompounds)

# compounds per phenotype
# compounds per pixel

############
##  Extracting and normalizing relative intensities according to pixels
############
normalize_intensity_by_compoundname <- function(x){sum(x$Relative.Intensity....)/dim(x)[1]} # extracting the relative intensity divided by the num of pixels 

from_matrix_2_vector_of_Compound_names_normalizes_by_Retalive_intensity <- function( amatrix){
  # This function takes a matrix and it returns a vector of Compound "Names" that extracted the relative intensity normalized by pixel,   
  list_from_amatrix_splitted_by_name <- split( amatrix, amatrix$Name ) # splitting the given matrix
  list_pheno_name_relative_int <- lapply( list_from_amatrix_splitted_by_name, normalize_intensity_by_compoundname ) # extracting the relative intensity divided by the num of pixels 
  vec_pheno_name_relative_int <- as.vector( unlist( list_pheno_name_relative_int ) ); names(vec_pheno_name_relative_int) <- names(list_pheno_name_relative_int)
  return(vec_pheno_name_relative_int)
}

list_vecs_name_relative_int <- lapply(list_subtables, from_matrix_2_vector_of_Compound_names_normalizes_by_Retalive_intensity)

# Creating the reference name vector of compounds
vector_of_names <- vector()
for(k in 1:length(list_vecs_name_relative_int_test)){
  vector_of_names <- c(vector_of_names, names(list_vecs_name_relative_int_test[[k]]) )
}
union_compund_names <- sort( unique( vector_of_names) , decreasing =  TRUE)  # here are the rownames

############## 
# From the list of intensity-pixel normalized vectors to a matrix with all the compounds and phenotypes 
##############  
vector_of_reference_names <- union_compund_names

list_of_titrated_vectors <- list()
for(zz in 1:length(list_vecs_name_relative_int)){
  titrated_vector<- matrix( rep( NA, length( vector_of_reference_names ) ) , nrow= length( vector_of_reference_names ) )
  colnames(titrated_vector) <- names(list_vecs_name_relative_int)[zz] ; rownames(titrated_vector) <- vector_of_reference_names
  for(k in vector_of_reference_names ){
    titrated_vector[k, colnames(titrated_vector)] <- list_vecs_name_relative_int[[colnames(titrated_vector)]][k]  
  }
  list_of_titrated_vectors[[zz]] <- titrated_vector[order( rownames(titrated_vector) , decreasing = TRUE),]
  which(titrated_vector[,colnames(titrated_vector)] %in% NA)
}
names(list_of_titrated_vectors) <- names( list_vecs_name_relative_int)

titrated_matix <- matrix( unlist(list_of_titrated_vectors) , ncol= length(list_of_titrated_vectors), byrow = FALSE)
colnames(titrated_matix) <- names(list_of_titrated_vectors)
rownames(titrated_matix) <- union_compund_names

# Compounds per phenotype
apply( titrated_matix,2 , function(x) { max(unlist(lapply( list_vecs_name_relative_int , length)))- sum(is.na(x)) } )
# vs 
# lapply( list_vecs_name_relative_int , length)

write.table( file= paste0(results_path,"/",basename(intable_path),"_Relative_intesity_Normalized_by_pixel.tsv" ), titrated_matix, sep="\t", col.names = TRUE)


