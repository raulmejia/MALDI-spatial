rangeExpMz <- function( somedf){ 
  range(somedf[,"Experimental.m.z"])
}


Slide1 <- read.csv(file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/Slide1/raw_summary_table.csv",
            header=TRUE)

Slide1_full_summary <- read.csv( file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/Slide1/full_sumary_table.csv",
                   header=TRUE)

Slide2 <- read.csv( file="/media/rmejia/mountme88/Projects/Phosoholipidosis/Proteomic-MALDI/Data/Slide_2/Full_summary_table_slide_6.csv",
                                 header=TRUE)





rangeExpMz(Slide1)
rangeExpMz(Slide1_full_summary)
rangeExpMz(Slide2)
