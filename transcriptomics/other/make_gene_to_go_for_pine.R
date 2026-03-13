## make gene_to_go.txt for pine transcriptome annotations

# generated annotation for pine transcriptome in EggNog mapper
# out put = xlsx files with many results
# converted to csv in Microsoft Excel

# columns of interest for GO analysis 
# query = "ProteinID"
# GOs = "GOTerm"

# have two output files as needed to split fastas for eggnog
# need to extract the columns from each and concatonate

# load in files

df1 <- read.csv("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Bioassay/EggNog/Part1_results/out.emapper.annotations1.csv",
                header = T)

df2 <- read.csv("~/Desktop/Stomentosus_Zn_Transcriptomics/Suitomnscriptome_8_PINE/Bioassay/EggNog/Part2_results/out.emapper.annotations2.csv",
                header = T)

# extract only query and GOs columns
df1_ex <- df1[c("query", "GOs")]
df2_ex <- df2[c("query", "GOs")]

# join 
df <- rbind(df1_ex, df2_ex)

# there are 8445 observations (genes)
# many of these do not have GO annotations

noAnnotation <- df[df$GOs=="-",]
# there are 3760 with no associated GO term

annotated <- df[!df$GOs=="-",]
# there are 4685 with GO annotations

# rename the columns in annotated
colnames(annotated) <- c("ProteinID", "GOTerm")


# this is the GO mapping file - write it
write.table(annotated, file = "~/Desktop/Stomentosus_Zn_Transcriptomics/Pinuscontorta_gene_to_go.txt", sep = "\t",
            row.names = FALSE)
