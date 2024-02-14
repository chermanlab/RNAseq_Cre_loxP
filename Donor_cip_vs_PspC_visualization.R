setwd('/***')
name = 'Donor_cip_versus_PspC'

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(stringr)
library(tidyHeatmap)
library(openxlsx)
library(EnhancedVolcano)



#Load LFC data
Donor_cip_FC <- read.csv("2.4_Donor_cip_Gene_LFCshrink.csv")
rownames(Donor_cip_FC) <- Donor_cip_FC[,1]
Donor_cip_FC <- Donor_cip_FC[,-1]

Donor_pspC_FC <- read.csv("2.4_Donor_pspC_Gene_LFCshrink.csv")
rownames(Donor_pspC_FC) <- Donor_pspC_FC[,1]
Donor_pspC_FC <- Donor_pspC_FC[,-1]


#Filter basemean and create list of top padj genes for cip donors
Donor_cip_depth20 <- Donor_cip_FC[Donor_cip_FC$baseMean>20,]
cip_top <- tail(arrange(Donor_cip_depth20, desc(padj)), n = 50)


#Make matrices for heatmaps
Psp_vector <- c("pspA","pspB","pspC", "pspD","pspG", "pspE")

Cell_shape_vector <- c("murE","murF","murD","ftsW","murG","murC","mrdA","mreD","mreC","mreB")

Protein_insertion_export_vector <- c("secY","tatC","yidD","yidC")

SOS_vector <- c("cho", "dinB", "dinD", "dinF", "dinG", "dinI", "dinQ", "dinS", "ftsK", "hokE", "lexA", "molR", "polB", "recA", "recN", "recX", "rmuC", "ruvA", "ruvB", "sbmC", "ssb", "sulA", "umuD", "umuC", "uvrA", "uvrB", "uvrD", "ybfE", "ydjM", "yebG", "symE", "tisA", "tisB", "rmuC", "ymfE", "ymfI", "ydeO", "ydeS", "yoaB", "intE", "ogrK", "yqgC", "yhiL", "glvB", "ibpA")


#Make linear regression of all overlapped genes
Donor_cip_FC_overlap <- subset(Donor_cip_FC, rownames(Donor_cip_FC) %in% Overlapped_genes_all)
Donor_pspC_FC_overlap <- subset(Donor_pspC_FC, rownames(Donor_pspC_FC) %in% Overlapped_genes_all)

Donor_cip_pspC_FC <- merge(Donor_cip_FC_overlap, Donor_pspC_FC_overlap, by = "row.names", all.x=TRUE)
rownames(Donor_cip_pspC_FC) <- Donor_cip_pspC_FC[,1]
Donor_cip_pspC_FC <- Donor_cip_pspC_FC[,-1]

Donor_cip_pspC_FC_SigCip <- rownames(subset(Donor_cip_pspC_FC, Donor_cip_pspC_FC$padj.x < 1e-5))
Donor_cip_pspC_FC_SigPspC <- rownames(subset(Donor_cip_pspC_FC, Donor_cip_pspC_FC$padj.y < 1e-5))
Donor_cip_pspC_FC_SigBoth <- Reduce(intersect, list(Donor_cip_pspC_FC_SigCip, Donor_cip_pspC_FC_SigPspC))

model <- lm(log2FoldChange.y~log2FoldChange.x, data = Donor_cip_pspC_FC)
summary(model)

Donor_cip_FC_top <- rownames(tail(arrange(Donor_cip_pspC_FC, desc(padj.x)), n = 50))
Donor_pspC_FC_top <- rownames(tail(arrange(Donor_cip_pspC_FC, desc(padj.y)), n = 50))
Overlapped_genes_top <- Reduce(intersect, list(Donor_cip_FC_top, Donor_pspC_FC_top))

Psp_vector_gg_pspC_top <- Reduce(intersect, list(Donor_pspC_FC_top, Psp_vector))
Psp_vector_gg <- Reduce(intersect, list(Donor_cip_pspC_FC_SigBoth, Psp_vector))

Protein_insertion_export_vector_gg <- Reduce(intersect, list(Donor_cip_pspC_FC_SigBoth, Protein_insertion_export_vector))

Cell_shape_vector_gg <- Reduce(intersect, list(Donor_cip_pspC_FC_SigBoth, Cell_shape_vector))

Heat_shock_vector_gg <- Reduce(intersect, list(Donor_cip_pspC_FC_SigBoth, Heat_shock_vector))

ToxAnti_vector_gg <- Reduce(intersect, list(Donor_cip_pspC_FC_SigBoth, ToxAnti_vector))

SOS_vector <- c("cho", "dinB", "dinD", "dinF", "dinG", "dinI", "dinQ", "dinS", "ftsK", "hokE", "lexA", "molR", "polB", "recA", "recN", "recX", "rmuC", "ruvA", "ruvB", "sbmC", "ssb", "sulA", "umuD", "umuC", "uvrA", "uvrB", "uvrD", "ybfE", "ydjM", "yebG", "symE", "tisA", "tisB", "rmuC", "ymfE", "ymfI", "ydeO", "ydeS", "yoaB", "intE", "ogrK", "yqgC", "yhiL", "glvB", "ibpA")


keyvalues.color <- 
  ifelse(rownames(Donor_cip_pspC_FC) %in% Psp_vector_gg_pspC_top, "black",
      ifelse(rownames(Donor_cip_pspC_FC) %in% Psp_vector_gg, "black",
         ifelse(rownames(Donor_cip_pspC_FC) %in% Protein_insertion_export_vector_gg, "black",
                ifelse(rownames(Donor_cip_pspC_FC) %in% Cell_shape_vector_gg, "black",
#                       ifelse(rownames(Donor_cip_pspC_FC) %in% Heat_shock_vector_gg, "gray60",
                            ifelse(rownames(Donor_cip_pspC_FC) %in% ToxAnti_vector_gg, "gray60",
                                  ifelse(rownames(Donor_cip_pspC_FC) %in% Overlapped_genes_top, "black",
                                      ifelse(rownames(Donor_cip_pspC_FC) %in% SOS_vector, "black",
                                        ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_FC_top, "gray60",
                                              ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_pspC_FC_top, "gray60",
                                                      ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigBoth, "gray60",
                                                          ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigCip, "gray60",
                                                              ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigPspC, "gray60",
                                                                  "gray60")
                                                                )
                                                            )
                                                                                                  
                                                      )
                                                  )
                                              )
                                          )
                                    )
                              )
                      )
              )
         )
#  )

keyvalues.shape <- 
  ifelse(rownames(Donor_cip_pspC_FC) %in% Psp_vector_gg_pspC_top, 21,
      ifelse(rownames(Donor_cip_pspC_FC) %in% Psp_vector_gg, 21, 
         ifelse(rownames(Donor_cip_pspC_FC) %in% Protein_insertion_export_vector_gg, 21,
                ifelse(rownames(Donor_cip_pspC_FC) %in% Cell_shape_vector_gg, 21,
#                       ifelse(rownames(Donor_cip_pspC_FC) %in% Heat_shock_vector_gg, 21,
#                              ifelse(rownames(Donor_cip_pspC_FC) %in% ToxAnti_vector_gg, 21,
                                     ifelse(rownames(Donor_cip_pspC_FC) %in% cip_top_SOS, 21,
                                     ifelse(rownames(Donor_cip_pspC_FC) %in% cip_pspC_SOS, 21,
                                        ifelse(rownames(Donor_cip_pspC_FC) %in% Overlapped_genes_top, 23,
                                                ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_FC_top, 24,
                                                      ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_pspC_FC_top, 25,
                                                            ifelse(rownames(Donor_cip_pspC_FC) %in% SOS_vector, 3,
                                                              ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigBoth, 21,
                                                                    ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigCip, 21,
                                                                          ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigPspC, 21,
                                                                                  "21")
                                                                              )
                                                                          )
                                                                     )
                                                                    )
                                                                )
                                                          )
                                                  )
                                            )
                                    )
                              )
                        )
                  )
#            )
#        )

keyvalues.fill <- 
  #                  ifelse(rownames(res_pspC_matched) %in% PMF_vector, "darkred",
  ifelse(rownames(Donor_cip_pspC_FC) %in% Psp_vector_gg_pspC_top, "darkcyan",
      ifelse(rownames(Donor_cip_pspC_FC) %in% Psp_vector_gg, "darkcyan",
         ifelse(rownames(Donor_cip_pspC_FC) %in% Protein_insertion_export_vector_gg, "deeppink",
                ifelse(rownames(Donor_cip_pspC_FC) %in% Cell_shape_vector_gg, "darkblue",
#                       ifelse(rownames(Donor_cip_pspC_FC) %in% Heat_shock_vector_gg, "mediumorchid3",
#                              ifelse(rownames(Donor_cip_pspC_FC) %in% ToxAnti_vector_gg, "darkred",
                                     ifelse(rownames(Donor_cip_pspC_FC) %in% cip_top_SOS, "darkorange",
                                        ifelse(rownames(Donor_cip_pspC_FC) %in% cip_pspC_SOS, "darkorange",
                                            ifelse(rownames(Donor_cip_pspC_FC) %in% Overlapped_genes_top, "gray60",
                                                ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_FC_top, "gray60",
                                                   ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_pspC_FC_top, "gray60",
                                                        ifelse(rownames(Donor_cip_pspC_FC) %in% SOS_vector, "darkorange",
                                                          ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigBoth, "gray60",
                                                                 ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigCip, "gray60",
                                                                        ifelse(rownames(Donor_cip_pspC_FC) %in% Donor_cip_pspC_FC_SigPspC, "gray60",
                                                                                "gray60")
                                                                              )
                                                                        )
                                                        
                                                                    )
                                                          
                                                                )
                                                            )
                                                        )
                                                  )
                                                )
                                          )
                                    )
                              )
                        )
#                  )
#  )
          


colorcode <- data.frame(keyvalues.color, keyvalues.shape, keyvalues.fill)
rownames(colorcode) <- Donor_cip_pspC_FC$target.x
Donor_cip_pspC_FC_arranged <- cbind(Donor_cip_pspC_FC, colorcode)
Donor_cip_pspC_FC_arranged <- Donor_cip_pspC_FC_arranged %>% arrange(factor(keyvalues.fill, levels = c('gray60', 'azure', 'gold', 'darkorange', 'darkred', 'mediumorchid3', 'darkblue', 'deeppink', 'darkcyan')))
Donor_cip_pspC_FC_arranged <- Donor_cip_pspC_FC_arranged %>% arrange(factor(keyvalues.color, levels = c('gray60','black')))

Donor_cip_pspC_FC_arranged$keyvalues.alpha <- Donor_cip_pspC_FC_arranged$keyvalues.fill

Donor_cip_pspC_FC_arranged$keyvalues.alpha[Donor_cip_pspC_FC_arranged$keyvalues.alpha == "gray60"] <- 0.2
Donor_cip_pspC_FC_arranged$keyvalues.alpha[Donor_cip_pspC_FC_arranged$keyvalues.alpha == "gray60" | Donor_cip_pspC_FC_arranged$keyvalues.alpha == "gold" | Donor_cip_pspC_FC_arranged$keyvalues.alpha == "deeppink" | Donor_cip_pspC_FC_arranged$keyvalues.alpha == "darkblue"] <- 0.4
Donor_cip_pspC_FC_arranged$keyvalues.alpha[Donor_cip_pspC_FC_arranged$keyvalues.alpha != 0.2 & Donor_cip_pspC_FC_arranged$keyvalues.alpha != 0.4] <- 1

Donor_cip_pspC_FC_arranged$keyvalues.alpha[Donor_cip_pspC_FC_arranged$keyvalues.shape == 3] <- 1
Donor_cip_pspC_FC_arranged$keyvalues.alpha[Donor_cip_pspC_FC_arranged$keyvalues.shape == 23] <- 1



ggplot(Donor_cip_pspC_FC_arranged, aes(log2FoldChange.x, log2FoldChange.y)) + 
  geom_point(size=3.0, color = Donor_cip_pspC_FC_arranged$keyvalues.color, shape = as.numeric(Donor_cip_pspC_FC_arranged$keyvalues.shape), fill = Donor_cip_pspC_FC_arranged$keyvalues.fill, alpha = as.numeric(Donor_cip_pspC_FC_arranged$keyvalues.alpha)) + 
#  geom_smooth(method='lm', linetype = "dashed", linewidth = 0.5, color = "black", se = FALSE)   +
  xlab("+Cip   log2FoldChange") + 
  ylab("+PspC   log2FoldChange") +
  theme_bw(base_size = 20,
           base_line_size = 0.5,
           base_rect_size = 0.5) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "black", size = 20),
        axis.text.y = element_text(color = "black", size = 20)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) 


#List SOS genes that are significantly up or down in cip or pspC. Note that there is no overlap between these two lists.
cip_SOS <- Reduce(intersect, list(Donor_cip_pspC_FC_SigCip, SOS_vector))
pspC_SOS <- Reduce(intersect, list(Donor_cip_pspC_FC_SigPspC, SOS_vector))
cip_pspC_SOS <- Reduce(intersect, list(Donor_cip_pspC_FC_SigBoth, SOS_vector))
cip_top_SOS <- Reduce(intersect, list(Donor_cip_FC_top, SOS_vector))

#Save lists of genes highlighted on scatterplot to csv
write.csv(Donor_cip_pspC_FC, "Donor_cip_pspC_FC.csv", row.names=FALSE)

Donor_cip_FC_top_padj <- Donor_cip_pspC_FC[Donor_cip_FC_top,]
Donor_cip_FC_top_padj <- Donor_cip_FC_top_padj[-c(5:8)]

Donor_pspC_FC_top_padj <- Donor_cip_pspC_FC[Donor_pspC_FC_top,]
Donor_pspC_FC_top_padj <- Donor_pspC_FC_top_padj[-c(1:4)]

write.csv(Donor_cip_FC_top_padj, "Donor_cip_FC_top_padj.csv", row.names=FALSE)
write.csv(Donor_pspC_FC_top_padj, "Donor_pspC_FC_top_padj.csv", row.names=FALSE)

Donor_cip_pspC_FC_HighSig <- Donor_cip_pspC_FC[Donor_cip_pspC_FC_SigBoth,]
write.csv(Donor_cip_pspC_FC_HighSig, "Donor_cip_pspC_FC_HighSig.csv", row.names=FALSE)


Donor_cip_pspC_FC_Protein_insertion <- Donor_cip_pspC_FC[Protein_insertion_export_vector_gg,]
write.csv(Donor_cip_pspC_FC_Protein_insertion, "Donor_cip_pspC_FC_Protein_insertion.csv", row.names=FALSE)

Donor_cip_pspC_FC_Cell_wall <- Donor_cip_pspC_FC[Cell_shape_vector_gg,]
write.csv(Donor_cip_pspC_FC_Cell_wall, "Donor_cip_pspC_FC_Cell_wall.csv", row.names=FALSE)

Donor_cip_pspC_FC_Heat_shock <- Donor_cip_pspC_FC[Heat_shock_vector_gg,]
write.csv(Donor_cip_pspC_FC_Heat_shock, "Donor_cip_pspC_FC_Heat_shock.csv", row.names=FALSE)

Donor_cip_pspC_FC_ToxAnti <- Donor_cip_pspC_FC[ToxAnti_vector_gg,]
write.csv(Donor_cip_pspC_FC_ToxAnti, "Donor_cip_pspC_FC_ToxAnti.csv", row.names=FALSE)

Donor_cip_pspC_FC_SOS <- Donor_cip_pspC_FC[SOS_vector,]
Donor_cip_pspC_FC_SOS <- na.omit(Donor_cip_pspC_FC_SOS)
write.csv(Donor_cip_pspC_FC_SOS, "Donor_cip_pspC_FC_SOS.csv", row.names=FALSE)





