#Author: Joern Pezoldt
#Date: 16.02.2020
#Function:
# 1) Species Mixing QC LiquidDrop v4

#PATH
#PATH_mac
#@Romina: Set the directory according to your the path on your system
PATH_Input <- "/home/pezoldt/NAS2/pezoldt/Analysis/074_LiquidDrop_WellCode_SM/RCMs"
#Table of Reads per Sample provided by the sequencing facility
reads_per_sample <- read.delim("/home/pezoldt/NAS2/pezoldt/Data/074_LiquidDrop_WellCode_SM/074_ReadsPerSample.txt")

#Samples
sample_Ids <- c("SMcy16_10c_w1","SMcy16_10c_w2","SMcy16_10c_w3","SMcy16_10c_w4","SMcy16_10c_w5","SMcy16_10c_w6","SMcy16_10c_w7","SMcy16_onlyIBA_10cell_w8",
                "SMcy20_10c_w2","SMcy20_10c_w3","SMcy20_10c_w4","SMcy20_10c_w5","SMcy20_onlyIBA_10c_w6","SMcy20_10c_w7","SMcy20_10c_w8",
                "SMcy20_5c_w1","SMcy20_5c_w2","SMcy20_5c_w3","SMcy20_5c_w4","SMcy20_5c_w5","SMcy20_5c_w6","SMcy20_5c_w7","SMcy20_5c_w8",
                "SMcy24_5c_w1","SMcy24_5c_w2","SMcy24_5c_w3","SMcy24_5c_w4","SMcy24_5c_w5","SMcy24_5c_w6","SMcy24_5c_w7","SMcy24_5c_w8",
                "SMcy24_1c_w1","SMcy24_onlyHEK_1c_w2","SMcy24_1c_w3","SMcy24_1c_w4","SMcy24_1c_w5","SMcy24_1c_w6","SMcy24_1c_w7","SMcy24_1c_w8",
                "SMcy28_1c_w1","SMcy28_1c_w2","SMcy28_1c_w3","SMcy28_1c_w4","SMcy28_1c_w5","SMcy28_1c_w6","SMcy28_1c_w7","SMcy28_1c_w8") 
no_reads <- c("SMcy20_10c_w1")

#####
#Process BAM log
#####
#Compiles intel for each mapping statistics-log file for each of the sequencing libraries prosessed
path <- paste0(PATH_Input,"/")
files <- paste0(sample_Ids,"/BAM/Log.final.out")


library(tidyverse)
merge_out <- function (files) {
  df <- df <- read.delim(paste0(path, files[1]), header= F) %>% 
    filter(grepl("Number of input reads", V1) |
             grepl("Uniquely mapped reads", V1) |
             grepl("Uniquely mapped reads %", V1)) %>% 
    set_names("Var", "value")
}
l_results <- lapply(files, merge_out)
l_results <- lapply(l_results, function(x){
  #x = l_results[[1]]
  Var_names <- c("n_Reads_input","n_Reads_uniqueMapped","%_Reads_mapped")
  percent_x <- as.numeric(as.character(unlist(strsplit(as.character(x$value[3]),"%")[[1]][1])))
  values_x <- c(as.numeric(as.character(x[1,2])),as.numeric(as.character(x[2,2])),percent_x)
  table_x <- data.frame(Variable = Var_names, values = values_x)
})
names(l_results) <- sample_Ids



#####
#Process Reads
#####
#Compiles the information from the RCMs per library and store ane table for each RCM in a list
Cell_number <- sample_Ids
#Store UMI tables in list
l_UMI_tables <- list()
for(i in 1:length(sample_Ids)){
  #i = 1
  sample_Id_i <- sample_Ids[i]
  print(sample_Id_i)
  t_UMI_i <- read.table(paste(PATH_Input,"/",sample_Id_i,"/UMI_Matrix_",sample_Id_i,".txt/","output.dge.umis.detailed.txt",sep=""),
                      nrows = 112482, header = TRUE)
  print(head(t_UMI_i))
  #eliminate Unknown_Barcode column
  t_UMI_i <- t_UMI_i[ , !(names(t_UMI_i) %in% c("Unknown_Barcode"))]
  #New column names
  colnames(t_UMI_i)[3:ncol(t_UMI_i)] <- paste(colnames(t_UMI_i[3:ncol(t_UMI_i)]),Cell_number[i],sep="__")
  print(head(t_UMI_i))
  l_UMI_tables[[i]] <- t_UMI_i
}
names(l_UMI_tables) <- sample_Ids

#Calculate number of reads/UMIs per barcode and thus the respective species
# Make empty matrix to store 
t_UMI_contamination <- matrix(nrow = 0, ncol = 2)
colnames(t_UMI_contamination) <- c("UMI_mouse","UMI_human")
t_nGene_contamination <- matrix(nrow = 0, ncol = 2)
colnames(t_nGene_contamination) <- c("nGene_mouse","nGene_human")

for(i in 1:length(l_UMI_tables)){
  #i = 1
  t_UMI_i <- l_UMI_tables[[i]]
  sample_i <- Cell_number[i]
  print(sample_i)
  head(t_UMI_i)
  #ratios for mouse
  t_UMI_mouse_i <- t_UMI_i[grep("ENSMUS",t_UMI_i$Gene_id),]
  head(t_UMI_mouse_i)
  #UMI
  mouse_UMI_mouse_i <- sum(t_UMI_mouse_i[,grep("Mouse",colnames(t_UMI_mouse_i))])
  human_UMI_mouse_i <- sum(t_UMI_mouse_i[,grep("Human",colnames(t_UMI_mouse_i))])
  mouse_UMI <- c(mouse_UMI_mouse_i,human_UMI_mouse_i)
  #Gene Number
  mouse_nGene_mouse_i <- sum(t_UMI_mouse_i[,grep("Mouse",colnames(t_UMI_mouse_i))] > 0)
  human_nGene_mouse_i <- sum(t_UMI_mouse_i[,grep("Human",colnames(t_UMI_mouse_i))] > 0)
  mouse_nGene <- c(mouse_nGene_mouse_i,human_nGene_mouse_i)
  
  
  t_UMI_human_i <- t_UMI_i[grep("ENSG",t_UMI_i$Gene_id),]
  #UMI
  human_UMI_human_i <- sum(t_UMI_human_i[,grep("Human",colnames(t_UMI_human_i))])
  mouse_UMI_human_i <- sum(t_UMI_human_i[,grep("Mouse",colnames(t_UMI_human_i))])
  human_UMI <- c(mouse_UMI_human_i,human_UMI_human_i)
  #Gene Number
  human_nGene_human_i <- sum(t_UMI_human_i[,grep("Human",colnames(t_UMI_human_i))] > 0)
  mouse_nGene_human_i <- sum(t_UMI_human_i[,grep("Mouse",colnames(t_UMI_human_i))] > 0)
  human_nGene <- c(mouse_nGene_human_i,human_nGene_human_i)
  
  #rbind nGene
  temp_t_nGene <- rbind(mouse_nGene,human_nGene)
  colnames(temp_t_nGene) <- c("nGene_mouse","nGene_human")
  rownames(temp_t_nGene) <- paste(rownames(temp_t_nGene), sample_i, sep="__")
  t_nGene_contamination <- rbind(t_nGene_contamination,temp_t_nGene)
  
  #rbing UMI
  temp_t_UMI <- rbind(mouse_UMI,human_UMI)
  colnames(temp_t_UMI) <- c("UMI_mouse","UMI_human")
  rownames(temp_t_UMI) <- paste(rownames(temp_t_UMI), sample_i, sep="__")
  t_UMI_contamination <- rbind(t_UMI_contamination,temp_t_UMI)
}

#Genes detected human MDA cells compare detergents
mouse_nGene <- t_nGene_contamination[grep("mouse",rownames(t_nGene_contamination)),]
human_nGene <- t_nGene_contamination[grep("human",rownames(t_nGene_contamination)),]

#Select by UMI
mouse_UMI <- t_UMI_contamination[grep("mouse",rownames(t_UMI_contamination)),]
human_UMI <- t_UMI_contamination[grep("human",rownames(t_UMI_contamination)),]

#Rough Plot
plot(t_nGene_contamination[,1],t_nGene_contamination[,2])
plot(t_UMI_contamination[,1],t_UMI_contamination[,2])

#####
#Build Pivot table
#####
t_pivot <- cbind(t_nGene_contamination,t_UMI_contamination)
#Add general ID
ID_general <- unlist(lapply(strsplit(rownames(t_pivot), "__"), function(x){x[[2]]}))
t_pivot <- data.frame(t_pivot, ID_general)
#Number UMIs
l_split_by_ID <- split(t_pivot, t_pivot$ID_general)
l_split_by_ID <- l_split_by_ID[sample_Ids]
number_UMIs <- unlist(lapply(l_split_by_ID, function(x){
                                        n_UMI_i <- x$UMI_mouse[1] + x$UMI_mouse[2] + x$UMI_human[1]  + x$UMI_human[2]
                                        c(n_UMI_i, n_UMI_i)
                                        }))
t_pivot <- cbind(t_pivot, number_UMIs)
#add organism
ID_organism <- rep(c("mouse", "human"), nrow(t_pivot)/2)
t_pivot <- data.frame(t_pivot, ID_organism)
#Add sequenced & mapped reads
l_values <- lapply(l_results, function(x){
  values_x <- t(x$value)
})
reads_seq_map_per <- do.call(rbind, l_values)
rownames(reads_seq_map_per) <- sample_Ids
colnames(reads_seq_map_per) <- c("read_n_seq","read_n_map","Percent_mapped_Read")
reads_seq_map_per <- reads_seq_map_per[rep(seq_len(nrow(reads_seq_map_per)), each = 2), ]
t_pivot <- data.frame(t_pivot, reads_seq_map_per)

#Add cell_number
cell_number <- unlist(lapply(strsplit(as.character(t_pivot$ID_general),"_"), function(x){x}[2]))
#@Romina: some of the samples only contain on species
cell_number[cell_number == "onlyIBA"] <- "10c_onlyIBA" 
cell_number[cell_number == "onlyHEK"] <- "10c_onlyHEK" 
t_pivot <- data.frame(t_pivot, cell_number)
#Add Percent UMI of Read_n
t_pivot$Percent_UMI_among_mapped_reads <- (t_pivot$number_UMIs / t_pivot$read_n_seq) * 100

#####
#Plot
#####
#Species Mixing ------------
theme_set(theme_bw())
ggplot(t_pivot, aes(x=UMI_mouse, y=UMI_human)) + 
  geom_point(aes(col=cell_number),size=3) +   # draw points
  xlim(c(0, 20000)) + 
  ylim(c(0, 180000)) +  
  scale_color_manual(values = c("black","red","orange","deepskyblue","green")) +
  labs(y="n_UMI_human", 
       x="n_UMI_mouse", 
       title="Species Mixing, LiD_v4")

#Reads sequenced/mapped/UMI-----------
t_pivot_half <- subset(t_pivot, ID_organism == "mouse")
ggplot(t_pivot_half, aes(x=read_n_seq, y=read_n_map)) + 
  geom_point(aes(size=Percent_mapped_Read),shape=1) +   # draw points
  xlim(c(0, 300000)) + 
  ylim(c(0, 300000)) +   # draw smoothing line
  #scale_color_manual(values = c("black","red","orange","deepskyblue","green")) +
  labs(y="[n] reads mapped", 
       x="[n] reads sequenced", 
       title="Species Mixing, LiD_v4")

t_pivot_half <- subset(t_pivot, ID_organism == "mouse")
ggplot(t_pivot_half, aes(x=read_n_seq, y=UMI_mouse)) + 
  geom_point(aes(col=cell_number),size=4) +   # draw points
  xlim(c(0, 300000)) + 
  ylim(c(0, 20000)) +   # draw smoothing line
  scale_color_manual(values = c("black","red","orange","deepskyblue","green")) +
  labs(y="[n] UMIs mouse", 
       x="[n] reads sequenced", 
       title="Species Mixing, LiD_v4")


t_pivot_half <- subset(t_pivot, ID_organism == "mouse")
ggplot(t_pivot_half, aes(x=UMI_mouse, y=nGene_mouse)) + 
  geom_point(aes(col=cell_number),size=4) +   # draw points
  xlim(c(0, 20000)) + 
  ylim(c(0, 5000)) +   # draw smoothing line
  scale_color_manual(values = c("black","red","orange","deepskyblue","green")) +
  labs(y="[n] Genes mouse", 
       x="[n] UMIs mouse", 
       title="Species Mixing, LiD_v4")

t_pivot_half <- subset(t_pivot, ID_organism == "human")
ggplot(t_pivot_half, aes(x=UMI_human, y=nGene_human)) + 
  geom_point(aes(col=cell_number),size=4) +   # draw points
  xlim(c(0, 200000)) + 
  ylim(c(0, 10000)) +   # draw smoothing line
  scale_color_manual(values = c("black","red","orange","deepskyblue","green")) +
  labs(y="[n] Genes human", 
       x="[n] UMIs human", 
       title="Species Mixing, LiD_v4")
