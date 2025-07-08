.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
library(ggtree)
library(tidyverse)
library(ape)
library(jsonlite)
library(GenomicRanges)
library(taxize)
library(RColorBrewer)
library(ggnewscale)
library(Biostrings)
############################### Read in log files
directory <- "/g/data/te53/zc4688/honours/analyses/chordates/new"

# List all JSON files in the directory
json_files <- list.files(directory, pattern = "\\.json$", full.names = TRUE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each file and read the JSON into a data frame
for (file in json_files) {
  tryCatch({
    # Read the file as plain text
    json_text <- paste(readLines(file, warn = FALSE), collapse = "")
    
    # Replace invalid `NaN` with `null`
    json_text <- gsub("\\bNaN\\b", "null", json_text)
    
    # Parse the JSON
    data <- fromJSON(json_text, flatten = TRUE)
    
    # Ensure data is a data frame with one row
    if (!is.data.frame(data)) {
      data <- as.data.frame(t(unlist(data)), stringsAsFactors = FALSE)
    }
    
    data_list[[file]] <- data
  })
}


# Combine all data frames into one large data frame
combined_df <- bind_rows(data_list)


#Remove unsuccessful assemblies, select only required cols and rename. Convert total and unit lengths to numeric and calculate approx copy number.
combined_df <- combined_df %>% 
  filter(is.na(`Errors.Alignment filtering`) & is.na(`Errors.Morph identification`) & is.na(`Errors.Non-ambiguous morph identification`)) %>% 
  filter(is.na(`Errors.Primary alignments`)) %>% 
  dplyr::select(`Sample_details.Sample ID`, `Sample_details.Species`, `Sample_details.TaxID`, 
                `rDNA_details.Unit length`, `rDNA_details.Total length`, `rDNA_details.Number of contigs`,
                `rDNA_details.Minimum unit length`, `rDNA_details.Maximum unit length`, `rDNA_details.Number of morphs`) %>% 
  dplyr::rename(
    `Sample ID` = `Sample_details.Sample ID`,
    Species = `Sample_details.Species`,
    TaxID = `Sample_details.TaxID`,
    `Unit length` = `rDNA_details.Unit length`,
    `Total length` = `rDNA_details.Total length`,
    `Contigs` = `rDNA_details.Number of contigs`,
    `Minimum unit length` = `rDNA_details.Minimum unit length`,
    `Maximum unit length` = `rDNA_details.Maximum unit length`
  ) %>% 
  mutate(`Total length` = as.numeric(`Total length`)) %>% 
  mutate(`Unit length` = as.numeric(`Unit length`)) %>% 
  mutate(`Minimum unit length` = as.numeric(`Minimum unit length`)) %>% 
  mutate(`Maximum unit length` = as.numeric(`Maximum unit length`)) %>% 
  mutate(CN = round(`Total length`/`Unit length`, digits=2))


stats <- read.delim("zc4688/honours/analyses/chordates/new/refmorph.stats.txt", sep = "\t", header = T)

stats <- stats %>% mutate(`Sample ID` = str_replace(file, ".rDNA.morphs.fasta", ""))

combined_df <- left_join(combined_df, stats)

rm(stats)
####################Inspect fails

allspecies <-read.tree("/g/data/te53/zc4688/honours/trees/allspecies.phy")
allspecies$edge.length <- NULL

allspeciestree <- combined_df %>% 
  mutate(success = ifelse(is.na(`Errors.Alignment filtering`) & is.na(`Errors.Morph identification`) & is.na(`Errors.Non-ambiguous morph identification`) & is.na(`Errors.Primary alignments`), "yes", "no")) %>% 
  select(`Sample_details.Species`, success) %>% 
  dplyr::rename("label" = "Sample_details.Species")

allspecies$tip.label <- str_replace_all(allspecies$tip.label, " ", "_")
allspecies$tip.label <- str_replace_all(allspecies$tip.label, "'", "")

common_organisms <- intersect(allspeciestree$label, allspecies$tip.label)

allspecies <- ape::drop.tip(allspecies, setdiff(allspecies$tip.label, common_organisms))

allspeciestree <- allspeciestree %>% filter(label %in% common_organisms) %>% left_join(classifications %>% dplyr::rename("label" = "Species") %>% mutate(label = str_replace_all(label, " ", "_")))

ggtree(allspecies) %<+% allspeciestree +
  geom_tiplab(size=1, aes(color = success)) +
  scale_color_manual(values = c("yes" = "black", "no" = "red")) + 
  new_scale_color() + 
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group", na.translate = FALSE)

  
rm(allspecies)
rm(allspeciestree)

############################### Create 'lengths' df - metadata df containing lengths + taxonomy info for all successful species with duplicates removed

lengths <- combined_df %>% 
  group_by(Species) %>%
  dplyr::summarize(across(everything(), ~ dplyr::last(.))) %>% 
  filter(!is.na(`Unit length`))

#classifications <- classification(lengths$label, db="ncbi")
#classifications <- cbind(classifications)
#classifications <- classifications %>% 
#select(kingdom, phylum, subphylum, superclass, class, subclass, superorder, order, suborder, family, subfamily, genus, species, query)
#write.table(classifications, "zc4688/honours/analyses/classifications.tsv", sep = "\t", quote = F, row.names = F)
classifications <- read.table("zc4688/honours/analyses/classifications.tsv", sep = "\t", header = TRUE)


classifications <- classifications %>% 
  dplyr::select(class, order, species) %>% 
  dplyr::rename("Species" = "species")  

classifications <- classifications %>% 
  group_by(class, order, Species) %>% summarize(across(everything(), ~ last(.))) #Take refseq if both refseq and genbank are available

lengths <- left_join(lengths, classifications %>% mutate(Species = str_replace_all(Species, " ", "_")), by="Species") %>% 
  mutate(Species = str_replace_all(Species, "_", " ")) %>%
  mutate(
    class = case_when(
      Species %in% c("Arripis georgianus", "Cyprinodon nevadensis mionectes", "Gymnocypris eckloni scoliostomus", "Distoechodon macrophthalmus",
                   "Pempheris klunzingeri", "Neosynchiropus ocellatus", "Anguilla bicolor pacifica", 
                   "Brachymystax lenok tsinlingensis", "Oncorhynchus clarkii lewisi", "Oncorhynchus masou masou", "Tachysurus vachellii") ~ "Actinopteri",
      Species %in% c("Bubalus carabanensis", "Canis lupus dingo", "Diceros bicornis minor", "Elephas maximus indicus", 
                   "Hippopotamus amphibius kiboko", "Lagenorhynchus acutus", "Ovis ammon polii", 
                   "Perognathus longimembris pacificus", "Gorilla gorilla gorilla", "Glossophaga mutica", "Molossus nigricans",
                   "Rangifer tarandus platyrhynchus", "Giraffa camelopardalis rothschildi", "mus musculus") ~ "Mammalia",
      Species %in% c("Elgaria multicarinata webbii", "Coluber constrictor foxii", "Spondylurus nitidus") ~ "Lepidosauria",
      Species %in% c("Melospiza melodia melodia", "Motacilla alba alba", "Pithys albifrons albifrons", "Ammospiza maritima maritima") ~ "Aves",
      Species %in% c("Mixophyes fleayi") ~ "Amphibia",
      TRUE ~ class  # Keep existing values if no match
    ),
    order = case_when(
      Species %in% c("Chrysemys picta bellii", "Malaclemys terrapin pileata", "Macrochelys suwanniensis", "Mauremys mutica") ~ "Testudines",
      Species %in% c("Latimeria chalumnae") ~ "Coelacanthiformes",
      Species %in% c("Neoceratodus forsteri") ~ "Ceratodontiformes",
      TRUE ~ order  # Keep existing values if no match
    ),
    class = case_when(
      order %in% c("Coelacanthiformes", "Ceratodontiformes") ~ "Sarcopterygii",
      order %in% c("Crocodylia", "Testudines") ~ "Testudines and Crocodylia",
      TRUE ~ class
    ),
    group = case_when(
      class %in% c("Actinopteri", "Cladistia") ~ "Actinopterygii",
      class %in% c("Ascidiacea", "Appendicularia") ~"Tunicata",
      class %in% c("Myxini", "Hyperoartia") ~ "Cyclostomata",
      TRUE ~ class
    ),
   Species = case_when(Species == "Branchiostoma floridae x Branchiostoma belcheri" ~ "B. floridae x B. belcheri",
                                 Species == "Branchiostoma floridae x Branchiostoma japonicum" ~ "B. floridae x B. japonicum",
                                 TRUE ~ Species #names too long for tree
   ),
  )

classifications <- lengths %>% select(Species, group, `Sample ID`) 


metadata <- lengths

lengths <- lengths %>% filter(`Unit length` < 120000)


groups <- unique(classifications$group)
colours <- setNames(brewer.pal(n = 11, name = "Paired"), groups)




########################################### Read in and prepare species tree

tree <- read.tree("/g/data/te53/zc4688/honours/trees/species_ncbi.phy")


#only needed for timetree names to save some of the failed tips
#treenames <- as.data.frame(tree$tip.label)
#colnames(treenames) <- c("Timetreename")
#replacementnames <- read.delim("zc4688/honours/trees/speciestimetree.tsv", header = T, sep = ",")

#replacementnames <- replacementnames %>% mutate(Species = str_replace_all(Species, " ", "_"), Timetreename = str_replace_all(Replacement, " ", "_") ) %>% select(-Replacement)

#treenames <- treenames %>% left_join(replacementnames) %>% mutate(Species = ifelse(is.na(Species), Timetreename, Species)) %>% select(Species)
#tree$tip.label <- as.character(treenames$Species)

tree$edge.length <- NULL
tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")
tree$tip.label <- str_replace_all(tree$tip.label, "'", "")

tree$tip.label <- ifelse(
  tree$tip.label == "Branchiostoma floridae x Branchiostoma belcheri", 
  "B. floridae x B. belcheri",
  ifelse(
    tree$tip.label == "Branchiostoma floridae x Branchiostoma japonicum",
    "B. floridae x B. japonicum",
    tree$tip.label
  )
)



common_organisms <- intersect(lengths$Species, tree$tip.label)


tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))


tree$tip.label <- str_replace_all(tree$tip.label, "'", "")


########################################### Summary tree
plotdata <- lengths %>% 
  dplyr::select(Species, `Unit length`, `CN`, `Total length`,  `Sample ID`, Contigs, group)%>% 
  dplyr::rename("label" = "Species") %>%  
  relocate(label, .before = everything())  %>% 
  filter(label %in% common_organisms)

ggplot(plotdata %>% filter(group %in% common_class)) + geom_density_ridges(aes(x = `Unit length`, fill = group, y = group), alpha = 0.8) + scale_fill_manual(values = colours, name = "Taxonomic group") + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
ggsave("zc4688/honours/analyses/chordates/figures/lengths.ridges.pdf",)
ggplot(plotdata %>% filter(group %in% common_class)) + geom_density_ridges(aes(x = `Unit length`, fill = group, y = group), alpha = 0.8) + scale_fill_manual(values = colours, name = "Taxonomic group") + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + xlim(0, 60000) 
ggsave("zc4688/honours/analyses/chordates/figures/lengths.ridges.zoomed.pdf",)
ggplot(plotdata) + geom_density(aes(x = `Unit length`, fill = group), alpha = 0.8) + scale_fill_manual(values = colours, name = "Taxonomic group") + theme_bw() + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)))
ggsave("zc4688/honours/analyses/chordates/figures/lengths.density.pdf",)

p <- ggtree(tree, size=0.5) + 
  geom_tiplab(size=1) 



p <- p + geom_facet(
  data = plotdata,
  mapping = aes(x = ifelse(`Unit length` < 60000, `Unit length`/1000, 60), fill = group),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
)  +
  geom_facet(
    data = plotdata,
    mapping = aes(x = log(CN, base = 10), fill = group),
    geom = geom_col,      
    panel = "log10(CN)",
    width = 0.6
  ) +
  geom_facet(
    data = plotdata,
    mapping = aes(x = `Total length`/1000000, fill = group),
    geom = geom_col,      
    panel = "Total length (Mb)",
    width = 0.6
  )  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 23), "Tree") + xlim_expand(c(0, 60), "Unit length (Kb)") +
  scale_fill_manual(values = colours) + labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(3, 1, 1, 1, 1))

#trees from timetree are saved at:
#ggsave("zc4688/honours/analyses/chordates/figures/tree_current.png", height = 20, width = 15)
#ggsave("zc4688/honours/analyses/chordates/figures/tree_current.pdf", height = 20, width = 15)


ggsave("zc4688/honours/analyses/chordates/figures/ncbitree.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/ncbitree.pdf", height = 20, width = 15)



########################### Read in and prepare structure data
structure_files <- list.files(directory, pattern = "\\.structure.tsv$", full.names = TRUE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each file and read the JSON into a data frame
for (file in structure_files) {
  data <- read.delim(file, header = T, sep="\t")
  samplename <- str_replace(basename(file), ".refmorph.structure.tsv", "")
  data <- data %>% 
    select(Type, HMMbegin, HMMend, Envstart, Envend, `X..HMM`, `Target.length`, Strand) %>% 
    mutate(`Sample ID` = samplename) 
  data_list[[file]] <- data
}

structure <- bind_rows(data_list)

# Get a list of all FASTA files in the directory
fasta_files <- list.files(directory, pattern = "\\.rDNA.refmorph.fasta$", full.names = TRUE)

# Initialize an empty list to store the sequence names from all files
all_sequence_names <- list()

# Loop through each FASTA file and extract sequence names
for (file in fasta_files) {
  samplename <- str_replace(basename(file), ".rDNA.refmorph.fasta", "")
  # Read the sequence names from the current file
  sequence_names <- grep("^>", readLines(file), value = TRUE)
  
  # Remove the '>' character from the sequence names
  sequence_names <- sub("^>", "", sequence_names)
  
  data <- data_frame(`Sample ID` = samplename, refmorph = sequence_names)
  
  # Add the sequence names to the list, with the filename as the key
  all_sequence_names[[file]] <- data
}

refmorphs <- bind_rows(all_sequence_names)

structure <- structure %>% left_join(refmorphs)  %>% 
  mutate(morphstart = as.numeric(str_replace_all(str_extract(refmorph, "(:\\d*)"), ":", "")),
         morphend = as.numeric(str_replace_all(str_extract(refmorph, "(-\\d*)"), "-", ""))) %>%  
  left_join(lengths) %>% 
  filter(!is.na(`Unit length`)) %>% 
  mutate(Start = ifelse(Strand == "+", Envstart - morphstart, morphend - Envend + 1),
         End = ifelse(Strand == "+", Envend - morphstart, morphend - Envstart + 1) )%>% 
  dplyr::rename("Hmmpalign" = "X..HMM") %>% 
  relocate(Species, .before = everything())%>% 
  mutate(Species = str_replace_all(Species, "_", " "))





################################################################################ Create fasta files for msa

msadf <- structure %>% 
  group_by(`Sample ID`, Type) %>% 
  filter((End - Start) == max(End - Start)) %>% 
  arrange(`Sample ID`, Type, Start) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(`Sample ID`, Type, Start, End, refmorph) %>% 
  mutate(length = End - Start, 
         name = paste(`Sample ID`, refmorph, sep = "_"))



eighteen <- msadf %>% 
  filter(Type == "18S_rRNA") %>% 
  #mutate(outside = length > (mean(length) + 3*sd(length)) | 
           #length < (mean(length) - 3*sd(length)))   %>% 
  #filter(outside == FALSE) %>% 
  select(name, Start, End)

write.table(eighteen, "zc4688/honours/analyses/chordates/new/msa/eighteen.csv", sep = ",", quote = F, row.names = F, col.names = F)
twoeight <- msadf %>% filter(Type == "28S_rRNA") %>% 
  select(name, Start, End)
write.table(twoeight, "zc4688/honours/analyses/chordates/new/msa/twoeight.csv", sep = ",", quote = F, row.names = F, col.names = F)


fiveeight <- msadf %>% 
  filter(Type == "5_8S_rRNA") %>% 
  select(name, Start, End)
write.table(fiveeight, "zc4688/honours/analyses/chordates/new/msa/fiveeight.csv", sep = ",", quote = F, row.names = F, col.names = F)


################################################################################ Structure plot


structureplot <- structure %>%
  filter(Species %in% common_organisms) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())

#structureplot <- structureplot %>% mutate(Hmmpalign = Hmmpalign * 100)







tmp <- structureplot %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + geom_facet(
  data = structureplot,
  mapping = aes(x = 0, xend = ifelse(`Unit length` < 35000, `Unit length`, 35000)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey80",
  alpha = 0.3
) + new_scale_color() + geom_facet(
  data = structureplot,
  mapping = aes(x = Start, xend = End, color = Type),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 6), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

p

ggsave("zc4688/honours/analyses/chordates/figures/structure_worepeats.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/structure_worepeats.pdf", height = 20, width = 15)


p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

test <- structureplot %>% select(label, `Unit length`, group) %>% dplyr::rename("testgroup" = "group") %>% group_by(label) %>% 
  slice_head(n=1) %>% ungroup()
p <- p  + new_scale_color() + geom_facet(
  data = structureplot,
  mapping = aes(x = Start, xend = End, color = Type, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = "Unit Structure"
)  + geom_facet(
  data = test,
  mapping = aes(x = ifelse(`Unit length`/1000 < 60, `Unit length`/1000, 60), fill = testgroup),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 6), "Tree")+labs(alpha = "% HMM Match") + scale_fill_manual(values = colours) + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1),
    fill = guide_none()) + coord_cartesian(clip = 'off')

facet_widths(p, widths = c(3, 2, 1))


ggsave("zc4688/honours/analyses/chordates/figures/structure_seplengths.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/structure_seplengths.pdf", height = 20, width = 15)


################################################################ Trees with repeats
repeatmasker_files <- list.files("/g/data/te53/zc4688/honours/analyses/chordates/repeatmasker", 
                                 pattern = "rDNA.refmorph.fasta.out.gff", 
                                 full.names = TRUE, 
                                 recursive = TRUE)
rm <- list()
for (i in 1:length(repeatmasker_files)) {
  x <- read.table(repeatmasker_files[i],header = F, sep = "\t")
  x <- x %>%
    mutate(`Sample ID` = str_replace(basename(repeatmasker_files[i]), ".rDNA.refmorph.fasta.out.gff", ""))
  rm[[i]] <- x
}

rm <- bind_rows(rm)

rm <- rm %>% 
  separate(V9, into = c("id", "motif", "class", "length"), sep = ";")  %>% 
  mutate(type = str_replace(motif, "Target Motif:", ""), class = str_replace(class, "RepeatType ClassFamily:", "")) %>% 
  select(`Sample ID`, V4, V5, V7, type, class) %>% dplyr::rename("start" = "V4", "end" = "V5", "strand" = "V7") %>% 
  left_join(classifications) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything()) %>% 
  filter(!(class == "rRNA" | str_detect(class, "Simple") | str_detect(class, "Satellite") | str_detect(class, "Low_complexity")))

rm <- rm %>%   
  filter(!(str_detect(class, "Simple|Satellite|RNA|complexity|Unknown")))%>% 
  mutate(repeatgroup = 
           case_when(str_detect(class, "LINE") ~ "LINE", 
                     str_detect(class, "SINE") ~ "SINE", 
                     str_detect(class, "Retro") ~ "Retrotransposon", 
                     class == "Simple_repeat" ~ "Simple", 
                     str_detect(class, "LTR") ~ "LTR", 
                     str_detect(class, "DNA|Helitron") ~ "DNA transposon", 
                     str_detect(class, "Satellite") ~ "Satellite", 
                     class == "Low_complexity" ~ "Low complexity", 
                     TRUE ~ "Other"))

trf_files <- list.files("/g/data/te53/zc4688/honours/analyses/chordates/trf", 
                        pattern = ".rDNA.refmorph.fasta.2.7.7.80.10.50.2000.dat", 
                        full.names = TRUE, 
                        recursive = TRUE)
trf <- list()
for (i in 1:length(trf_files)) {
  x <- read_table(trf_files[i],col_names  = F, skip=13)
  x <- x %>%
    mutate(`Sample ID` = str_replace(basename(trf_files[i]), ".rDNA.refmorph.fasta.2.7.7.80.10.50.2000.dat", ""))
  if (nrow(x) > 0) {
    trf[[length(trf) + 1]] <- x  # Only append if not empty
  }
}

trf <- bind_rows(trf)


####################### Summary plot

#trf gives overlapping alignments
trf_ranges <- GRanges(seqnames = trf$`Sample ID`, ranges = IRanges(start = trf$X1, end = trf$X2))
trf_ranges <- GenomicRanges::reduce(trf_ranges, min.gapwidth = 10)
trf_ranges <- data.frame(`Sample ID` = as.character(seqnames(trf_ranges)), start = start(trf_ranges), end = end(trf_ranges))
colnames(trf_ranges) <- c("Sample ID", "start", "end")

rm_ranges <- GRanges(seqnames = rm$`Sample ID`, ranges = IRanges(start = rm$start, end = rm$end))
rm_ranges <- GenomicRanges::reduce(rm_ranges, min.gapwidth = 10)
rm_ranges <- data.frame(`Sample ID` = as.character(seqnames(rm_ranges)), start = start(rm_ranges), end = end(rm_ranges))
colnames(rm_ranges) <- c("Sample ID", "start", "end")

#only IGS

trf_ranges <- structuremetadata %>% 
  left_join(trf_ranges) %>% 
  filter(start > `28S_rRNA_End`) %>% 
  group_by(`Sample ID`) %>% 
  summarise(`Simple repeats` = sum(end - start)/IGS_length, .groups = "drop") %>% 
  unique()

rm_ranges <- structuremetadata %>% 
  left_join(rm_ranges) %>% 
  filter(start > `28S_rRNA_End`) %>% 
  group_by(`Sample ID`) %>% 
  summarise(TEs = sum(end - start)/IGS_length, .groups = "drop") %>% 
  unique()

classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>% 
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>% 
  mutate(prop_repeat = TEs + `Simple repeats`) %>% 
  select(group, TEs, `Simple repeats`) %>% 
  filter(group %in% common_class) %>% 
  pivot_longer(cols = c(TEs, `Simple repeats`), names_to = "variable", values_to = "value") %>% view %>% 
  ggplot() + 
  geom_boxplot(aes(x = group, y = value, fill = variable)) + 
  theme_bw() + 
  labs(x = NULL, y = "Proportion of IGS", fill = NULL) + 
  scale_fill_brewer(palette = "Pastel1") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 



####################### Repeat tree

trf <- trf %>% left_join(classifications) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())

repeats <- bind_rows(rm %>% 
                       select(label, `Sample ID`, group, start, end, repeatgroup), 
                     trf %>% 
                       mutate(start = as.integer(X1), end = as.integer(X2)) %>% 
                       select(label, `Sample ID`, group, start, end) %>% 
                       mutate(repeatgroup = "Simple repeat")
                     ) %>% 
  left_join(structure %>% filter(Type == "28S_rRNA") %>% select(`Sample ID`, Start, End)) %>% 
  filter(start > End | end < Start)
#for trees with timetree, add xlim_expand(c(0, 40), "Tree")


repeat_colors <- c("DNA transposon" = "#7570B3", "LINE" = "#E7298A", "LTR" = "#D95F02", "Retrotransposon" = "#E6AB02", "SINE" = "#1B9E77", "Simple repeat" = "grey3")
test <- structureplot %>% select(label, group)

p <- ggtree(tree) %<+% test +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + geom_facet(
  data = structureplot,
  mapping = aes(x = 0, xend = ifelse(`Unit length` < 35000, `Unit length`, 35000)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey80",
  alpha = 0.3
) +
  new_scale_color() + geom_facet(
    data = structureplot,
    mapping = aes(x = Start, xend = End, color = Type),
    geom = geom_segment,       # Horizontal bar chart
    panel = " "
  )  +  
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  new_scale_color() +
  geom_facet(
    data = repeats %>% filter(end < 35000),
    mapping = aes(x = start, xend = end, color = repeatgroup),
    geom = geom_segment, 
    panel = " ")+ scale_color_manual(
      values = repeat_colors,
      name = "Repeat Type"
    ) +
  scale_y_discrete() +  
  theme_tree2() + scale_alpha(range = c(0.1, 1))+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

#for trees with timetree, add xlim_expand(c(0, 40), "Tree")

p


#structure on old timetree tree saved at:
#ggsave("zc4688/honours/analyses/chordates/figures/barrnaptree.png", height = 20, width = 15)
#ggsave("zc4688/honours/analyses/chordates/figures/barrnaptree.pdf", height = 20, width = 15)

ggsave("zc4688/honours/analyses/chordates/figures/structure_wrepeats.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/structure_wrepeats.pdf", height = 20, width = 15)


############################### Read in primary structure file for per group detailed plots

primarystructure <- read_table("zc4688/honours/analyses/chordates/new/structures.txt", col_names = F)


primarystructure <- primarystructure %>% select(X1, X4, X6, X7, X10, X11, X13) %>% 
  setNames(c("Sample ID", "type", "hmmstart", "hmmend", "envstart", "envend", "strand"))  %>% 
  left_join(refmorphs) %>% 
  mutate(morphstart = as.numeric(str_replace_all(str_extract(refmorph, "(:\\d*)"), ":", "")),
         morphend = as.numeric(str_replace_all(str_extract(refmorph, "(-\\d*)"), "-", "")))
  
primarystructure <- left_join(lengths %>% select(-type), primarystructure) %>% 
  mutate(Start = ifelse(strand == "+", envstart - morphstart, morphend - envstart + 1),
         End = ifelse(strand == "+", envend - morphstart, morphend - envend + 1) ) 



primarystructure <- primarystructure %>%
  mutate(Species = str_replace_all(Species, "_", " ")) %>% 
  left_join(classifications) %>% 
  filter(Species %in% common_organisms) 


primarystructure <- primarystructure %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())




tmp <- primarystructure %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + geom_facet(
  data = primarystructure,
  mapping = aes(x = 0, xend = ifelse(`Unit length` < 35000, `Unit length`, 35000)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey80",
  alpha = 0.3
) +
  new_scale_color() + geom_facet(
    data = primarystructure,
    mapping = aes(x = Start, xend = End, color = type),
    geom = geom_segment,       # Horizontal bar chart
    panel = " "
  )  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 6), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

p

ggsave("zc4688/honours/analyses/chordates/figures/structure_primary.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/structure_primary.pdf", height = 20, width = 15)





make_classplots <- function(plotclass){
  x <- primarystructure %>% filter(group == plotclass)
  height <- case_when(length(unique(x$label)) < 100 ~ 8, 
                      length(unique(x$label)) > 100 & length(unique(x$label)) < 200 ~ 12,
                      TRUE ~ 15)
  
  linewidth <- case_when(length(unique(x$label)) < 100 ~ 5, 
                      length(unique(x$label)) > 100 & length(unique(x$label)) < 200 ~ 2,
                      TRUE ~ 1)
  
  
  xtree <- ape::drop.tip(tree, setdiff(tree$tip.label, x$label))
  
  p <- ggtree(xtree)  +
    geom_tiplab(size = case_when(length(unique(x$label)) < 100 ~ 3, 
                                 length(unique(x$label)) > 100 & length(unique(x$label)) < 200 ~ 2,
                                 TRUE ~ 1.5))
  
  p <- p +  geom_facet(
    data = x,
    mapping = aes(x = 0, xend = ifelse(`Unit length` < 15000, `Unit length`, 15000)),
    geom = geom_segment,       # Horizontal bar chart
    panel = "Unit Structure",
    color = "grey80",
    alpha = 0.3, linewidth = linewidth
  ) +
    geom_facet(
    data = x %>% filter(End < 15000),
    mapping = aes(x = Start, xend = End, color = type),
    geom = geom_segment, linewidth = linewidth,
    panel = "Unit Structure"
  )  +
    theme(
      legend.key.size = unit(0.5, "in")) +
    scale_y_discrete() +  
    theme_tree2() +xlim_expand(c(0, 28), "Tree") + 
    scale_color_manual(values = c("18S_rRNA" = "#F75F86", 
                                  "28S_rRNA" = "cornflowerblue", 
                                  "5_8S_rRNA" = "#50C878"
                                  ), name = "rRNA")+
    guides(
      color = guide_legend(order = 1),
      fill = "none") + 
    coord_cartesian(clip = 'off') + scale_fill_manual(values=colours)
  
  facet_widths(p, widths = c(3, 3, 1))
  
  ggsave(glue::glue("zc4688/honours/analyses/chordates/figures/{plotclass}_norepeats.pdf"), height = height, width = 8, units = "in")
  ggsave(glue::glue("zc4688/honours/analyses/chordates/figures/{plotclass}_norepeats.png"), height = height, width = 8, units = "in")
  
  
  
  p <- ggtree(xtree)  +
    geom_tiplab(size = case_when(length(unique(x$label)) < 100 ~ 3, 
                                 length(unique(x$label)) > 100 & length(unique(x$label)) < 200 ~ 2,
                                 TRUE ~ 1.5))
  
  p <- p + geom_facet(
    data = x,
    mapping = aes(x = 0, xend = `Unit length` ),
    geom = geom_segment,       # Horizontal bar chart
    panel = "Unit Structure",
    color = "grey80",
    alpha = 0.3, linewidth = linewidth
  )+ geom_facet(
    data = x,
    mapping = aes(x = Start, xend = End, color = type),
    geom = geom_segment, linewidth = linewidth, 
    panel = "Unit Structure"
  )  +  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "rRNA")+
    theme(
      legend.key.size = unit(0.5, "cm")) +
    new_scale_color() +
    geom_facet(
      data = repeats %>% filter(group == plotclass),
      mapping = aes(x = start, xend = end, color = repeatgroup),
      geom = geom_segment, linewidth = linewidth, 
      panel = "Unit Structure") +
    scale_color_manual(values = repeat_colors, name = "Repeat type") +
    scale_y_discrete() +  
    theme_tree2() +xlim_expand(c(0, 28), "Tree") + 
     coord_cartesian(clip = 'off') + 
    scale_fill_manual(values=colours) + guides(fill = "none") 
    
  
  facet_widths(p, widths = c(3, 3, 1))
  
  ggsave(glue::glue("zc4688/honours/analyses/chordates/figures/{plotclass}tree.pdf"), height = height, width = 8, units = "in")
  ggsave(glue::glue("zc4688/honours/analyses/chordates/figures/{plotclass}tree.png"), height = height, width = 8, units = "in")
  
}
  
make_classplots("Mammalia")
make_classplots("Aves")
make_classplots("Amphibia")
make_classplots("Actinopterygii")
make_classplots("Chondrichthyes")
make_classplots("Lepidosauria")
make_classplots("Testudines and Crocodylia")
make_classplots("Tunicata")


#ind plots
ggplot() + geom_segment(data = structure %>% filter(Species == "Eospalax fontanierii"), aes(x = 0, xend = `Unit length`, y = 0), linewidth = 5, color = "grey")+ 
  geom_segment(data = structure %>% filter(Species == "Eospalax fontanierii"), aes(x = Start, xend = End, y = 0, color = Type), linewidth = 5)+  
  geom_segment(data = structure %>% filter(Species == "Homo sapiens"), aes(x = 0, xend = `Unit length`, y = 0.5), linewidth = 5, color = "grey")+ 
  geom_segment(data = structure %>% filter(Species == "Homo sapiens"), aes(x = Start, xend = End, y = 0.5, color = Type), linewidth = 5)+ 
  geom_segment(data = structure %>% filter(Species == "Antrozous pallidus"), aes(x = 0, xend = `Unit length`, y = 1), linewidth = 5, color = "grey")+ 
  geom_segment(data = structure %>% filter(Species == "Antrozous pallidus"), aes(x = Start, xend = End, y = 1, color = Type), linewidth = 5)+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "rRNA") + 
  new_scale_color() + 
  geom_segment(data = repeats %>% filter(label == "Eospalax fontanierii"), aes(x = start, xend = end, y = 0, color = repeatgroup), linewidth = 5) + 
  geom_segment(data = repeats %>% filter(label == "Homo sapiens"), aes(x = start, xend = end, y = 0.5, color = repeatgroup), linewidth = 5) +
  geom_segment(data = repeats %>% filter(label == "Antrozous pallidus"), aes(x = start, xend = end, y = 1, color = repeatgroup), linewidth = 5) +
  scale_color_manual(
    values = repeat_colors,
    name = "Repeat Type"
  )   + coord_fixed(ratio = 5000, clip = "off")+
  theme_void() 
############################### Boxplots

boxplots <- structure %>% 
  group_by(`Sample ID`, Type) %>% 
  filter((End - Start) == max(End - Start)) %>% 
  arrange(`Sample ID`, Type, Start) %>% 
  slice_head(n = 1) %>%  #if there are multiple of same length, take the first
  ungroup() %>% 
  select(`Sample ID`, Type, Start, End, group, `Unit length`) %>% pivot_wider(
    names_from = Type,
    values_from = c(Start, End),
    names_glue = "{Type}_{.value}"
  ) %>% 
  mutate(ITS1 = `5_8S_rRNA_Start` - `18S_rRNA_End`,
         `ITS2` = `28S_rRNA_Start` - `5_8S_rRNA_End`,
         `18S_length` = `18S_rRNA_End` - `18S_rRNA_Start`,
         `28S_length` = `28S_rRNA_End` - `28S_rRNA_Start`,
         `IGS_length` = `Unit length` - `28S_rRNA_End`)

#5 species are missing ITS measuremnts due to lack of 5.8S. 4 of these don't have 5.8s alignments in refmorph. the remaining (GCA_963675345, tunciate) seems to genuinely be missing part of the 5.8S unit - consistently only cover 75% of the model so don't pass filtering

its1 <- boxplots %>% 
  left_join(refmorphs) %>% 
  mutate(name = paste(`Sample ID`, refmorph, sep = "_")) %>% 
  select(name, `18S_rRNA_End`, `5_8S_rRNA_Start`) %>% 
  na.omit()
write.table(its1, "zc4688/honours/analyses/chordates/motif/its1.csv", sep = ",", quote = F, row.names = F, col.names = F)

its2 <- boxplots %>% 
  left_join(refmorphs) %>% 
  mutate(name = paste(`Sample ID`, refmorph, sep = "_")) %>% 
  select(name, `5_8S_rRNA_End`, `28S_rRNA_Start`) %>% 
  na.omit()
write.table(its2, "zc4688/honours/analyses/chordates/motif/its2.csv", sep = ",", quote = F, row.names = F, col.names = F)

structuremetadata <- boxplots

boxplots <- boxplots %>% 
  group_by(group) %>% mutate(n = n()) %>% ungroup() %>%  filter(n >=5) %>% select(-n) %>% 
  pivot_longer(cols = c(`ITS1`, ITS2, `18S_length`, `28S_length`, `IGS_length`, `Unit length`), names_to = "Variable", values_to = "Value") %>% 
  dplyr::rename(label = group) %>% 
  select(label, Variable, Value)



ggplot(boxplots, aes(y = Value, fill = label)) +
  geom_boxplot(coef = 1.5, outliers = FALSE, notch = TRUE) +
  facet_wrap(~ Variable, scales = "free_y")  + theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = colours) +
  labs(fill = "Taxonomic Group")

ggsave("zc4688/honours/analyses/chordates/figures/boxplot.png", height = 8, width = 12)


############################### Boxplots with tree
classtree <- read.tree("/g/data/te53/zc4688/honours/trees/classes_ncbi.phy")


classtree$tip.label[4] <- "Testudines and Crocodylia"
common_class <- intersect(boxplots$label, classtree$tip.label)


boxplots <- boxplots %>%
  filter(label %in% common_class)



classtree$edge.length <- NULL
classtree <- ape::drop.tip(classtree, setdiff(classtree$tip.label, common_class))

############################### ITS only
p <- ggtree(classtree) + 
  geom_tiplab(size=3) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours) + theme(legend.position = "none")



p <- p + geom_facet(
  data = boxplots %>% filter(Variable == "ITS1"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS1 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = boxplots %>% filter(Variable == "ITS2"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS2 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position = "none") + labs(color = NULL, fill = NULL) + scale_fill_brewer(palette = "Set1") 

facet_widths(p, widths = c(1, 1, 1))
ggsave("zc4688/honours/analyses/chordates/figures/itsboxplot.pdf", height = 8, width = 12)


############################### All lengths

p <- ggtree(classtree) + 
  geom_tiplab(size=3, offset = 0.5) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours) + theme(legend.position = "none")


p <- p + new_scale_color() + geom_facet(
  data = boxplots %>% filter(Variable == "Unit length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = label),
  geom = geom_violin,       
  panel = "Unit Length (Kb)",
  notch = TRUE, outliers = F, coef = 2.5
) + geom_facet(
  data = boxplots %>% filter(Variable == "ITS1"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = label),
  geom = geom_boxplot,       
  panel = "ITS1 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = boxplots %>% filter(Variable == "ITS2"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = label),
  geom = geom_boxplot,       
  panel = "ITS2 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +  geom_facet(
  data = boxplots %>% filter(Variable == "IGS_length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = label),
  geom = geom_boxplot,       
  panel = "IGS Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +
  geom_facet(
    data = boxplots %>% filter(Variable == "18S_length"),
    mapping = aes(x = Value/1000, fill = label, group = interaction(label, Variable)),
    geom = geom_boxplot,       
    panel = "18S Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
    data = boxplots %>% filter(Variable == "28S_length"),
    mapping = aes(x = Value/1000, fill = label, group = interaction(label, Variable)),
    geom = geom_boxplot,       
    panel = "28S Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position =  "none")  + scale_fill_manual(values = colours) +
  coord_cartesian(clip = "off")

facet_widths(p, widths = c(1, 1, 1))
ggsave("zc4688/honours/analyses/chordates/figures/combinedboxplots.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/combinedboxplots.png", height = 8, width = 12)


############################### Summary stats
boxplots %>% 
  group_by(label, Variable) %>% 
  summarise(Median = median(Value, na.rm = TRUE), Min = min(Value, na.rm = TRUE),  Max = max(Value, na.rm = TRUE), Species = n()) %>% 
  view


############################### Metadata table

metadata <- metadata %>% 
  left_join(structuremetadata)

metadata <- metadata %>% 
  select(-"18S_rRNA_Start", -"28S_rRNA_Start", -"5_8S_rRNA_Start", -"18S_rRNA_End", -"28S_rRNA_End", -"5_8S_rRNA_End") %>% 
  dplyr::rename("Number of morphs" = "rDNA_details.Number of morphs")

write.table(metadata, "zc4688/honours/metadata/metadata.tsv", sep = "\t", quote = F, row.names = F)

################################correlations
metadata <- read.delim("zc4688/honours/metadata/metadata.tsv", sep = "\t")
moremetadata <- read.delim("zc4688/honours/metadata/assembly_metadata.tsv")
metadata <- left_join(metadata, moremetadata %>% separate(accession, into = c("Sample.ID", "number"), sep = "\\."))
summary(lm(`Unit.length` ~ X28S_length, data = metadata))

ggplot(metadata, aes(x = X28S_length, y = Unit.length, color = group)) + 
  geom_point() + 
  scale_color_manual(values = colours) + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 120000))

################################Conservation

eighteen <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/fiveeight_clusters/cleaned.alignment.txt")

shannon_entropy <- function(column) {
  column <- column[column != "-"]
  number <- length(column)
  freqs <- table(column) / length(column)  # Frequency of each base
  -sum(freqs * log2(freqs), na.rm = TRUE)  # Shannon entropy formula
}

# Convert alignment to matrix
msa_matrix <- as.matrix(eighteen)

# Apply function to each column
entropy_values <- apply(msa_matrix, 2, shannon_entropy)

counts <- sapply(as.data.frame(msa_matrix), function(col) sum(col != "-"))


tmp <- data.frame(Position = 1:length(entropy_values), Score = entropy_values, count = counts)



cp_result <- cpt.mean(entropy_values, method = "PELT", penalty = "MBIC")
plot(cp_result, main = "Change Points in Shannon Entropy")
ggplot()+  
  geom_area(data = tmp, aes(x = Position, y = Score), linewidth = 0.4, fill = "cornflowerblue", alpha = 0.4) + 
  geom_point(data = tmp %>% filter(count > 50), aes(x = Position, y = Score, color = Score), size = 0.6) + 
  scale_color_gradient2(low = "deeppink",  mid = "#F4A582", high = "limegreen", na.value = "grey", midpoint = 1) + theme_bw() +
  labs(title = "Shannon entropy for 5.8S rRNA in chordates", y = "Shannon entropy")

ggplot()+  
  geom_area(data = tmp, aes(x = Position, y = Score), linewidth = 0.4, fill = "cornflowerblue", alpha = 0.4) + 
  geom_point(data = tmp %>% filter(count > 50), aes(x = Position, y = Score, color = Score), size = 0.6) + 
  scale_color_gradient(high = "deeppink", low = "limegreen", na.value = "grey") + theme_bw() +
  labs(title = "Shannon entropy for 5.8S rRNA in chordates", y = "Shannon entropy")


############################### 18S tree
eighteentree <- read.tree("/g/data/te53/zc4688/honours/analyses/chordates/new/msa/eighteen_clusters/eighteen.treefile")

treenames <- as.data.frame(eighteentree$tip.label)


treenames <- treenames %>% 
  mutate(`Sample ID` = str_extract(`eighteentree$tip.label`, "^GCA_\\d+|^GCF_\\d+|hg002"))

treenames <- treenames %>% left_join(classifications %>% select(`Sample ID`, Species))

eighteentree$tip.label <- as.character(treenames$label)
common_organisms <- intersect(plotdata$label, eighteentree$tip.label)

eighteentree <- ape::drop.tip(eighteentree, setdiff(eighteentree$tip.label, common_organisms))
eighteentree$edge.length <- NULL
#eighteentree$tip.label <- str_replace_all(tree$tip.label, "_", " ")


tmp <- plotdata %>% select(label, group)

ggtree(eighteentree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") +
  coord_cartesian(clip = 'off') 

################################


common_cotree <- intersect(eighteentree$tip.label, tree$tip.label)

eighteentree <- ape::drop.tip(eighteentree, setdiff(eighteentree$tip.label, common_cotree))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_cotree))
A <- cbind(tree$tip.label, tree$tip.label)
x <- A %>% as.data.frame(.)  %>% dplyr::rename("Species" = "V1") %>% left_join(classifications)
x$color <- colours[x$group]

tree <- ladderize(tree)
eighteentree <- root(eighteentree, outgroup = c("Branchiostoma lanceolatum", "B. floridae x B. belcheri", "B. floridae x B. japonicum"), resolve.root = T)

eighteentree <- ladderize(eighteentree)

dist.topo(tree, eighteentree)
comparePhylo(tree, eighteentree, plot = TRUE)

cophyloplot(tree, eighteentree, assoc = A, show.tip.label = F, space=500, col = x$color) 


co <- eighteentree$tip.label[order(match(eighteentree$tip.label, tree$tip.label))]
newtr2 <- rotateConstr(eighteentree, co)



############################################ Motif results

tmp <-  read.table("zc4688/honours/analyses/chordates/motif/its1_remasked/meme.txt", sep = ",") %>% 
      dplyr::rename(raw = V1)


#need to improve parsing (some positions etc not working). perhaps include length. create plot showing all significant motifs in diff colours

#include analysis of no. of each motif
df.temp <- tmp %>%
    dplyr::mutate(
        row.num = rownames(.),
        nchar.1 = nchar(raw),
        nchar.2 = 16,
        start = nchar.1 - nchar.2,
        end = nchar.1
      ) %>%
    dplyr::mutate(
        five = stringr::str_sub(raw, start, end),
        end.sym = stringr::str_sub(raw, 1, 2)
      ) %>%
    dplyr::filter(five == " in BLOCKS format" | end.sym == "//") %>%
    dplyr::mutate(row.num = ifelse(raw == "//", as.numeric(row.num) - 1, as.numeric(row.num) + 3)) %>%
    dplyr::select(1, 2)
    


df.motif <- NULL

for (i in seq(1, nrow(df.temp), 2)) {
  start.row <- df.temp$row.num[i]
  
  end.row <- df.temp$row.num[i + 1]
  
  df.motif.temp <- tmp[start.row:end.row, ] %>%
    as.data.frame()
  
  colnames(df.motif.temp) <- "raw"
  
  df.motif.temp <- df.motif.temp %>%
    dplyr::mutate(
      motif.num = paste0("Motif.", (i + 1) / 2),
      input.seq.name = "",
      input.seq.motif = "",
      input.seq.pos = ""
    )
  
  for (j in 1:nrow(df.motif.temp)) {
    df.motif.temp$input.seq.name[j] <- stringr::str_split(df.motif.temp$raw[j], " ")[[1]][1]
    df.motif.temp$input.seq.motif[j] <- stringr::str_split(df.motif.temp$raw[j], " ")[[1]][(length(stringr::str_split(df.motif.temp$raw[j], " ")[[1]]) - 3)]
    df.motif.temp$input.seq.pos[j] <- as.numeric(str_extract(df.motif.temp$raw[j], "(( \\d+))"))
    
  }
  
  df.motif <- rbind(df.motif, df.motif.temp)
}

#look at all major motifs
motif1 <- df.motif %>% filter(motif.num %in% c("Motif.1", "Motif.2", "Motif.3", "Motif.4", "Motif.5")) %>% 
  mutate(`Sample ID` = str_extract(input.seq.name, "^GCA_\\d+")) %>% 
  left_join(classifications) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything()) %>% 
  filter(label %in% common_organisms) %>% 
  left_join(boxplots %>% select(`Sample ID`, ITS1))


#look at no.1 motif
motif1 <- df.motif %>% filter(motif.num %in% c("Motif.3")) %>% 
  mutate(`Sample ID` = str_extract(input.seq.name, "^GCA_\\d+|^GCF_\\d+")) %>% 
  left_join(classifications) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything()) %>% 
  filter(label %in% common_organisms) %>% 
  left_join(boxplots %>% select(`Sample ID`, ITS1, ITS2)) %>% 
  group_by(`Sample ID`) %>% 
  mutate(number = n()) %>%
  ungroup()


#all motifs and no.
motif1 <- df.motif  %>% 
  filter(motif.num %in% c("Motif.1", "Motif.2", "Motif.3", "Motif.4", "Motif.5", "Motif.6", "Motif.7", "Motif.8")) %>% 
  #filter(motif.num %in% c("Motif.9", "Motif.10", "Motif.11", "Motif.12", "Motif.13", "Motif.14", "Motif.15", "Motif.16")) %>% 
  mutate(`Sample ID` = str_extract(input.seq.name, "^GCA_\\d+|^GCF_\\d+")) %>% 
  left_join(classifications) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything()) %>% 
  filter(label %in% common_organisms) %>% 
  left_join(boxplots %>% select(`Sample ID`, ITS1, ITS2)) %>% 
  group_by(`Sample ID`, motif.num) %>% 
  mutate(number = n()) %>%
  ungroup()
  
p <- ggtree(tree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") 
  
#for first 2  
its1 <- motif1 %>% select(label, ITS2, number, group) %>% unique
p <- p + new_scale_color() + 
  geom_facet(
    data = its1,
    mapping = aes(x = ITS2),
    geom = geom_col,
    color = "grey",
    panel = "Motif",
    size = 2 
  ) + geom_facet(
    data = motif1,
    mapping = aes(x = as.numeric(input.seq.pos), xend = as.numeric(input.seq.pos) + 50, color = motif.num),
    geom = geom_segment,       # Horizontal bar chart
    panel = "Motif",
    size = 1 
  )     + 
  geom_facet(
    data = its1,
    mapping = aes(x = number),
    geom = geom_col,
    color = "pink",
    panel = "number",
    size = 2 
  )  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 40), "Tree") +
  labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

p



#for last
its1 <- motif1 %>% select(label, ITS2, motif.num, number, group) %>% unique

p <- p + new_scale_color() + 
  geom_facet(
    data = its1 %>% select(label, ITS2) %>% unique,
    mapping = aes(x = ITS2),
    geom = geom_col,
    color = "grey",
    panel = "Motif",
    size = 2 
  ) + geom_facet(
    data = motif1,
    mapping = aes(x = as.numeric(input.seq.pos), xend = as.numeric(input.seq.pos) + 50, color = motif.num),
    geom = geom_segment,       # Horizontal bar chart
    panel = "Motif",
    size = 1 
  )     + 
  geom_facet(
    data = its1,
    mapping = aes(x = number, color = motif.num),
    geom = geom_col,
    panel = "Number of motifs",
    size = 2 
  )  +
    scale_y_discrete() +  scale_color_brewer(palette = "Dark2", name = "Motif name") +
    theme_tree2() + xlim_expand(c(0, 40), "Tree") +
    labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 
  
p

  
  
ggsave("zc4688/honours/analyses/chordates/figures/its2_motifs1-8.pdf", height = 20, width = 15)



###########################Related genes
server <- "https://rest.ensembl.org"
ext <- "/homology/symbol/human/"
type <- "?type=orthologues"
gene_ids <- c("FBL", "NOB1")

data_list <- list()

for (gene in gene_ids){
  r <- GET(paste(server, ext, gene, type, sep = ""), content_type("application/json"))

  orthologues <- content(r)
  tmp <- orthologues$data[[1]]$homologies

  tmp <- do.call(rbind, lapply(tmp, function(homology) {
    data.frame(
      species = homology$species,
      protein_id = homology$protein_id,
      id = homology$id,
      type = homology$type
    )
  }))
  
  tmp <- tmp %>% mutate(humangene = gene)
  
  data_list[[gene]] <- tmp
  
}

tmp <- bind_rows(data_list)

tmp$species <- str_replace_all(tmp$species, "_", " ")

common_organisms <- intersect(tolower(tmp$species), tolower(plotdata$label))


tmp <- read.delim("zc4688/honours/analyses/chordates/motif/its1.self.blastn.out", header = F)

mat <- tmp %>% 
  mutate(`Sample ID` = str_extract(V1, "^GCA_\\d+|^GCF_\\d+")) %>% 
  left_join(boxplots %>% select(`Sample ID`, ITS1)) %>% 
  mutate(prop = V4/ITS1) %>% 
  dplyr::select(V1, V2, prop) %>% 
  group_by(V1, V2) %>% 
  summarise(prop = max(prop), .groups = "drop") %>% 
  pivot_wider(names_from = V2, values_from = prop, values_fill = 0) %>% 
  column_to_rownames("V1") 

mat <- tmp %>% 
  mutate(`Sample ID` = str_extract(V1, "^GCA_\\d+|^GCF_\\d+")) %>% 
  left_join(boxplots %>% select(`Sample ID`, `18S_length`)) %>% 
  mutate(prop = V4/`18S_length`) %>% 
  dplyr::select(V1, V2, prop) %>% 
  group_by(V1, V2) %>% 
  summarise(prop = max(prop), .groups = "drop") %>% 
  pivot_wider(names_from = V2, values_from = prop, values_fill = 0) %>% 
  column_to_rownames("V1") 


mat <- mat[intersect(rownames(mat), colnames(mat)), intersect(rownames(mat), colnames(mat))]

