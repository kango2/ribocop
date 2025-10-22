.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
library(ggtree)
library(tidyverse)
library(phangorn)
library(ape)
library(jsonlite)
library(GenomicRanges)
library(taxize)
library(RColorBrewer)
library(ggnewscale)
library(Biostrings)
library(ggridges)
library(rstatix)
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

############################### Create 'lengths' df - metadata df containing lengths + taxonomy info for all successful species with duplicates removed

lengths <- combined_df %>% 
  filter(`Unit length` < 120000) %>% 
  filter(Species != "Antrozous_pallidus" & Species != "Leptobrachium_leishanense" & Species != "Choloepus_didactylus") %>% 
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
#  dplyr::select(class, order, species) %>% 
  dplyr::rename("Species" = "species")  

classifications <- classifications %>% 
  group_by(class, order, Species) %>% summarize(across(everything(), ~ last(.))) #Take refseq if both refseq and genbank are available

lengths[409,1] <- "Mus_musculus"

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
                   "Rangifer tarandus platyrhynchus", "Giraffa camelopardalis rothschildi", "Mus musculus") ~ "Mammalia",
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



groups <- unique(classifications$group)
colours <- setNames(brewer.pal(n = 11, name = "Paired"), groups)




########################################### Read in and prepare species tree

tree <- read.tree("/g/data/te53/zc4688/honours/trees/species_ncbi.phy")

nodes <- (tree$node.label)

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
  dplyr::select(Species, `Unit length`, `CN`, `Total length`,  `Sample ID`, Contigs, group, `Minimum unit length`, `Maximum unit length`, Q1, Q3)%>% 
  mutate(Q1 = as.numeric(Q1), Q3 = as.numeric(Q3)) %>% 
  dplyr::rename("label" = "Species") %>%  
  relocate(label, .before = everything())  %>% 
  filter(label %in% common_organisms)

common_class <- lengths %>% group_by(group) %>% summarise(n = n()) %>% filter(n > 5) %>% select(group) %>% unlist()

ggplot(plotdata %>% filter(group %in% common_class)) + 
  geom_density_ridges(aes(x = `Unit length`, fill = group, y = group), alpha = 0.8) + 
  scale_fill_manual(values = colours, name = "Taxonomic group") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
  labs(x = "Unit length (bp)")
ggsave("zc4688/honours/analyses/chordates/figures/lengths.ridges.pdf",)
ggplot(plotdata %>% filter(group %in% common_class)) + 
  geom_density_ridges(aes(x = `Unit length`, fill = group, y = group), alpha = 0.8) + 
  scale_fill_manual(values = colours, name = "Taxonomic group") + theme_bw() + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  xlim(0, 60000) +
  labs(x = "Unit length (bp)")
ggsave("zc4688/honours/analyses/chordates/figures/lengths.ridges.zoomed.pdf",)
ggplot(plotdata) + 
  geom_density(aes(x = `Unit length`, fill = group), alpha = 0.8) + 
  scale_fill_manual(values = colours, name = "Taxonomic group") + 
  theme_bw() + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + 
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)))+
  labs(x = "Unit length (bp)")
ggsave("zc4688/honours/analyses/chordates/figures/lengths.density.pdf",)

pairwise_res <- lengths %>% filter(group %in% common_class) %>% 
  pairwise_wilcox_test(`Unit length` ~ group, p.adjust.method = "BH")


pairwise_res <- pairwise_res %>% 
  select(group1, group2, p.adj) %>% 
  pivot_wider(names_from = group2, values_from = p.adj, values_fill = NA) %>% 
  column_to_rownames("group1")
pheatmap(pairwise_res, 
         cluster_rows = F, cluster_cols = F, 
         breaks = c(0, 0.005, 0.05, 1), color = c("firebrick", "darkorange", "grey"), 
         border_color = "white", na_col = "white", legend=F)


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

facet_widths(p, widths = c(2, 1, 1, 1, 1))



ggsave("zc4688/honours/analyses/chordates/figures/ncbitree.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/ncbitree.pdf", height = 20, width = 15)


p <- ggtree(tree, size=0.5) 


p <- p + geom_facet(
  data = plotdata,
  mapping = aes(x = ifelse(`Unit length` < 60000, `Unit length`/1000, 60), fill = group),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
)  +
  geom_facet(
    data = plotdata,
    mapping = aes(x = ifelse(`Unit length` >= 60000, 60, NA)),
    geom = geom_point,
    panel = "Unit length (Kb)",
    size = 1
  )+ 
  geom_facet(
    data = plotdata,
    mapping = aes(x = ifelse(`Q1` < 60000, `Q1`/1000, 60), xend = ifelse(`Q3` < 60000, `Q3`/1000, 60)),
    geom = geom_segment,      
    panel = "Unit length (Kb)",
    linewidth = 0.2
  )  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 15), "Tree") + xlim_expand(c(0, 60), "Unit length (Kb)") +
  scale_fill_manual(values = colours) + labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(1, 3))

ggsave("zc4688/honours/analyses/chordates/figures/lengths_quartiles.pdf", height = 20, width = 15)

p <- ggtree(tree, size=0.5) 



p <- p + geom_facet(
  data = plotdata,
  mapping = aes(x = ifelse(`Unit length` < 60000, `Unit length`/1000, 60), fill = group),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
)  +
  geom_facet(
    data = plotdata,
    mapping = aes(x = ifelse(`Unit length` >= 60000, 60, NA), color = group),
    geom = geom_point,
    panel = "Unit length (Kb)",
    size = 1
  )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 15), "Tree") + xlim_expand(c(0, 60), "Unit length (Kb)") +
  scale_fill_manual(values = colours) + scale_color_manual(values = colours) +labs(fill = "Taxonomic Group", color = "Taxonomic Group") +coord_cartesian(clip = 'off') 
facet_widths(p, widths = c(1, 3))

ggsave("zc4688/honours/analyses/chordates/figures/ncbi_lengthsonly_nolabs.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/ncbi_lengthsonly_nolabs.png", height = 20, width = 15)

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
    mapping = aes(x = ifelse(`Unit length` >= 60000, 60, NA), color = group),
    geom = geom_point,
    panel = "Unit length (Kb)",
    size = 1
  )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 25), "Tree") + xlim_expand(c(0, 60), "Unit length (Kb)") +
  scale_fill_manual(values = colours) + scale_color_manual(values = colours) +labs(fill = "Taxonomic Group", color = "Taxonomic Group") +coord_cartesian(clip = 'off') 
facet_widths(p, widths = c(1, 3))

ggsave("zc4688/honours/analyses/chordates/figures/ncbi_lengthsonly_withlabs.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/ncbi_lengthsonly_withlabs.png", height = 20, width = 15)

########################### Read in and prepare structure data
structure_files <- list.files(directory, pattern = "\\.structure.tsv$", full.names = TRUE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each file and read the JSON into a data frame
for (file in structure_files) {
  data <- read.delim(file, header = T, sep="\t", comment.char = "#")
  samplename <- str_replace(basename(file), ".refmorph.structure.tsv", "")
  data <- data %>% 
    select(Type, HMMbegin, HMMend, Envstart, Envend, `X..HMM`, `Target.length`, Strand) %>% 
    mutate(`Sample ID` = samplename) 
  data_list[[file]] <- data
}

structure <- bind_rows(data_list)

# Get list of refmorph fasta files
fasta_files <- list.files(directory, pattern = "\\.rDNA.refmorph.fasta$", full.names = TRUE)

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

tmp <- structureplot %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + geom_facet(
  data = structureplot,
  mapping = aes(x = 0, xend = ifelse(`Unit length` < 60000, `Unit length`, 60000)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey40",
  alpha = 0.3
) + new_scale_color() +   geom_facet(
  data = plotdata,
  mapping = aes(x = ifelse(`Unit length` >= 60000, 60000, NA)),
  geom = geom_point,
  color = "grey40",
  panel = " ",
  size = 1
)+ geom_facet(
  data = structureplot,
  mapping = aes(x = Start, xend = End, color = Type),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 23), "Tree")+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "rRNA")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')
facet_widths(p, widths = c(1, 3))


ggsave("zc4688/honours/analyses/chordates/figures/structure_worepeats.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/structure_worepeats.pdf", height = 20, width = 15)


p <- ggtree(tree) %<+% tmp +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

test <- structureplot %>% select(label, `Unit length`, group) %>% dplyr::rename("testgroup" = "group") %>% group_by(label) %>% 
  slice_head(n=1) %>% ungroup()
p <- p  +   geom_facet(
  data = test,
  mapping = aes(x = ifelse(`Unit length`/1000 < 60, `Unit length`/1000, 60), fill = testgroup),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
)  + geom_facet(
  data = test,
  mapping = aes(x = ifelse(`Unit length` >= 60000, 60, NA), color = testgroup),
  geom = geom_point,
  panel = "Unit length (Kb)",
  size = 1
)+ scale_fill_manual(values = colours, name = "Taxonomic Group") + 
  scale_color_manual(values = colours, name = "Taxonomic Group") +
      new_scale_color() +
  geom_facet(
  data = structureplot,
  mapping = aes(x = Start, xend = End, color = Type),
  geom = geom_segment,       # Horizontal bar chart
  panel = "Unit Structure") +
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 23), "Tree")+labs(alpha = "% HMM Match") +  scale_alpha(range = c(0.1, 1))+
  guides(
    color = guide_legend(order = 1),
    fill = guide_none()) + coord_cartesian(clip = 'off')

facet_widths(p, widths = c(1, 2, 3))


ggsave("zc4688/honours/analyses/chordates/figures/structure_seplengths.png", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/structure_seplengths.pdf", height = 20, width = 15)

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
  geom = geom_boxplot,       
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
  view %>% 
  write.table("zc4688/honours/metadata/summary_stats.tsv", sep = "\t", quote = F, row.names = F)



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

rm <- bind_rows(rm) %>% 
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


####################### Summary plot - igs

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
  summarise(`Simple repeats` = sum(end - start), .groups = "drop") %>% 
  unique()

rm_ranges <- structuremetadata %>% 
  left_join(rm_ranges) %>% 
  filter(start > `28S_rRNA_End`) %>% 
  group_by(`Sample ID`) %>% 
  summarise(TEs = sum(end - start), .groups = "drop") %>% 
  unique()

classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(total = `Simple repeats` + TEs, remain = IGS_length - TEs) %>% 
  filter(group %in% common_class)  %>% 
  ggplot(aes(x = remain, fill = group)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = colours) + theme_bw() + 
  labs(x = "IGS length (without TEs)", y = "Density", fill = "Taxonomic group")

ggsave("zc4688/honours/analyses/chordates/figures/igs.wo.tes.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/igs.wo.tes.png", height = 8, width = 12)

classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(total = `Simple repeats` + TEs, remain = IGS_length - TEs) %>% 
  filter(group %in% common_class) %>% pairwise_wilcox_test(remain ~ group, p.adj.method = "BH")


classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(`Length of repeats in IGS` = `Simple repeats` + TEs, `Length of IGS without TEs` = IGS_length - TEs, `Length of IGS without repeats` = IGS_length - `Length of repeats in IGS`) %>% 
  filter(group %in% common_class)  %>% 
  select(`Length of repeats in IGS`, `Length of IGS without TEs`, `Length of IGS without repeats`, IGS_length, group, `Sample ID`) %>% 
  pivot_longer(-c(group, `Sample ID`)) %>% 
  ggplot(aes(y = value, fill = group)) + 
  geom_boxplot(coef = 1.5, outliers = FALSE) + 
  labs(y = "Length (bp)") +
  scale_fill_manual(values = colours, name = "Taxonomic group") +
  facet_wrap(~name) + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_blank())
ggsave("zc4688/honours/analyses/chordates/figures/igs.repeat.summary.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/igs.repeat.summary.png", height = 8, width = 12)


classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(total = `Simple repeats` + TEs, remain = IGS_length - total) %>% 
  filter(group %in% common_class)  %>% 
  ggplot(aes(x = remain, fill = group)) + 
  geom_density(alpha = 0.5) + 
  scale_fill_manual(values = colours) + theme_bw() + 
  labs(x = "IGS length (without repeats)", y = "Density", fill = "Taxonomic group")

ggsave("zc4688/honours/analyses/chordates/figures/igs.wo.repeats.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/igs.wo.repeats.png", height = 8, width = 12)

classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(total = `Simple repeats` + TEs, remain = IGS_length - total) %>% 
  filter(group %in% common_class) %>% pairwise_wilcox_test(remain ~ group, p.adj.method = "BH")

classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>% 
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>% 
  mutate(prop_repeat = TEs + `Simple repeats`) %>% 
  select(group, TEs, `Simple repeats`) %>% 
  filter(group %in% common_class) %>% 
  pivot_longer(cols = c(TEs, `Simple repeats`), names_to = "variable", values_to = "value") %>% 
  ggplot() + 
  geom_boxplot(aes(x = group, y = value, fill = variable), notch = F) + 
  theme_bw() + 
  labs(x = NULL, y = "Length of repeats in IGS (bp)", fill = NULL) + 
  scale_fill_brewer(palette = "Pastel1") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 

ggsave("zc4688/honours/analyses/chordates/figures/igsrepeats.absolute.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/igsrepeats.absolute.png", height = 8, width = 12)



trf_ranges <- trf_ranges %>% 
  left_join(structuremetadata) %>% 
  mutate(`Simple repeats` = `Simple repeats`/IGS_length)

rm_ranges <- rm_ranges %>% 
  left_join(structuremetadata) %>% 
  mutate(TEs = TEs/IGS_length)

classifications %>% 
  left_join(rm_ranges %>% select(TEs, `Sample ID`)) %>% 
  left_join(trf_ranges %>% select(`Simple repeats`, `Sample ID`)) %>% 
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>% 
  mutate(prop_repeat = TEs + `Simple repeats`) %>% 
  select(group, TEs, `Simple repeats`) %>% 
  filter(group %in% common_class) %>% 
  pivot_longer(cols = c(TEs, `Simple repeats`), names_to = "variable", values_to = "value") %>%
  ggplot() + 
  geom_boxplot(aes(x = group, y = value, fill = variable), notch = F) + 
  theme_bw() + 
  labs(x = NULL, y = "Proportion of IGS", fill = NULL) + 
  scale_fill_brewer(palette = "Pastel1") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 

ggsave("zc4688/honours/analyses/chordates/figures/igsrepeats.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/igsrepeats.png", height = 8, width = 12)


####################### 28S repeats

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
  filter(start > `28S_rRNA_Start` & start < `28S_rRNA_End`) %>% 
  group_by(`Sample ID`) %>% 
  summarise(`Simple repeats` = sum(end - start), .groups = "drop") %>% 
  unique()

rm_ranges <- structuremetadata %>% 
  left_join(rm_ranges) %>% 
  filter(start > `28S_rRNA_Start` & start < `28S_rRNA_End`) %>% 
  group_by(`Sample ID`) %>% 
  summarise(TEs = sum(end - start), .groups = "drop") %>% 
  unique()



classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>% 
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>% 
  select(group, TEs, `Simple repeats`) %>% 
  filter(group %in% common_class) %>% 
  pivot_longer(cols = c(TEs, `Simple repeats`), names_to = "variable", values_to = "value") %>% view %>% 
  ggplot() + 
  geom_boxplot(aes(x = group, y = value, fill = variable), notch = F) + 
  theme_bw() + 
  labs(x = NULL, y = "Length of repeats in 28S", fill = NULL) + 
  scale_fill_brewer(palette = "Pastel1") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 

ggsave("zc4688/honours/analyses/chordates/figures/twoeightrepeats.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/twoeightrepeats.png", height = 8, width = 12)
classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(`Length of repeats in 28S` = `Simple repeats` + TEs, `Length of 28S without TEs` = `28S_length` - TEs, `Length of 28S without repeats` = `28S_length` - `Length of repeats in 28S`) %>% 
  filter(group %in% common_class)  %>% 
  select(`Length of repeats in 28S`, `Length of 28S without TEs`, `Length of 28S without repeats`, `28S_length`, group, `Sample ID`) %>% 
  pivot_longer(-c(group, `Sample ID`)) %>% 
  ggplot(aes(y = value, fill = group)) + 
  geom_boxplot(coef = 1.5, outliers = FALSE) + 
  labs(y = "Length (bp)") +
  scale_fill_manual(values = colours, name = "Taxonomic group") +
  facet_wrap(~name) + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_blank())
ggsave("zc4688/honours/analyses/chordates/figures/twoeightrepeats.summary.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/twoeightrepeats.summary.png", height = 8, width = 12)

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
  filter(start > `18S_rRNA_End` & start < `28S_rRNA_Start`) %>% 
  group_by(`Sample ID`) %>% 
  summarise(`Simple repeats` = sum(end - start), .groups = "drop") %>% 
  unique()

rm_ranges <- structuremetadata %>% 
  left_join(rm_ranges) %>% 
  filter(start > `18S_rRNA_End` & start < `28S_rRNA_Start`) %>% 
  group_by(`Sample ID`) %>% 
  summarise(TEs = sum(end - start), .groups = "drop") %>% 
  unique()



classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>% 
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>% 
  select(group, TEs, `Simple repeats`) %>% 
  filter(group %in% common_class) %>% 
  pivot_longer(cols = c(TEs, `Simple repeats`), names_to = "variable", values_to = "value") %>% view %>% 
  ggplot() + 
  geom_boxplot(aes(x = group, y = value, fill = variable), notch = F, outliers = F) + 
  theme_bw() + 
  labs(x = NULL, y = "Length of repeats in ITS1", fill = NULL) + 
  scale_fill_brewer(palette = "Pastel1") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 

ggsave("zc4688/honours/analyses/chordates/figures/its1repeats.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/its1repeats.png", height = 8, width = 12)


classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(total = `Simple repeats` + TEs) %>% 
  filter(group %in% common_class) %>% pairwise_wilcox_test(total ~ group, p.adj.method = "BH") %>% view

classifications %>% 
  left_join(rm_ranges) %>% 
  left_join(trf_ranges) %>%  
  left_join(structuremetadata) %>%
  replace_na(list(TEs = 0, `Simple repeats` = 0)) %>%
  mutate(`Length of repeats in ITS1` = `Simple repeats` + TEs, `Length of ITS1 without TEs` = ITS1 - TEs, `Length of ITS1 without repeats` = ITS1 - `Length of repeats in ITS1`) %>% 
  filter(group %in% common_class)  %>% 
  select(`Length of repeats in ITS1`, `Length of ITS1 without TEs`, `Length of ITS1 without repeats`, ITS1, group, `Sample ID`) %>% 
  pivot_longer(-c(group, `Sample ID`)) %>% 
  ggplot(aes(y = value, fill = group)) + 
  geom_boxplot(coef = 1.5, outliers = FALSE) + 
  labs(y = "Length (bp)") +
  scale_fill_manual(values = colours, name = "Taxonomic group") +
  facet_wrap(~name) + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x=element_blank())
ggsave("zc4688/honours/analyses/chordates/figures/its1repeats.summary.pdf", height = 8, width = 12)
ggsave("zc4688/honours/analyses/chordates/figures/its1repeats.summary.png", height = 8, width = 12)



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
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + geom_facet(
  data = structureplot,
  mapping = aes(x = 0, xend = ifelse(`Unit length` < 60000, `Unit length`, 60000)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey50",
  alpha = 0.3
) + geom_facet(
  data = plotdata,
  mapping = aes(x = ifelse(`Unit length` >= 60000, 60000, NA)),
  geom = geom_point,
  color = "grey50",
  panel = " ",
  size = 1
)+
  new_scale_color() + geom_facet(
    data = structureplot,
    mapping = aes(x = Start, xend = End, color = Type),
    geom = geom_segment,       # Horizontal bar chart
    panel = " "
  )  +  
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  new_scale_color() +
  geom_facet(
    data = repeats %>% filter(end < 60000),
    mapping = aes(x = start, xend = end, color = repeatgroup),
    geom = geom_segment, 
    panel = " ")+ scale_color_manual(
      values = repeat_colors,
      name = "Repeat Type"
    ) +
  scale_y_discrete() +  
  theme_tree2() + scale_alpha(range = c(0.1, 1))+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off') + xlim_expand(c(0, 23), "Tree")

facet_widths(p, widths = c(1, 3))



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
         End = ifelse(strand == "+", envend - morphstart, morphend - envend + 1) )  %>%
  mutate(Species = str_replace_all(Species, "_", " ")) %>% 
  filter(Species %in% common_organisms)  %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())




p <- ggtree(tree) %<+% (primarystructure %>% select(label, group)) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + geom_facet(
  data = primarystructure,
  mapping = aes(x = 0, xend = ifelse(`Unit length` < 60000, `Unit length`, 60000)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey50",
  alpha = 0.3
) + geom_facet(
  data = plotdata,
  mapping = aes(x = ifelse(`Unit length` >= 60000, 60000, NA)),
  geom = geom_point,
  color = "grey50",
  panel = " ",
  size = 1
)+
  new_scale_color() + geom_facet(
    data = primarystructure,
    mapping = aes(x = Start, xend = End, color = type),
    geom = geom_segment,       # Horizontal bar chart
    panel = " "
  )  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 23), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "rRNA")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

facet_widths(p, widths = c(1, 3))


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
    color = "grey50",
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
    color = "grey50",
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



############################### Metadata table

metadata <- metadata %>% 
  left_join(structuremetadata)

metadata <- metadata %>% 
  select(-"18S_rRNA_Start", -"28S_rRNA_Start", -"5_8S_rRNA_Start", -"18S_rRNA_End", -"28S_rRNA_End", -"5_8S_rRNA_End", -(11:29)) %>% 
  dplyr::rename("Number of morphs" = "rDNA_details.Number of morphs")

write.table(metadata, "zc4688/honours/metadata/metadata.tsv", sep = "\t", quote = F, row.names = F)


################################Conservation
shannon_entropy <- function(column) {
  column <- column[column != "-"]
  freqs <- table(column) / length(column)  # Frequency of each base
  -sum(freqs * log2(freqs), na.rm = TRUE)  # Shannon entropy formula
}

fiveeight <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/fiveeight_clusters/cleaned.alignment.txt")


# Convert alignment to matrix
fiveeight <- as.data.frame(as.matrix(fiveeight))
fiveeight <- fiveeight[!rownames(fiveeight) %in% "GCF_015220235_NC_051310.1", ]
#overall:
median(apply(fiveeight, 2, shannon_entropy))


# Apply function to each column
fiveeight <- apply(fiveeight, 2, shannon_entropy)


twoeight <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/twoeight_clusters/cleaned.alignment.txt")
twoeight <- as.data.frame(as.matrix(twoeight))
twoeight <- twoeight[!rownames(twoeight) %in% "GCF_015220235_NC_051310.1", ]

#overall:
median(apply(twoeight, 2, shannon_entropy))


# Apply function to each column
twoeight <- apply(twoeight, 2, shannon_entropy)


eighteen <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/eighteen_clusters/cleaned.alignment.txt")

# Convert alignment to matrix
eighteen <- as.data.frame(as.matrix(eighteen))
eighteen <- eighteen[!rownames(eighteen) %in% "GCF_015220235_NC_051310.1", ]

#overall:
median(apply(eighteen, 2, shannon_entropy))


# Apply function to each column
eighteen <- apply(eighteen, 2, shannon_entropy)



eighteen <- data.frame(Position = 1:length(eighteen), Score = eighteen)
twoeight <- data.frame(Position = 1:length(twoeight), Score = twoeight)
fiveeight <- data.frame(Position = 1:length(fiveeight), Score = fiveeight)

fiveeight <- fiveeight %>% mutate(gene = "5.8S")
twoeight <- twoeight %>% mutate(gene = "28S")
eighteen <- eighteen %>% mutate(gene = "18S")
tmp <- bind_rows(fiveeight, eighteen, twoeight)
ggplot(tmp, aes(x = Score, fill = gene)) + 
  geom_density(alpha = 0.6) + theme_bw() + 
  scale_fill_manual(values = c("18S" = "#F75F86", "28S" = "cornflowerblue", "5.8S" = "#50C878"), name = "rRNA") + 
  labs(x = "Shannon entropy",
       y = "Density")
ggsave("zc4688/honours/analyses/chordates/figures/entropy.pdf", height = 6, width = 10)
ggsave("zc4688/honours/analyses/chordates/figures/entropy.png", height = 6, width = 10)

ggplot()+  
  #geom_line(data = tmp, aes(x = Position, y = Score), linewidth = 0.4) + 
  geom_point(data = eighteen, aes(x = Position, y = Score, color = Score), size = 2) + 
  geom_point(data = eighteen %>% filter(Score == 0), aes(x = Position, y = 00), size = 2, color = "black") + 
  scale_color_gradient2(low = "forestgreen", mid="gold", high = "red", na.value = "grey", midpoint = 1.16) + theme_bw() +
  labs(y = "Shannon entropy")
ggsave("zc4688/honours/analyses/chordates/figures/eighteen.entropy.pdf", height = 6, width = 10)

ggplot()+  
  #geom_line(data = tmp, aes(x = Position, y = Score), linewidth = 0.4) + 
  geom_point(data = fiveeight, aes(x = Position, y = Score, color = Score), size = 2) + 
  geom_point(data = fiveeight %>% filter(Score == 0), aes(x = Position, y = 00), size = 2, color = "black") + 
  scale_color_gradient2(low = "forestgreen", mid="gold", high = "red", na.value = "grey", midpoint = 1.16) + theme_bw() +
  labs(y = "Shannon entropy")
ggsave("zc4688/honours/analyses/chordates/figures/fiveeight.entropy.pdf", height = 6, width = 10)

ggplot()+  
  #geom_line(data = tmp, aes(x = Position, y = Score), linewidth = 0.4) + 
  geom_point(data = twoeight, aes(x = Position, y = Score, color = Score), size = 2) + 
  geom_point(data = twoeight %>% filter(Score == 0), aes(x = Position, y = 00), size = 2, color = "black") + 
  scale_color_gradient2(low = "forestgreen", mid="gold", high = "red", na.value = "grey", midpoint = 1.16) + theme_bw() +
  labs(y = "Shannon entropy")
ggsave("zc4688/honours/analyses/chordates/figures/twoeight.entropy.pdf", height = 6, width = 10)


############################### invariant pos

fiveeight <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/fiveeight_clusters/final_alignment_with_all.fasta")


# Convert alignment to matrix
fiveeight <- as.data.frame(as.matrix(fiveeight))
fiveeight <- fiveeight[!rownames(fiveeight) %in% "GCF_015220235_NC_051310.1", ]

hg002_row <- grep("hg002", rownames(fiveeight))

fiveeight <- fiveeight %>%
  select(which(fiveeight[hg002_row, ] != "-"))

# Apply function to each column
fiveeight <- apply(fiveeight, 2, shannon_entropy)


twoeight <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/twoeight_clusters/final_alignment_with_all.fasta")
twoeight <- as.data.frame(as.matrix(twoeight))
twoeight <- twoeight[!rownames(twoeight) %in% "GCF_015220235_NC_051310.1", ]



hg002_row <- grep("hg002", rownames(twoeight))

twoeight <- twoeight %>%
  select(which(twoeight[hg002_row, ] != "-"))

# Apply function to each column
twoeight <- apply(twoeight, 2, shannon_entropy)


eighteen <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/eighteen_clusters/final_alignment_with_all.fasta")

# Convert alignment to matrix
eighteen <- as.data.frame(as.matrix(eighteen))
eighteen <- eighteen[!rownames(eighteen) %in% "GCF_015220235_NC_051310.1", ]


hg002_row <- grep("hg002", rownames(eighteen))

eighteen <- eighteen %>%
  select(which(eighteen[hg002_row, ] != "-"))

# Apply function to each column
eighteen <- apply(eighteen, 2, shannon_entropy)



eighteen <- data.frame(Position = 1:length(eighteen), Score = eighteen)
twoeight <- data.frame(Position = 1:length(twoeight), Score = twoeight)
fiveeight <- data.frame(Position = 1:length(fiveeight), Score = fiveeight)

eighteen <- eighteen %>% mutate(Score = ifelse(Score == 0, 0, 1))
twoeight <- twoeight %>% mutate(Score = ifelse(Score == 0, 0, 1))
fiveeight <- fiveeight %>% mutate(Score = ifelse(Score == 0, 0, 1))


write.csv(eighteen, 
          "zc4688/honours/analyses/chordates/new/msa/eighteen.pos.status.csv", quote = F, row.names = F)

write.csv(twoeight,"zc4688/honours/analyses/chordates/new/msa/twoeight.pos.status.csv", quote = F, row.names = F)
write.csv(fiveeight,"zc4688/honours/analyses/chordates/new/msa/fiveeight.pos.status.csv", quote = F, row.names = F)




############################### 18S tree

tree <- read.tree("/g/data/te53/zc4688/honours/trees/species_ncbi.phy")

nodes <- (tree$node.label)

#tree$edge.length <- NULL
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

tree$tip.label <- str_replace_all(tree$tip.label, "'", "")

genetree <- read.tree("/g/data/te53/zc4688/honours/analyses/chordates/new/msa/eighteen_clusters/eighteen.treefile")

treenames <- as.data.frame(genetree$tip.label)


treenames <- treenames %>% 
  mutate(`Sample ID` = str_extract(`genetree$tip.label`, "^GCA_\\d+|^GCF_\\d+|hg002"))

treenames <- treenames %>% left_join(classifications %>% select(`Sample ID`, Species))

genetree$tip.label <- as.character(treenames$Species)
common_organisms <- intersect(plotdata$label, genetree$tip.label)

genetree <- ape::drop.tip(genetree, setdiff(genetree$tip.label, common_organisms))


genetree <- ladderize(genetree)


tmp <- plotdata %>% select(label, group)

ggtree(genetree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") +
  coord_cartesian(clip = 'off') 

ggsave("zc4688/honours/analyses/chordates/figures/eighteentree_branches.pdf", height = 20, width = 15)

common_cotree <- intersect(genetree$tip.label, tree$tip.label)

genetree <- ape::drop.tip(genetree, setdiff(genetree$tip.label, common_cotree))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_cotree))
A <- cbind(tree$tip.label, tree$tip.label)
x <- A %>% as.data.frame(.)  %>% dplyr::rename("Species" = "V1") %>% left_join(classifications)
x$color <- colours[x$group]

tree <- ladderize(tree)


dist.topo(tree, genetree)




tree_un <- unroot(tree)
dist.topo(tree_un, genetree)
treedist(tree_un, genetree)

genetree$edge.length <- NULL
tree$edge.length <- NULL

genetree <- root(genetree, outgroup = c("Branchiostoma lanceolatum", "B. floridae x B. belcheri", "B. floridae x B. japonicum"), resolve.root = T)

cophyloplot(tree, genetree, assoc = A, show.tip.label = F, space=500, col = x$color) 

ggtree(genetree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") +
  coord_cartesian(clip = 'off') 
ggsave("zc4688/honours/analyses/chordates/figures/eighteentree.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/eighteentree.png", height = 20, width = 15)

############################### 28S tree

tree <- read.tree("/g/data/te53/zc4688/honours/trees/species_ncbi.phy")

nodes <- (tree$node.label)

#tree$edge.length <- NULL
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

tree$tip.label <- str_replace_all(tree$tip.label, "'", "")

genetree <- read.tree("/g/data/te53/zc4688/honours/analyses/chordates/new/msa/twoeight_clusters/twoeight.treefile")

treenames <- as.data.frame(genetree$tip.label)


treenames <- treenames %>% 
  mutate(`Sample ID` = str_extract(`genetree$tip.label`, "^GCA_\\d+|^GCF_\\d+|hg002"))

treenames <- treenames %>% left_join(classifications %>% select(`Sample ID`, Species))

genetree$tip.label <- as.character(treenames$Species)
common_organisms <- intersect(plotdata$label, genetree$tip.label)

genetree <- ape::drop.tip(genetree, setdiff(genetree$tip.label, common_organisms))


genetree <- ladderize(genetree)


tmp <- plotdata %>% select(label, group)

ggtree(genetree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") +
  coord_cartesian(clip = 'off') 

ggsave("zc4688/honours/analyses/chordates/figures/twoeighttree_branches.pdf", height = 20, width = 15)

common_cotree <- intersect(genetree$tip.label, tree$tip.label)

genetree <- ape::drop.tip(genetree, setdiff(genetree$tip.label, common_cotree))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_cotree))
A <- cbind(tree$tip.label, tree$tip.label)
x <- A %>% as.data.frame(.)  %>% dplyr::rename("Species" = "V1") %>% left_join(classifications)
x$color <- colours[x$group]

tree <- ladderize(tree)


dist.topo(tree, genetree)




tree_un <- unroot(tree)
dist.topo(tree_un, genetree)
treedist(tree_un, genetree)

genetree$edge.length <- NULL
tree$edge.length <- NULL

genetree <- root(genetree, outgroup = c("Branchiostoma lanceolatum", "B. floridae x B. belcheri", "B. floridae x B. japonicum"), resolve.root = T)

cophyloplot(tree, genetree, assoc = A, show.tip.label = F, space=500, col = x$color) 

ggtree(genetree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") +
  coord_cartesian(clip = 'off') 
ggsave("zc4688/honours/analyses/chordates/figures/twoeighttree.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/twoeighttree.png", height = 20, width = 15)

############################### 5.8S tree

tree <- read.tree("/g/data/te53/zc4688/honours/trees/species_ncbi.phy")

nodes <- (tree$node.label)

#tree$edge.length <- NULL
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

tree$tip.label <- str_replace_all(tree$tip.label, "'", "")

genetree <- read.tree("/g/data/te53/zc4688/honours/analyses/chordates/new/msa/fiveeight_clusters/fiveeight.treefile")

treenames <- as.data.frame(genetree$tip.label)


treenames <- treenames %>% 
  mutate(`Sample ID` = str_extract(`genetree$tip.label`, "^GCA_\\d+|^GCF_\\d+|hg002"))

treenames <- treenames %>% left_join(classifications %>% select(`Sample ID`, Species))

genetree$tip.label <- as.character(treenames$Species)
common_organisms <- intersect(plotdata$label, genetree$tip.label)

genetree <- ape::drop.tip(genetree, setdiff(genetree$tip.label, common_organisms))


genetree <- ladderize(genetree)


tmp <- plotdata %>% select(label, group)

ggtree(genetree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") +
  coord_cartesian(clip = 'off') 

ggsave("zc4688/honours/analyses/chordates/figures/fiveeighttree_branches.pdf", height = 20, width = 15)

common_cotree <- intersect(genetree$tip.label, tree$tip.label)

genetree <- ape::drop.tip(genetree, setdiff(genetree$tip.label, common_cotree))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_cotree))
A <- cbind(tree$tip.label, tree$tip.label)
x <- A %>% as.data.frame(.)  %>% dplyr::rename("Species" = "V1") %>% left_join(classifications)
x$color <- colours[x$group]

tree <- ladderize(tree)


dist.topo(tree, genetree)




tree_un <- unroot(tree)
dist.topo(tree_un, genetree)
treedist(tree_un, genetree)

genetree$edge.length <- NULL
tree$edge.length <- NULL

genetree <- root(genetree, outgroup = c("Branchiostoma lanceolatum", "B. floridae x B. belcheri", "B. floridae x B. japonicum"), resolve.root = T)

cophyloplot(tree, genetree, assoc = A, show.tip.label = F, space=500, col = x$color) 

ggtree(genetree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") +
  coord_cartesian(clip = 'off') 
ggsave("zc4688/honours/analyses/chordates/figures/fiveeightntree.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/fiveeighttree.png", height = 20, width = 15)

##############

tmp <-  read.table("zc4688/honours/analyses/chordates/motif/its1_single/meme.txt", sep = ",") %>% 
      dplyr::rename(raw = V1)


#need to improve parsing (some positions etc not working). perhaps include length. create plot showing all significant motifs in diff colours

#include analysis of no. of each motif
motifs.tmp <- tmp %>%
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
    


motifs <- NULL

for (i in seq(1, nrow(motifs.tmp), 2)) {
  start.row <- motifs.tmp$row.num[i]
  
  end.row <- motifs.tmp$row.num[i + 1]
  
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
  
  motifs <- rbind(motifs, df.motif.temp)
}




motifs <- motifs  %>% 
  mutate(`Sample ID` = str_extract(input.seq.name, "^GCA_\\d+|^GCF_\\d+")) %>% 
  left_join(classifications) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything()) %>% 
  filter(label %in% common_organisms) %>% 
  left_join(structuremetadata %>% select(`Sample ID`, ITS1, ITS2, IGS_length)) %>% 
  group_by(`Sample ID`, motif.num) %>% 
  mutate(number = n()) %>%
  ungroup()


tmp <- motifs %>% filter(motif.num == "Motif.1") %>% select(label, input.seq.motif, `Sample ID`, group) 
dna <- DNAStringSet(tmp$input.seq.motif)
motif_mat <- as.matrix(dna)


motif_df <- data.frame(
  label = rep(tmp$label, each = width(dna)[1]),
  tax_group = rep(tmp$group, each = width(dna)[1]),
  position = rep(1:width(dna)[1], times = length(dna)),
  base = as.vector(t(motif_mat))
)

motif_df <- motif_df %>%
  mutate(
    tax_group = factor(tax_group, levels = unique(tmp$group)),
    label = factor(label, levels = tmp %>%
                         arrange(group, label) %>%
                         pull(label))
  )

ggplot(motif_df, aes(x = position, y = label, fill = base)) +
  geom_tile() +
  scale_fill_manual(values = c(A = "#1b9e77", C = "#d95f02", G = "#7570b3", T = "#e7298a")) +
  facet_grid(rows = vars(tax_group), scales = "free_y", space = "free_y") +
  theme_minimal(base_size = 10) +
  labs(x = "Position (bp)", y = "Sample", fill = "Base") +
  theme(
    axis.text.y = element_text(size = 6),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold")
  )



p <- ggtree(tree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") 


p <- p + new_scale_color() + 
  geom_facet(
    data = motif_df,
    mapping = aes(x = position, fill = base),
    geom = geom_tile,
    panel = "Motif",
    size = 2 
  ) + scale_y_discrete()  +
  theme_tree2() + coord_cartesian(clip = 'off') +
  scale_fill_manual(values = c(A = "#1b9e77", C = "#d95f02", G = "#7570b3", T = "#e7298a"), name = "Base") 

facet_widths(p, widths = c(1, 3))

ggsave("zc4688/honours/analyses/chordates/figures/its1.motif.1.seq.pdf", height = 20, width = 15)


motifs %>% group_by(motif.num, group) %>% 
  summarise(n = n()) %>% 
  left_join(classifications %>% 
              group_by(group) %>% 
              summarise(total = n())) %>% mutate(n = n/total) %>% 
  ggplot(aes(x = group, y = n, fill = group)) + geom_boxplot() + facet_wrap(~motif.num)





#for last
aves <- motifs %>% 
  mutate(num = as.numeric(str_replace_all(motif.num, "Motif.", ""))) %>% 
  filter(num %in% c(1,2,4,5,6, 10, 11,16,17,19,20,21,25,27,28,32))
motif_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(its1$motif.num)))

aves_labels <- classifications %>% filter(group == "Aves") %>%  pull(Species)

# Prune tree to Aves only
aves_tree <- drop.tip(tree, setdiff(tree$tip.label, aves_labels))

p <- ggtree(aves_tree) %<+% classifications +
  geom_tiplab(size=3) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group") 
p <- p + new_scale_color() + 
  geom_facet(
    data = structuremetadata %>% left_join(classifications) %>% dplyr::rename("label" = "Species") %>% 
      relocate(label, .before = everything()),
    mapping = aes(x = ITS1),
    geom = geom_col,
    color = "grey",
    fill = "grey",
    panel = "Motif",
    size = 2 
  ) + geom_facet(
    data = aves,
    mapping = aes(x = as.numeric(input.seq.pos), xend = as.numeric(input.seq.pos) + 50, color = motif.num),
    geom = geom_segment,       # Horizontal bar chart
    panel = "Motif",
    size = 2 
  )     +
    scale_y_discrete() +  scale_color_manual(values = motif_colors) +
    theme_tree2() + xlim_expand(c(0, 25), "Tree") +
    labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 
  
p
ggsave("zc4688/honours/analyses/chordates/figures/its1.ave.motifs.pdf", height = 20, width = 15)



p <- ggtree(tree) %<+% classifications +
  geom_tippoint(aes(color = group), size = 2) +
  scale_color_manual(values = colours, name = "Taxonomic Group") 
p <- p  +
  new_scale_color() + 
  geom_facet(
    data = motif_df,
    mapping = aes(x = position, fill = base),
    geom = geom_tile,
    panel = "Motif sequence",
    size = 2 
  )  + new_scale_color() + 
  geom_facet(
    data = structuremetadata %>% left_join(classifications) %>% dplyr::rename("label" = "Species") %>% 
      relocate(label, .before = everything()),
    mapping = aes(x = ITS1 ),
    geom = geom_col,
    color = "grey",
    panel = "Position of motif in ITS1",
    size = 2 
  ) + geom_facet(
    data = motifs %>% filter(motif.num == "Motif.1")  ,
    mapping = aes(x = as.numeric(input.seq.pos), xend = as.numeric(input.seq.pos) + 50),
    geom = geom_segment,       # Horizontal bar chart
    color = "black",
    panel = "Position of motif in ITS1",
    size = 1 
  )   + scale_y_discrete()  +
  theme_tree2() + coord_cartesian(clip = 'off') +
  scale_fill_manual(values = c(A = "#1b9e77", C = "#d95f02", G = "#7570b3", T = "#e7298a"), name = "Base") +
  scale_y_discrete() +  scale_color_brewer(palette = "Dark2", name = "Motif name") +
  theme_tree2() + xlim_expand(c(0, 25), "Tree") +
  labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(1, 3, 2))


ggsave("zc4688/honours/analyses/chordates/figures/its1.motif.1.pdf", height = 20, width = 15)




p <- ggtree(tree) %<+% classifications +
  geom_tippoint(aes(color = group), size = 2) +
  scale_color_manual(values = colours, name = "Taxonomic Group") 
p <- p + new_scale_color() + 
  geom_facet(
    data = structuremetadata %>% left_join(classifications) %>% dplyr::rename("label" = "Species") %>% 
      relocate(label, .before = everything()),
    mapping = aes(x = ITS1 ),
    geom = geom_col,
    color = "grey",
    panel = "Position in ITS1",
    size = 2 
  ) + geom_facet(
    data = motifs %>% 
      mutate(num = as.numeric(str_replace_all(motif.num, "Motif.", ""))) %>% 
      filter(num %in% c(1, 7, 9)),
    mapping = aes(x = as.numeric(input.seq.pos), xend = as.numeric(input.seq.pos) + 50, color = factor(num)),
    geom = geom_segment,       # Horizontal bar chart
    panel = "Position in ITS1",
    size = 1 
  )    +
  scale_y_discrete() +  scale_color_brewer(palette = "Dark2", name = "Motif name") +
  theme_tree2() + xlim_expand(c(0, 25), "Tree") +
  labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(1, 3))
ggsave("zc4688/honours/analyses/chordates/figures/its1.motifs.conserved.pdf", height = 20, width = 15)

######################### ITS2
#include analysis of no. of each motif
motifs.tmp <- tmp %>%
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



motifs <- NULL

for (i in seq(1, nrow(motifs.tmp), 2)) {
  start.row <- motifs.tmp$row.num[i]
  
  end.row <- motifs.tmp$row.num[i + 1]
  
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
  
  motifs <- rbind(motifs, df.motif.temp)
}




motifs <- motifs  %>% 
  mutate(`Sample ID` = str_extract(input.seq.name, "^GCA_\\d+|^GCF_\\d+")) %>% 
  left_join(classifications) %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything()) %>% 
  filter(label %in% common_organisms) %>% 
  left_join(structuremetadata %>% select(`Sample ID`, ITS1, ITS2, IGS_length)) %>% 
  group_by(`Sample ID`, motif.num) %>% 
  mutate(number = n()) %>%
  ungroup()




p <- ggtree(tree) %<+% classifications +
  geom_tippoint(aes(color = group), size = 2) +
  scale_color_manual(values = colours, name = "Taxonomic Group") 
p <- p + new_scale_color() + 
  geom_facet(
    data = structuremetadata %>% left_join(classifications) %>% dplyr::rename("label" = "Species") %>% 
      relocate(label, .before = everything()),
    mapping = aes(x = ITS2 ),
    geom = geom_col,
    color = "grey",
    panel = "Position in ITS2",
    size = 2 
  ) + geom_facet(
    data = motifs %>% 
      mutate(num = as.numeric(str_replace_all(motif.num, "Motif.", ""))) %>% 
      filter(num %in% c(5, 7, 10, 1, 2, 3)),
    mapping = aes(x = as.numeric(input.seq.pos), xend = as.numeric(input.seq.pos) + 50, color = factor(num)),
    geom = geom_segment,       # Horizontal bar chart
    panel = "Position in ITS2",
    size = 1 
  )    +
  scale_y_discrete() +  scale_color_brewer(palette = "Dark2", name = "Motif name") +
  theme_tree2() + xlim_expand(c(0, 25), "Tree") +
  labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(1, 3))
ggsave("zc4688/honours/analyses/chordates/figures/its2.motifs.conserved.pdf", height = 20, width = 15)


################ GC content
tmp <- read.delim("zc4688/honours/analyses/chordates/new/msa/twoeight_gc.txt", header = F)

tmp <- tmp %>% 
  mutate(`Sample ID` = str_extract(V1, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  mutate(start = str_extract(V1, "(?<=sliding:)\\d+")) %>% 
  left_join(classifications)

tmp <- tmp %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())


p <- ggtree(tree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")
p + new_scale_color() + geom_facet(
  data = tmp,
  mapping = aes(x = as.numeric(start), color = V2),
  geom = geom_point,       # Horizontal bar chart
  panel = "Position in 288S"
) +
  scale_y_discrete()+scale_color_viridis_b(option = "viridis", n.breaks=10) + labs(color = "GC content (%)")+ coord_cartesian(clip = 'off')+ xlim_expand(c(0, 25), "Tree")


ggsave("zc4688/honours/analyses/chordates/figures/twoeight_gc.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/twoeight_gc.png", height = 20, width = 15)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = n, fill = group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = colours) + 
  theme_bw() + 
  labs(x = NULL, y = "Number of windows with >80% GC content", fill = "Taxonomic Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 

ggsave("zc4688/honours/analyses/chordates/figures/gc.length.pdf", height = 6, width = 10)
ggsave("zc4688/honours/analyses/chordates/figures/gc.length.png", height = 6, width = 10)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(n ~ group, p.adj.method = "BH") %>% view

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = max, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "Maximum GC % in 28S", fill = "Taxonomic Group")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 
ggsave("zc4688/honours/analyses/chordates/figures/gc.max.pdf", height = 6, width = 10)
ggsave("zc4688/honours/analyses/chordates/figures/gc.max.png", height = 6, width = 10)

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(max ~ group, p.adj.method = "BH") %>% view


classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  left_join(
    structuremetadata
  ) %>%
  mutate(remain = `28S_length` - n*50) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = remain, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "28S length without GC rich regions", fill = "Taxonomic Group") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 
ggsave("zc4688/honours/analyses/chordates/figures/gc.removed.pdf", height = 6, width = 10)
ggsave("zc4688/honours/analyses/chordates/figures/gc.removed.png", height = 6, width = 10)


tmp <- read.delim("zc4688/honours/analyses/chordates/motif/its1_gc.txt", header = F)

tmp <- tmp %>% 
  mutate(`Sample ID` = str_extract(V1, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  mutate(start = str_extract(V1, "(?<=sliding:)\\d+")) %>% 
  left_join(classifications)

tmp <- tmp %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())


p <- ggtree(tree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")
p + new_scale_color() + geom_facet(
  data = tmp,
  mapping = aes(x = as.numeric(start), color = V2),
  geom = geom_point,       # Horizontal bar chart
  panel = "Position in 288S"
) +
  scale_y_discrete()+scale_color_viridis_b(option = "viridis", n.breaks=10) + labs(color = "GC content (%)")+ coord_cartesian(clip = 'off')+ xlim_expand(c(0, 25), "Tree")


ggsave("zc4688/honours/analyses/chordates/figures/its1_gc.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1_gc.png", height = 20, width = 15)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = n, fill = group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = colours) + 
  theme_bw() + 
  labs(x = NULL, y = "Number of windows with >80% GC content", fill = "Taxonomic Group")
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.length.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.length.png", height = 20, width = 15)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(n ~ group, p.adj.method = "BH") %>% view

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = max, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "Maximum GC % in ITS1", fill = "Taxonomic Group")
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.max.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.max.png", height = 20, width = 15)

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(max ~ group, p.adj.method = "BH") %>% view


classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  left_join(
    structuremetadata
  ) %>%
  mutate(remain = ITS1 - n*50) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = remain, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "ITS1 length without GC rich regions", fill = "Taxonomic Group") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.removed.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.removed.png", height = 20, width = 15)



