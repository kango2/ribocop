.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))
library(ggtree)
library(tidyverse)
library(ape)
library(jsonlite)
library(GenomicRanges)
library(taxize)
library(RColorBrewer)
library(ggnewscale)
############################### Figure 1: length tree
directory <- "/g/data/te53/zc4688/honours/ribocop/results/new"

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
new <- bind_rows(data_list)
tmp <- combined_df


combined_df <- combined_df %>% 
  filter(is.na(`Errors.Alignment filtering`) & is.na(`Errors.Morph identification`) & is.na(`Errors.Non-ambiguous morph identification`)) %>% 
  filter(!is.na(`rDNA_details.Number of primary alignments`)) %>% 
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
  mutate(`Minimum unut length` = as.numeric(`Minimum unit length`)) %>% 
  mutate(`Maximum unit length` = as.numeric(`Maximum unit length`)) %>% 
  mutate(CN = round(`Total length`/`Unit length`, digits=2))


lengths <- combined_df

lengths <- lengths %>% 
  filter(!is.na(`Unit length`))

tree <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (5).nwk") #ORIGINAL tree before filtering, just for checking. most current results are in refiltered directory (with current filtering), all has more because filtering was lest strict

treenames <- as.data.frame(tree$tip.label)
colnames(treenames) <- c("Timetreename")
replacementnames <- read.delim("zc4688/honours/ribocop/speciestimetree.tsv", header = T, sep = ",")

replacementnames <- replacementnames %>% mutate(Species = str_replace_all(Species, " ", "_"), Timetreename = str_replace_all(Replacement, " ", "_") ) %>% select(-Replacement)

treenames <- treenames %>% left_join(replacementnames) %>% mutate(Species = ifelse(is.na(Species), Timetreename, Species)) %>% select(Species)
tree$tip.label <- as.character(treenames$Species)

plotdata <- lengths %>% select(Species, `Unit length`, `CN`, `Total length`, `Minimum unit length`, `Maximum unit length`, `Sample ID`, Contigs) %>% 
  group_by(Species) %>% dplyr::summarize(across(everything(), ~ dplyr::last(.))) %>% 
  filter(`Unit length` < 120000)





plotdata <- plotdata %>% 
  dplyr::rename("label" = "Species")

metadata <- plotdata
tree$tip.label <- str_replace_all(tree$tip.label, "'", "")

common_organisms <- intersect(plotdata$label, tree$tip.label)


plotdata <- plotdata %>%
  filter(label %in% common_organisms)



#classifications <- classification(plotdata$label, db="ncbi")
#classifications <- cbind(classifications)
#classifications <- classifications %>% 
#select(kingdom, phylum, subphylum, superclass, class, subclass, superorder, order, suborder, family, subfamily, genus, species, query)
#write.table(classifications, "zc4688/honours/ribocop/classifications.tsv", sep = "\t", quote = F, row.names = F)
classifications <- read.table("zc4688/honours/ribocop/classifications.tsv", sep = "\t", header = TRUE)




classifications <- classifications %>% 
  dplyr::select(class, order, species) %>% 
  dplyr::rename("label" = "species") 

classifications <- classifications %>% 
  group_by(class, order, label) %>% summarize(across(everything(), ~ last(.))) #Take refseq if both refseq and genbank are available
plotdata <- plotdata %>% mutate(label = str_replace_all(label, "_", " "))
plotdata <- left_join(plotdata, classifications, by="label")


tree$edge.length <- NULL
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))


tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")

plotdata <- plotdata %>% mutate(class = ifelse(is.na(class), "Other", class)) %>% relocate(label, .before = everything())


#Manual fixes for species where classification hasnt worked
plotdata <- plotdata %>%
  mutate(
    class = case_when(
      label %in% c("Arripis georgianus", "Cyprinodon nevadensis mionectes", "Gymnocypris eckloni scoliostomus", "Distoechodon macrophthalmus",
                   "Pempheris klunzingeri", "Neosynchiropus ocellatus") ~ "Actinopteri",
      label %in% c("Bubalus carabanensis", "Canis lupus dingo", "Diceros bicornis minor", "Elephas maximus indicus", 
                   "Hippopotamus amphibius kiboko", "Lagenorhynchus acutus", "Ovis ammon polii", 
                   "Perognathus longimembris pacificus", "Gorilla gorilla gorilla", "Glossophaga mutica", "Molossus nigricans",
                   "Rangifer tarandus platyrhynchus") ~ "Mammalia",
      label %in% c("Elgaria multicarinata webbii", "Coluber constrictor foxii", "Spondylurus nitidus") ~ "Lepidosauria",
      label %in% c("Melospiza melodia melodia", "Motacilla alba alba", "Pithys albifrons albifrons") ~ "Aves",
      label %in% c("Mixophyes fleayi") ~ "Amphibia",
      TRUE ~ class  # Keep existing values if no match
    ),
    order = case_when(
      label %in% c("Chrysemys picta bellii", "Malaclemys terrapin pileata", "Macrochelys suwanniensis", "Mauremys mutica") ~ "Testudines",
      label %in% c("Latimeria chalumnae") ~ "Coelacanthiformes",
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
      TRUE ~ class
    )
  )

classifications <- plotdata %>% select(label, group)
groups <- unique(plotdata$group)
colours <- setNames(brewer.pal(n = 11, name = "Paired"), groups)

p <- ggtree(tree) + 
  geom_tiplab(size=1) 



p <- p + geom_facet(
  data = plotdata,
  mapping = aes(x = `Unit length`/1000, fill = group),
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
  theme_tree2() + xlim_expand(c(0, 40), "Tree")+scale_fill_manual(values = colours) + labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(3, 1, 1, 1, 1))

ggsave("zc4688/honours/ribocop/results/figures/tree_current.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/tree_current.pdf", height = 20, width = 15)



########################### Figure 2: Barrnap with edits (marked Palign for every unit)
structure_files <- list.files(directory, pattern = "\\.structure.tsv$", full.names = TRUE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each file and read the JSON into a data frame
for (file in structure_files) {
  data <- read.delim(file, header = TRUE, sep="\t")
  samplename <- str_replace(basename(file), ".refmorph.structure.tsv", "")
  data <- data %>% 
    select(Type, HMMbegin, HMMend, Envstart, Envend, `X..HMM`, `Target.length`, Strand) %>% 
    mutate(`Sample ID` = samplename) 
  data_list[[file]] <- data
}

barrnap <- bind_rows(data_list)

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

barrnap <- barrnap %>% left_join(refmorphs) 

barrnap <- barrnap %>% 
  mutate(morphstart = as.numeric(str_replace_all(str_extract(refmorph, "(:\\d*)"), ":", "")),
         morphend = as.numeric(str_replace_all(str_extract(refmorph, "(-\\d*)"), "-", "")))

barrnap <- barrnap %>%  left_join(lengths) 


barrnap <- barrnap %>% 
  mutate(Start = ifelse(Strand == "+", Envstart - morphstart, morphend - Envend),
         End = ifelse(Strand == "+", Envend - morphstart, morphend - Envstart) )



tree <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (5).nwk") #ORIGINAL tree before filtering, just for checking. most current results are in refiltered directory (with current filtering), all has more because filtering was lest strict

treenames <- as.data.frame(tree$tip.label)
colnames(treenames) <- c("Timetreename")
replacementnames <- read.delim("zc4688/honours/ribocop/speciestimetree.tsv", header = T, sep = ",")

replacementnames <- replacementnames %>% mutate(Species = str_replace_all(Species, " ", "_"), Timetreename = str_replace_all(Replacement, " ", "_") ) %>% select(-Replacement)

treenames <- treenames %>% left_join(replacementnames) %>% mutate(Species = ifelse(is.na(Species), Timetreename, Species)) %>% select(Species)
tree$tip.label <- as.character(treenames$Species)




barrnap <- barrnap %>% group_by(Species) %>% arrange(`Sample ID`) %>% 
  filter(`Sample ID` == last(`Sample ID`)) %>%
  ungroup()

barrnap <- barrnap %>% 
  dplyr::rename("label" = "Species") %>% 
  dplyr::rename("Hmmpalign" = "X..HMM")

barrnapplot <- barrnap %>% relocate(label, .before = everything())
################################################################################

msadf <- barrnapplot %>% 
  group_by(`Sample ID`, Type) %>% 
  filter((End - Start) == max(End - Start)) %>% 
  arrange(`Sample ID`, Type, Start) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(`Sample ID`, Type, Start, End, refmorph) %>% 
  mutate(length = End - Start, 
         name = paste(`Sample ID`, refmorph, sep = "_"))



eighteen <- msadf %>% filter(Type == "18S_rRNA")
eighteen <- eighteen %>% 
  mutate(outside = length > (mean(eighteen$length) + 3*sd(eighteen$length)) | 
           length < (mean(eighteen$length) - 3*sd(eighteen$length)))  %>% 
  filter(outside == FALSE) %>% 
  select(name, Start, End)
write.table(eighteen, "zc4688/honours/ribocop/results/new/msa/eighteen.tsv", sep = ",", quote = F, row.names = F, col.names = F)
twoeight <- msadf %>% filter(Type == "28S_rRNA")
twoeight <- twoeight %>% 
  mutate(outside = length > (mean(twoeight$length) + 3*sd(twoeight$length)) | 
           length < (mean(twoeight$length) - 3*sd(twoeight$length)))%>% 
  filter(outside == FALSE)
write.table(twoeight, "zc4688/honours/ribocop/results/new/msa/twoeight.tsv", sep = ",", quote = F, row.names = F, col.names = F)


fiveeight <- msadf %>% filter(Type == "5_8S_rRNA")
fiveeight <- fiveeight %>% 
  mutate(outside = length > (mean(fiveeight$length) + 3*sd(fiveeight$length)) | 
           length < (mean(fiveeight$length) - 3*sd(fiveeight$length)))%>% 
  filter(outside == FALSE)
write.table(fiveeight, "zc4688/honours/ribocop/results/new/msa/fiveeight.tsv", sep = ",", quote = F, row.names = F, col.names = F)


barrnapmetadata <- barrnapplot %>% 
  group_by(Type, label, `Sample ID`) %>% 
  mutate(Hmmpalign = ifelse(Hmmpalign > 1, 1, Hmmpalign)) %>% 
  mutate(Hmmpalign = round(Hmmpalign, digits = 2)) %>% 
  filter(Hmmpalign == max((Hmmpalign))) %>% 
  arrange(Unit, Envstart) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(label, `Sample ID`, Unit, Envstart, Envend, `Unit length`) %>% pivot_wider(
    names_from = Unit,
    values_from = c(Envstart, Envend),
    names_glue = "{Unit}_{.value}"
  ) %>% 
  mutate(ITS1 = `5_8S_rRNA_Envstart` - `18S_rRNA_Envend`,
         `ITS2` = `28S_rRNA_Envstart` - `5_8S_rRNA_Envend`,
         `18S_length` = `18S_rRNA_Envend` - `18S_rRNA_Envstart`,
         `28S_length` = `28S_rRNA_Envend` - `28S_rRNA_Envstart`,
         `IGS_length` = `Unit length` - `28S_rRNA_Envend`) %>% 
  select(ITS1, ITS2, `18S_length`, `28S_length`, `IGS_length`, label, `Unit length`, `Sample ID`)

################################################################################

common_organisms <- intersect(barrnapplot$label, tree$tip.label)


barrnapplot <- barrnapplot %>%
  filter(label %in% common_organisms)





barrnapplot <- barrnapplot %>% mutate(label = str_replace_all(label, "_", " "))
barrnapplot <- left_join(barrnapplot, classifications, by="label")

#pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$edge.length <- NULL
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")




barrnapplot <- barrnapplot %>% mutate(Hmmpalign = Hmmpalign * 100)




tmp <- barrnapplot %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + geom_facet(
  data = barrnapplot,
  mapping = aes(x = 0, xend = ifelse(`Unit length` < 35000, `Unit length`, 35000)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey50",
  alpha = 0.3
) +
  new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Start, xend = End, color = Type, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
)  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 40), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

p

ggsave("zc4688/honours/ribocop/results/figures/barrnaptree.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/barrnaptree.pdf", height = 20, width = 15)

tmp <- barrnapplot %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 40), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

p

ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_nolength.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_nolength.pdf", height = 20, width = 15)


p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

test <- barrnapplot %>% select(label, `Unit length`, group) %>% rename("testgroup" = "group")
p <- p  + new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = "Unit Structure"
)  + geom_facet(
  data = test,
  mapping = aes(x = `Unit length`/1000, fill = testgroup),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 40), "Tree")+labs(alpha = "% HMM Match") + scale_fill_manual(values = colours) + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1),
    fill = guide_none()) + coord_cartesian(clip = 'off')

facet_widths(p, widths = c(3, 3, 1))


ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_lengths.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_lengths.pdf", height = 20, width = 15)


############################### Per group plots
mammals <- barrnapplot %>% filter(group == "Mammalia")

mammaltree <- ape::drop.tip(tree, setdiff(tree$tip.label, mammals$label))

p <- ggtree(mammaltree)  +
  geom_tiplab(size=2) 

p <- p + geom_facet(
  data = mammals %>% filter(Envend < 15000),
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = "Unit Structure"
)  +  geom_facet(
  data = plotdata %>% filter(group == "Mammalia"),
  mapping = aes(x = `Unit length`/1000),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 , fill = "#B2DF8A"
) +
  scale_y_discrete() +  
  theme_tree2() +xlim_expand(c(0, 28), "Tree") + labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

facet_widths(p, widths = c(3, 3, 1))

ggsave("zc4688/honours/ribocop/results/figures/mammaltree.pdf", height = 10, width = 8)




#inspect the merging
p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")
p <- p + new_scale_color() + geom_facet(
  data = test,
  mapping = aes(x = Envstart, xend = Envend, color = type),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
)  + 
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')
p
ggsave("zc4688/honours/ribocop/results/figures/gaptree.png", height = 20, width = 15)



############################### Boxplots

barnapboxplot <- barrnapplot %>% 
  group_by(`Sample ID`, Type) %>% 
  filter((End - Start) == max(Envend - Envstart)) %>% 
  arrange(`Sample ID`, Type, Start) %>% 
  slice_head(n = 1) %>% 
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



barnapboxplot <- barnapboxplot %>% 
  group_by(group) %>% mutate(n = n()) %>% ungroup() %>%  filter(n >=5) %>% select(-n) %>% 
  pivot_longer(cols = c(`ITS1`, ITS2, `18S_length`, `28S_length`, `IGS_length`, `Unit length`), names_to = "Variable", values_to = "Value") %>% 
  dplyr::rename(label = group) %>% 
  select(label, Variable, Value)



ggplot(barnapboxplot, aes(y = Value, fill = label)) +
  geom_boxplot(coef = 2.5, outliers = FALSE, notch = TRUE) +
  facet_wrap(~ Variable, scales = "free_y")  + theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = colours) +
  labs(fill = "Taxonomic Group")

ggsave("zc4688/honours/ribocop/results/figures/boxplot.png", height = 8, width = 12)


############################### Boxplots with tree
classtree <- read.tree("/g/data/te53/zc4688/honours/ribocop/phyliptree (2).phy")


classtree$tip.label <- c("Mammalia", "Aves", "Lepidosauria", "Testudines and Crocodylia", "Amphibia", "Actinopterygii", "Chondrichthyes")

common_class <- intersect(barnapboxplot$label, classtree$tip.label)


barnapboxplot <- barnapboxplot %>%
  filter(label %in% common_class)



classtree$edge.length <- NULL
classtree <- ape::drop.tip(classtree, setdiff(classtree$tip.label, common_class))

############################### ITS only
p <- ggtree(classtree) + 
  geom_tiplab(size=3) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours) + theme(legend.position = "none")



p <- p + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS1"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS1 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS2"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS2 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position = "none") + labs(color = NULL, fill = NULL) + scale_fill_brewer(palette = "Set1") 

facet_widths(p, widths = c(1, 1, 1))
ggsave("zc4688/honours/ribocop/results/figures/itsboxplot.pdf", height = 8, width = 12)


############################### All lengths

p <- ggtree(classtree) + 
  geom_tiplab(size=3, offset = 0.5) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours) + theme(legend.position = "none")


p <- p + new_scale_color() + geom_facet(
  data = barnapboxplot %>% filter(Variable == "Unit length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "Unit Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS1"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS1 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS2"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS2 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +  geom_facet(
  data = barnapboxplot %>% filter(Variable == "IGS_length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "IGS Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +
  geom_facet(
    data = barnapboxplot %>% filter(Variable == "18S_length"),
    mapping = aes(x = Value/1000, fill = Variable, group = interaction(label, Variable)),
    geom = geom_boxplot,       
    panel = "18S Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
    data = barnapboxplot %>% filter(Variable == "28S_length"),
    mapping = aes(x = Value/1000, fill = Variable, group = interaction(label, Variable)),
    geom = geom_boxplot,       
    panel = "28S Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position =  "none")  + scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  coord_cartesian(clip = "off")

facet_widths(p, widths = c(1, 1, 1))
ggsave("zc4688/honours/ribocop/results/figures/combinedboxplots.pdf", height = 8, width = 12)

############################### Summary stats
barnapboxplot %>% 
  group_by(label, Variable) %>% 
  summarise(Median = median(Value, na.rm = TRUE), Min = min(Value, na.rm = TRUE),  Max = max(Value, na.rm = TRUE), Species = n()) %>% 
  view



############################### Metadata table

metadata <- metadata %>% 
  left_join(barrnapmetadata)

metadata <- metadata %>% 
  rename("rDNA contigs" = "Contigs", "Species" = "label")

write.table(metadata, "zc4688/honours/ribocop/metadata.tsv", sep = "\t", quote = F, row.names = F)

################################Conservation

eighteen <- readDNAMultipleAlignment("zc4688/honours/ribocop/results/main/morph_fasta/msa/twoeight_filtered_msa.txt")

shannon_entropy <- function(column) {
  column <- column[column != "-"]
  freqs <- table(column) / length(column)  # Frequency of each base
  -sum(freqs * log2(freqs), na.rm = TRUE)  # Shannon entropy formula
}

# Convert alignment to matrix
msa_matrix <- as.matrix(eighteen)

# Apply function to each column
entropy_values <- apply(msa_matrix, 2, shannon_entropy)

human <- DNAStringSet(eighteen)["GCF_000001405_NC_000021.9:8392567-8436777:+:99-1967"]
human <- strsplit(as.character(human), "")[[1]]

tmp <- data.frame(Position = 1:length(entropy_values), Score = entropy_values, 
                  Sequence = human)

tmp <- tmp %>%
  arrange(Position) %>%
  mutate(Entropy_Smoothed = zoo::rollmean(Score, k = 10, fill = NA, align = "center"))

ggplot(tmp %>% filter(Sequence != "-")) + geom_point(aes(x = Position, y = Score, color = Score), size = 0.5) + 
  scale_color_gradient(low = "red", high = "green") + new_scale_colour() + 
  geom_rug(aes(x = Position, y = 0, color = Sequence), size = 2, sides = "b") +
  scale_color_manual(values = c("A" = "blue", "C" = "purple", "G" = "purple", "T" = "blue"))


############################### 18S tree
eighteentree <- read.tree("/g/data/te53/zc4688/honours/ribocop/results/main/morph_fasta/msa/twoeight_msa.txt.treefile")

treenames <- as.data.frame(eighteentree$tip.label)

treenames <- treenames %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", `eighteentree$tip.label`))

treenames <- treenames %>% left_join(plotdata %>% select(`Sample ID`, label))

eighteentree$tip.label <- as.character(treenames$label)
common_organisms <- intersect(plotdata$label, eighteentree$tip.label)

eighteentree <- ape::drop.tip(eighteentree, setdiff(eighteentree$tip.label, common_organisms))
eighteentree$edge.length <- NULL
eighteentree$tip.label <- str_replace_all(tree$tip.label, "_", " ")



ggtree(eighteentree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

################################
tree <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (5).nwk") 

treenames <- as.data.frame(tree$tip.label)
colnames(treenames) <- c("Timetreename")
replacementnames <- read.delim("zc4688/honours/ribocop/speciestimetree.tsv", header = T, sep = ",")

replacementnames <- replacementnames %>% mutate(Species = str_replace_all(Species, " ", "_"), Timetreename = str_replace_all(Replacement, " ", "_") ) %>% select(-Replacement)

treenames <- treenames %>% left_join(replacementnames) %>% mutate(Species = ifelse(is.na(Species), Timetreename, Species)) %>% select(Species)
tree$tip.label <- as.character(treenames$Species)
tree$edge.length <- NULL
tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")

eighteentree <- read.tree("/g/data/te53/zc4688/honours/ribocop/results/main/morph_fasta/msa/twoeight_filtered_msa.txt.bionj")

treenames <- as.data.frame(eighteentree$tip.label)


treenames <- treenames %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", `eighteentree$tip.label`))

treenames <- treenames %>% left_join(plotdata %>% select(`Sample ID`, label))

eighteentree$tip.label <- as.character(treenames$label)
eighteentree$edge.length <- NULL

eighteentree$tip.label <- str_replace_all(tree$tip.label, "_", " ")


common_organisms <- intersect(eighteentree$tip.label, tree$tip.label)

eighteentree <- ape::drop.tip(eighteentree, setdiff(eighteentree$tip.label, common_organisms))
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))
A <- cbind(tree$tip.label, tree$tip.label)

x$color <- colours[x$group]
dist.topo(tree, eighteentree)


 #cophyloplot(tree, eighteentree, assoc = A, show.tip.label = F, space=500, col = x$color, type = "cladogram") 



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

