library(ggtree)
library(tidyverse)
library(ape)
library(jsonlite)
library(GenomicRanges)
library(taxize)
library(RColorBrewer)
############################### Figure 1: length tree
directory <- "/g/data/te53/zc4688/honours/ribocop/results/chordata_filter2"

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

combined_df <- combined_df %>% 
  filter(is.na(`Errors.rDNA identification`) & is.na(`Errors.Morph identification`)) %>% 
  filter(!is.na(`rDNA_details.Number of primary alignments`)) %>% 
  dplyr::select(`Sample_details.Sample ID`, `Sample_details.Species`, `Sample_details.TaxID`, 
                `rDNA_details.Unit length`, `rDNA_details.Total length`, `rDNA_details.Number of contigs`,
                `rDNA_details.Minimum unit length`, `rDNA_details.Maximum unit length`) %>% 
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
  rename("label" = "Species")

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
                   "Pempheris klunzingeri") ~ "Actinopteri",
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
      label %in% c("Chrysemys picta bellii", "Malaclemys terrapin pileata", "Macrochelys suwanniensis") ~ "Testudines",
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
barrnap <- read.delim("zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/barrnap_newpalign.gff", header = F)
tree <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (5).nwk") #ORIGINAL tree before filtering, just for checking. most current results are in refiltered directory (with current filtering), all has more because filtering was lest strict

treenames <- as.data.frame(tree$tip.label)
colnames(treenames) <- c("Timetreename")
replacementnames <- read.delim("zc4688/honours/ribocop/speciestimetree.tsv", header = T, sep = ",")

replacementnames <- replacementnames %>% mutate(Species = str_replace_all(Species, " ", "_"), Timetreename = str_replace_all(Replacement, " ", "_") ) %>% select(-Replacement)

treenames <- treenames %>% left_join(replacementnames) %>% mutate(Species = ifelse(is.na(Species), Timetreename, Species)) %>% select(Species)
tree$tip.label <- as.character(treenames$Species)

colnames(barrnap) <- c("Sample", "Version", "Start", "End", "Envstart", "Envend", "Hmmstart", "Hmmend", "Score", "Bitscore", "Hmmpalign", "Targetpalign", "Strand", "Other", "Product")

barrnap <- barrnap %>% 
  filter(!is.na(Start)) %>% 
  separate(Product, into = c("Name", "Other", "Note"), sep = ";") %>% 
  separate(Name, into = c("Name", "Unit"), sep = "=") %>% 
  #separate(Other, into = c("Other", "Completeness"), sep = "\\(") %>% 
  # separate(Note, into = c("Note", "Palign"), sep = "=") %>% 
  #mutate(Completeness = str_replace_all(Completeness, "\\)", "")) %>% 
  select(-Note, -Name, -Other, -Version, -Strand) %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", Sample)) %>% 
  left_join(lengths) 


barrnap <- barrnap %>% group_by(Species) %>%
  filter(`Sample ID` == last(`Sample ID`)) %>%
  ungroup()

barrnap <- barrnap %>% 
  rename("label" = "Species")

barrnapplot <- barrnap %>% relocate(label, .before = everything())

barrnapmetadata <- barrnapplot %>% 
  group_by(Unit, label, `Sample ID`) %>% 
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

p <- p + new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
)  +  geom_facet(
  data = barrnapplot %>% filter(`Unit length` < 35000),
  mapping = aes(x = `Unit length`, xend = `Unit length` + 50),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "black"
) +  geom_facet(
  data = barrnapplot %>% filter(`Unit length` >= 35000),
  mapping = aes(x = 35000, xend = 35050),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "black"
)+
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


############################### Separate vs broken alignments
test <- barrnapplot %>% 
  arrange(Sample, Envstart) %>% 
  group_by(label, Sample, Unit) %>% 
  mutate(query_gap = Envstart - lag(Envend, default = first(Envstart)), 
         ref_gap = Hmmstart - lag(Hmmend, default = first(Hmmstart)), 
         diff = query_gap - ref_gap, 
         type = case_when((ref_gap < -100 & query_gap > 100) | query_gap > 5000 ~ "separate",
                          ref_gap == 0 & query_gap == 0 ~ "single", 
                          TRUE ~ "gap")) %>% 
  select(label, Sample, Unit, Envstart, Envend, type, `Unit length`, group) %>% 
  mutate(
    merged_Envstart = if_else(type == "gap", lag(Envstart, default = first(Envstart)), Envstart),
    merged_Envend = Envend)  %>% 
  filter(!lead(type, default = "none") == "gap") %>%
  ungroup() %>% 
  filter(Unit != "5S_rRNA") %>% 
  mutate(Envstart = merged_Envstart,
         Envend = merged_Envend)

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

barnapboxplot <- test %>% 
  group_by(Sample, Unit) %>% 
  filter((Envstart - Envend) == max(Envstart - Envend)) %>% 
  arrange(Sample, Unit, Envstart) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(Sample, Unit, Envstart, Envend, group, `Unit length`) %>% pivot_wider(
    names_from = Unit,
    values_from = c(Envstart, Envend),
    names_glue = "{Unit}_{.value}"
  ) %>% 
  mutate(ITS1 = `5_8S_rRNA_Envstart` - `18S_rRNA_Envend`,
         `ITS2` = `28S_rRNA_Envstart` - `5_8S_rRNA_Envend`,
         `18S_length` = `18S_rRNA_Envend` - `18S_rRNA_Envstart`,
         `28S_length` = `28S_rRNA_Envend` - `28S_rRNA_Envstart`,
         `IGS_length` = `Unit length` - `28S_rRNA_Envend`)



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
