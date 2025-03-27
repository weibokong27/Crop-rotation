


#############################################################################################################################################

library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(ape)

################################################################# 16S ######################################################################
load("DESeq2.CropSystem.microbiomeMarker_DADA2.rare10000.RData")
DESeq2.data = marker_table(deseq2.16s) %>% data.frame()
colnames(DESeq2.data)

load("DESeq2.CropSystem_DADA2.PREV0.01.rare5000.RData")
DESeq2.data = subset(deseq2.16s, padj < 0.05) %>% rownames_to_column(var = "OTUID") %>%
  dplyr::select(feature = "OTUID", group = "group", ef_logFC = "log2FoldChange", pvalue = "pvalue", padj = "padj") %>% 
  mutate(enrich_group = ifelse(group == "up-regulated", "RS", "CS"))


others <- c("p__Armatimonadetes", "p__Latescibacteria", "p__Nitrospirae", "p__Rokubacteria", "p__WPS-2")
phylum.color.16s <- tribble(~Phylum, ~color, 
                        "c__Alphaproteobacteria","#80cdc1",
                        "c__Gammaproteobacteria","#009e73",
                        "c__Deltaproteobacteria","#114662",
                        "p__Actinobacteria", "#B2DF8A", 
                        "p__Acidobacteria", "#1F78B4", 
                        "p__Bacteroidetes", "#FB9A99",
                        "p__Chloroflexi", "#E31A1C",
                        "p__Firmicutes","#6A3D9A",
                        "p__Gemmatimonadetes","#FDBF6F",
                        #"p__Nitrospirae","#075c89",
                        "p__Planctomycetes","#CAB2D6",
                        "p__Cyanobacteria","#A6CEE3",
                        "p__Verrucomicrobia","#75b831",
                        "p__Patescibacteria", "#b04380",
                        "Others", "#a6a6a6") %>% deframe()

OTUID = DESeq2.data$feature

node = which(!(phy_tree(ps.16s.norm)$tip.label %in% OTUID))
(tree = drop.tip(as.phylo(phy_tree(ps.16s.norm)), node))# 

tax_table = data.frame(OTUID = OTUID, tax_table(ps.16s.norm)[OTUID, ])# %>% dplyr::select(., -label)
tax_table = tax_table %>% mutate(., label = .$Phylum) %>% mutate(., label = ifelse(.$label == "p__Proteobacteria", .$Class, .$Phylum))
tax_table[tax_table$label %in% others, 'label'] <- "Others" 

groupInfo = split(tree$tip.label, as.vector(tax_table[tree$tip.label, "label"]))
tree = tidytree::groupOTU(tree, groupInfo)


logFC = data.frame(OTUID = DESeq2.data$feature,
                   Enrich = DESeq2.data$enrich_group, 
                   logFC = DESeq2.data$ef_logFC) %>% 
  mutate(logFC2 = abs(logFC)) %>% 
  mutate(logFC2 = scales::rescale(logFC2, c(0.2, 1))) %>% 
  mutate(logFC2 = ifelse(logFC > 0, logFC2, -logFC2)) %>%
  dplyr::left_join(tax_table %>% dplyr::select(OTUID, Phylum = label), by = "OTUID")

abundance = prune_taxa(OTUID, ps.16s.norm) %>% speedyseq::psmelt() %>% 
  dplyr::group_by(OTU, CropSystem) %>% 
  dplyr::summarise(Abundance = mean(Abundance)) %>% mutate(., Abundance = log1p(.$Abundance))

ggtree(tree, layout = "circular", size = 0.3, aes(color = group)) + #%<+% tax_table + 
  xlim(-0.5, NA) + 
  scale_fill_manual(values = phylum.color.16s) + scale_color_manual(values = phylum.color.16s) + 
  new_scale_fill() + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1), offset = -0.15, fill = "#ddf3ff", width = 0.35, alpha = 1) + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1), offset = 0.12, fill = "#fbf9d9", width = 0.35, alpha = 1) + 
  new_scale_color() + 
  geom_fruit(data = logFC, geom = geom_bar, aes(y = OTUID, x = logFC2, fill = Phylum), offset = -0.16, pwidth = 0.35, width = 0.9, orientation = "y", stat = "identity") + 
  scale_fill_manual(values = phylum.color.16s) + 
  new_scale_fill() + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1, fill = Enrich), offset = -0.15, width = 0.12, color = NA) + 
  scale_fill_manual(values = Palette.CropSystem) + 
  new_scale_fill() + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1, fill = Phylum), offset = -0.11, width = 0.02, color = NA) +
  scale_fill_manual(values = phylum.color.16s)


##############################################################################################################################################

abun.enrich = ps.16s.norm %>% 
  transform_sample_counts(., function(x) 100 * x / sum(x)) %>% 
  prune_taxa(OTUID, .) %>% 
  speedyseq::psmelt() %>% 
  mutate(., label = .$Phylum) %>% 
  mutate(., label = ifelse(.$label == "p__Proteobacteria", .$Class, .$Phylum)) %>% 
  dplyr::group_by(label) %>% 
  dplyr::summarise(RelAbun = mean(Abundance)) %>% 
  rbind(., data.frame(label = 'Overall', RelAbun = sum(.$RelAbun)))

num.enrich = ps.16s.norm %>% 
  transform_sample_counts(., function(x) 100 * x / sum(x)) %>% 
  prune_taxa(OTUID, .) %>% tax_table() %>% data.frame() %>% 
  mutate(., label = .$Phylum) %>% 
  mutate(., label = ifelse(.$label == "p__Proteobacteria", .$Class, .$Phylum)) %>% 
  dplyr::select(label) %>% rownames_to_column(var = "OTUID") %>%
  left_join(logFC %>% dplyr::select(OTUID, Enrich), by = "OTUID") %>% 
  dplyr::group_by(label, Enrich) %>% 
  dplyr::summarise(num = n()) %>% 
  dplyr::ungroup()

num.enrich %<>% dplyr::group_by(Enrich) %>% dplyr::summarize(num = sum(num)) %>% mutate(label = 'Overall') %>% rbind(num.enrich) 

num.all.enrich = num.enrich %>% dplyr::group_by(label) %>% dplyr::summarize(all = sum(num))

num.enrich %>% left_join(num.all.enrich, by = "label") %>% 
  mutate(prop = num / all) %>% 
  left_join(., abun.enrich, by = "label") %>%
  subset(!label %in% "Overall") %>%
  ggplot(., aes(x = RelAbun/2, y = prop, fill = Enrich, width = RelAbun)) + 
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = Palette.CropSystem) +
  coord_polar("y", start = 0) +
  facet_wrap(~ label, ncol = 5) +
  # theme_void() +
  theme(legend.position = "bottom")


################################################################# ITS ######################################################################
load("DESeq2.CropSystem.microbiomeMarker_DADA2.PREV0.01.rare5000.RData")
DESeq2.data = marker_table(deseq2.its) %>% data.frame()

load("DESeq2.CropSystem_DADA2.PREV0.01.rare5000.RData")
DESeq2.data = subset(deseq2.its, padj < 0.05) %>% rownames_to_column(var = "OTUID") %>%
  dplyr::select(feature = "OTUID", group = "group", ef_logFC = "log2FoldChange", pvalue = "pvalue", padj = "padj") %>% 
  mutate(enrich_group = ifelse(group == "up-regulated", "RS", "CS"))


others <- c("p__Mucoromycota", "p__unidentified", "p__Zoopagomycota", "p__Blastocladiomycota", "p__Aphelidiomycota", 'p__Rozellomycota')
phylum.color.its <- tribble(~Phylum, ~color, 
                        "c__Sordariomycetes","#80cdc1",
                        "c__Dothideomycetes","#1F78B4",
                        "c__Eurotiomycetes","#B2DF8A",
                        "c__Leotiomycetes", "#FDBF6F",
                        "c__Pezizomycetes", "#FB9A99",
                        "p__Basidiomycota","#E31A1C",
                        "p__Chytridiomycota","#7570B3",
                        "p__Glomeromycota","#D95F02",
                        "p__Mortierellomycota", "#A6761D", 
                        #"p__Rozellomycota", "#E6AB02", 
                        "Others", "#a6a6a6") %>% deframe()

Palette.16s

OTUID = DESeq2.data$feature

node = which(!(phy_tree(ps.its.norm)$tip.label %in% OTUID))
(tree = drop.tip(as.phylo(phy_tree(ps.its.norm)), node))# 

tax_table = data.frame(OTUID = OTUID, tax_table(ps.its.norm)[OTUID, ])# %>% dplyr::select(., -label)
tax_table = tax_table %>% mutate(., label = .$Phylum) %>% mutate(., label = ifelse(.$label == "p__Ascomycota", .$Class, .$Phylum))
tax_table[tax_table$label %in% others, 'label'] <- "Others" 


groupInfo = split(tree$tip.label, as.vector(tax_table[tree$tip.label, "label"]))
tree = tidytree::groupOTU(tree, groupInfo)


logFC = data.frame(OTUID = DESeq2.data$feature,
                   Enrich = DESeq2.data$enrich_group, 
                   logFC = DESeq2.data$ef_logFC) %>% 
  mutate(logFC2 = abs(logFC)) %>% 
  mutate(logFC2 = scales::rescale(logFC2, c(0.2, 1))) %>% 
  mutate(logFC2 = ifelse(logFC > 0, logFC2, -logFC2)) %>% 
  dplyr::left_join(tax_table %>% dplyr::select(OTUID, Phylum = label), by = "OTUID")

abundance = prune_taxa(OTUID, ps.its.norm) %>% speedyseq::psmelt() %>% 
  dplyr::group_by(OTU, CropSystem) %>% 
  dplyr::summarise(Abundance = mean(Abundance)) %>% mutate(., Abundance = log1p(.$Abundance))


ggtree(tree, layout = "circular", size = 0.3, aes(color = group), branch.length = "none") + #%<+% tax_table + 

  scale_fill_manual(values = phylum.color.its) + scale_color_manual(values = phylum.color.its) + 
  new_scale_fill() + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1), offset = -0.07, fill = "#ddf3ff", width = 6, alpha = 0.7) + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1), offset = 0.03, fill = "#fbf9d9", width = 6, alpha = 0.7) + 
  geom_fruit(data = logFC, geom = geom_bar, aes(y = OTUID, x = logFC2, fill = Phylum), offset = -0.11, pwidth = 0.25, width = 0.9, orientation = "y", stat = "identity") + 
  scale_fill_manual(values = phylum.color.its) +   
  new_scale_fill() + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1, fill = Enrich), offset = -0.15, width = 2, color = NA) + 
  scale_fill_manual(values = Palette.CropSystem) + 
  new_scale_fill() + 
  geom_fruit(data = logFC, geom = geom_tile, aes(y = OTUID, x = 1, fill = Phylum), offset = -0.12, width = 0.5, color = NA) +
  scale_fill_manual(values = phylum.color.its)


abun.enrich = ps.its.norm %>% 
  transform_sample_counts(., function(x) 100 * x / sum(x)) %>% 
  prune_taxa(OTUID, .) %>% 
  speedyseq::psmelt() %>% 
  mutate(., label = .$Phylum) %>% 
  mutate(., label = ifelse(.$label == "p__Ascomycota", .$Class, .$Phylum)) %>% 
  dplyr::group_by(label) %>% 
  dplyr::summarise(RelAbun = mean(Abundance)) %>% 
  rbind(., data.frame(label = 'Overall', RelAbun = sum(.$RelAbun)))

num.enrich = ps.its.norm %>% 
  transform_sample_counts(., function(x) 100 * x / sum(x)) %>% 
  prune_taxa(OTUID, .) %>% tax_table() %>% data.frame() %>% 
  mutate(., label = .$Phylum) %>% 
  mutate(., label = ifelse(.$label == "p__Ascomycota", .$Class, .$Phylum)) %>% 
  dplyr::select(label) %>% rownames_to_column(var = "OTUID") %>%
  left_join(logFC %>% dplyr::select(OTUID, Enrich), by = "OTUID") %>% 
  dplyr::group_by(label, Enrich) %>% 
  dplyr::summarise(num = n()) %>% 
  dplyr::ungroup()

num.enrich %<>% dplyr::group_by(Enrich) %>% dplyr::summarize(num = sum(num)) %>% mutate(label = 'Overall') %>% rbind(num.enrich) 

num.all.enrich = num.enrich %>% dplyr::group_by(label) %>% dplyr::summarize(all = sum(num))

num.enrich %>% left_join(num.all.enrich, by = "label") %>% 
  mutate(prop = num / all) %>% 
  left_join(., abun.enrich, by = "label") %>%
  subset(!label %in% "Overall") %>%
  ggplot(., aes(x = RelAbun/2, y = prop, fill = Enrich, width = RelAbun)) + 
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = Palette.CropSystem) +
  coord_polar("y", start = 0) +
  facet_wrap(~ label, ncol = 5) +
  # theme_void() +
  theme(legend.position = "bottom")





