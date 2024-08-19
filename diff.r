DefaultAssay(adhesion_MAC_MONO.sct) <- 'RNA'

Idents(adhesion_MAC_MONO.sct) <- 'celltype'
adhesion_MAC_MONO.sct.TREM2 <- subset(adhesion_MAC_MONO.sct,
                                            idents = "M3")
Idents(adhesion_MAC_MONO.sct.TREM2) <- 'group'
TREM2_diff_IUA <- adhesion_MAC_MONO.sct.TREM2@misc$TREM2_diff_IUA <- 
  FindMarkers(adhesion_MAC_MONO.sct.TREM2,
              ident.1 = c('IUA'), 
              ident.2 = c('Normal'), 
              verbose = FALSE,
              logfc.threshold = 0)


TREM2_diff_IUA <- TREM2_diff_IUA %>% 
  mutate(sig=case_when(avg_log2FC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_log2FC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_log2FC > -0.25 | avg_log2FC < 0.25 ~ 'non')) %>%
  mutate(celltype='M3 MAC') %>% 
  mutate(gene=row.names(TREM2_diff_IUA)) %>% 
  rownames_to_column("gene1")


TREM2_diff_IUA.gene.ego <- TREM2_diff_IUA  %>%
  dplyr::filter(sig != 'non')  %>%
  dplyr::pull(gene, name=avg_log2FC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


phagocytosis,cytokine-mediated signaling pathway,

cell chemotaxis, regulation of autophagy,reactive oxygen species metabolic process
chemokine 


phagocytosis <- list(unique(phagocytosis$V1))

PHAGO <- AddModuleScore(
  object = adhesion_MAC_MONO.sct.TREM2,
  features = phagocytosis,
  ctrl = 100, #默认值是100
  name = 'PHAGO'
)
colnames(PHAGO@meta.data)[17] <- 'phago' 
data<- FetchData(PHAGO, vars = c("celltype","group","phago"))

data$group <- factor(data$group, 
                     levels=c("Normal",
                              "IUA"))

ggplot(data , aes(x=group,y=phago), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = phago,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  )    



phagocytosis <- list(unique(phagocytosis$V1))

PHAGO <- AddModuleScore(
  object = adhesion_MAC_MONO.sct.TREM2,
  features = phagocytosis,
  ctrl = 100, #默认值是100
  name = 'PHAGO'
)
colnames(PHAGO@meta.data)[17] <- 'phago' 
data<- FetchData(PHAGO, vars = c("celltype","group","phago"))

data$group <- factor(data$group, 
                     levels=c("Normal",
                              "IUA"))

TREM2_phagocytosis <- ggplot(data , aes(x=group,y=phago), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = phago,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  )    



lipid_metabolic <- list(unique(lipid_metabolic$V1))

LIPID <- AddModuleScore(
  object = adhesion_MAC_MONO.sct.TREM2,
  features = lipid_metabolic,
  ctrl = 100, #默认值是100
  name = 'LIPID'
)
colnames(LIPID@meta.data)[17] <- 'lipid' 
data<- FetchData(LIPID, vars = c("celltype","group","lipid"))

data$group <- factor(data$group, 
                     levels=c("Normal",
                              "IUA"))

TREM2_lipid <- ggplot(data , aes(x=group,y=lipid), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = lipid,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  ) 

TREM2_phagocytosis + TREM2_lipid




#################


TREM2_diff_lipid <- TREM2_diff_IUA %>% dplyr::filter(sig != 'non') %>%
  dplyr::filter(gene %in% lipid_metabolic[[1]])



TREM2_diff_pha <- TREM2_diff_IUA %>% dplyr::filter(sig != 'non') %>%
  dplyr::filter(gene %in% phagocytosis[[1]])


FeatureStatPlot(
  srt = adhesion_MAC_MONO.sct.TREM2, group.by = "group", bg.by = "group",
  stat.by = c(TREM2_diff_lipid$gene1[1:10]), add_box = TRUE,
  comparisons = list(
    c("IUA", "Normal")
  )
)


FeatureStatPlot(adhesion_MAC_MONO.sct.TREM2, 
                stat.by = c(TREM2_diff_lipid$gene1[1:10]), 
                group.by = "group", plot.by = "group")

FeatureStatPlot(adhesion_MAC_MONO.sct.TREM2,  
                stat.by = c(TREM2_diff_lipid$gene1[1:5]),
                group.by = "group", 
                bg.by = "group", add_box = TRUE, stack = TRUE,
                comparisons = list(
                  c("IUA", "Normal")),
                flip = FALSE) + 


FeatureStatPlot(adhesion_MAC_MONO.sct.TREM2,  
                stat.by = c(TREM2_diff_pha$gene1[1:5]),
                group.by = "group", 
                bg.by = "group", add_box = TRUE, stack = TRUE,
                comparisons = list(
                  c("IUA", "Normal")),
                flip = FALSE)




##########  pttg1 epi

CellDimPlot(
  srt = adhesion_CILIA_UNCILIA.sct, group.by = c("celltype"),
  reduction = "UMAP", theme_use = "theme_blank"
)



DefaultAssay(adhesion_CILIA_UNCILIA.sct) <- 'RNA'

Idents(adhesion_CILIA_UNCILIA.sct) <- 'celltype'
adhesion_CILIA_UNCILIA.sct.PTTG1 <- subset(adhesion_CILIA_UNCILIA.sct,
                                      idents = "U4")
Idents(adhesion_CILIA_UNCILIA.sct.PTTG1) <- 'group'
PTTG1_diff_IUA <- adhesion_CILIA_UNCILIA.sct.PTTG1@misc$PTTG1_diff_IUA <- 
  FindMarkers(adhesion_CILIA_UNCILIA.sct.PTTG1,
              ident.1 = c('IUA'), 
              ident.2 = c('Normal'), 
              verbose = FALSE,
              logfc.threshold = 0)


PTTG1_diff_IUA <- PTTG1_diff_IUA %>% 
  mutate(sig=case_when(avg_log2FC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_log2FC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_log2FC > -0.25 | avg_log2FC < 0.25 ~ 'non')) %>%
  mutate(celltype='PTTG1 EPI') %>% 
  mutate(gene=row.names(PTTG1_diff_IUA)) %>% 
  rownames_to_column("gene1")

PTTG1_diff_IUA.gene.ego <- PTTG1_diff_IUA  %>%
  dplyr::filter(sig != 'non')  %>%
  dplyr::pull(gene, name=avg_log2FC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


#########  gene set scores  stem cell proliferation//response to oxidative stress//regeneration//
#########   epithelial to mesenchymal transition  stem cell fate commitment/// stem cell population maintenance
############# Wnt signaling pathway


#### diff genes
#自定义颜色：
mycol <- c("#3288BD","#d8d8d8","tomato")
#自定义主题：
mytheme <- theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(15,5.5,5.5,5.5))

p <- ggplot(data = PTTG1_diff_IUA,
       aes(x = avg_log2FC,
           y = -log10(p_val_adj),
           color = sig)) + 
  geom_point(size = 2.2) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + 
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-5, 5, by = 2)) + 
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 60),
                     breaks = seq(0, 60, by = 10))+ 
  geom_hline(yintercept = c(-log10(0.05)),size = 0.5,color = "black",lty = "dashed") + 
  geom_vline(xintercept = c(-0.25, 0.25),size = 0.5,color = "black",lty = "dashed") + 
  mytheme

up <- dplyr::filter(PTTG1_diff_IUA, sig == 'up') %>% dplyr::distinct(gene, .keep_all = T) %>%
  dplyr::top_n(3, avg_log2FC)
down <- dplyr::filter(PTTG1_diff_IUA, sig == 'down') %>% distinct(gene, .keep_all = T) %>%
  top_n(3, -avg_log2FC)

pttg1_vlo <- p + geom_text_repel(data = up,
                    aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_text_repel(data = down,
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene))




####################

DefaultAssay(adhesion_CILIA_UNCILIA.sct) <- 'RNA'

Idents(adhesion_CILIA_UNCILIA.sct) <- 'celltype'
adhesion_CILIA_UNCILIA.sct.MMP7 <- subset(adhesion_CILIA_UNCILIA.sct,
                                           idents = "U2")
Idents(adhesion_CILIA_UNCILIA.sct.MMP7) <- 'group'
MMP7_diff_IUA <- adhesion_CILIA_UNCILIA.sct.MMP7@misc$MMP7_diff_IUA <- 
  FindMarkers(adhesion_CILIA_UNCILIA.sct.MMP7,
              ident.1 = c('IUA'), 
              ident.2 = c('Normal'), 
              verbose = FALSE,
              logfc.threshold = 0)


MMP7_diff_IUA <- MMP7_diff_IUA %>% 
  mutate(sig=case_when(avg_log2FC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_log2FC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_log2FC > -0.25 | avg_log2FC < 0.25 ~ 'non')) %>%
  mutate(celltype='MMP7 EPI') %>% 
  mutate(gene=row.names(MMP7_diff_IUA)) %>% 
  rownames_to_column("gene1")


############

DefaultAssay(adhesion_CILIA_UNCILIA.sct) <- 'RNA'

Idents(adhesion_CILIA_UNCILIA.sct) <- 'celltype'
adhesion_CILIA_UNCILIA.sct.CCL5 <- subset(adhesion_CILIA_UNCILIA.sct,
                                          idents = "U7")
Idents(adhesion_CILIA_UNCILIA.sct.CCL5) <- 'group'
CCL5_diff_IUA <- adhesion_CILIA_UNCILIA.sct.CCL5@misc$CCL5_diff_IUA <- 
  FindMarkers(adhesion_CILIA_UNCILIA.sct.CCL5,
              ident.1 = c('IUA'), 
              ident.2 = c('Normal'), 
              verbose = FALSE,
              logfc.threshold = 0)


CCL5_diff_IUA <- CCL5_diff_IUA %>% 
  mutate(sig=case_when(avg_log2FC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_log2FC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_log2FC > -0.25 | avg_log2FC < 0.25 ~ 'non')) %>%
  mutate(celltype='CCL5 EPI') %>% 
  mutate(gene=row.names(CCL5_diff_IUA)) %>% 
  rownames_to_column("gene1")


p_mmp7 <- ggplot(data = MMP7_diff_IUA,
            aes(x = avg_log2FC,
                y = -log10(p_val_adj),
                color = sig)) + 
  geom_point(size = 2.2) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + 
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-5, 5, by = 2)) + 
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 60),
                     breaks = seq(0, 60, by = 10))+ 
  geom_hline(yintercept = c(-log10(0.05)),size = 0.5,color = "black",lty = "dashed") + 
  geom_vline(xintercept = c(-0.25, 0.25),size = 0.5,color = "black",lty = "dashed") + 
  mytheme

up <- dplyr::filter(MMP7_diff_IUA, sig == 'up') %>% dplyr::distinct(gene, .keep_all = T) %>%
  dplyr::top_n(3, avg_log2FC)
down <- dplyr::filter(MMP7_diff_IUA, sig == 'down') %>% distinct(gene, .keep_all = T) %>%
  top_n(3, -avg_log2FC)

mmp7_vlo <-p_mmp7 + geom_text_repel(data = up,
                    aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_text_repel(data = down,
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene))



p_ccl5 <- ggplot(data = CCL5_diff_IUA,
                 aes(x = avg_log2FC,
                     y = -log10(p_val_adj),
                     color = sig)) + 
  geom_point(size = 2.2) +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + 
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-5, 5, by = 2)) + 
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 60),
                     breaks = seq(0, 60, by = 10))+ 
  geom_hline(yintercept = c(-log10(0.05)),size = 0.5,color = "black",lty = "dashed") + 
  geom_vline(xintercept = c(-0.25, 0.25),size = 0.5,color = "black",lty = "dashed") + 
  mytheme


up <- dplyr::filter(CCL5_diff_IUA, sig == 'up') %>% dplyr::distinct(gene, .keep_all = T) %>%
  dplyr::top_n(3, avg_log2FC)
down <- dplyr::filter(CCL5_diff_IUA, sig == 'down') %>% distinct(gene, .keep_all = T) %>%
  top_n(3, -avg_log2FC)


ccl5_vlo <- p_ccl5 + geom_text_repel(data = up,
                         aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
  geom_text_repel(data = down,
                  aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene))


pttg1_vlo + mmp7_vlo + ccl5_vlo

cowplot::plot_grid(pttg1_vlo,
                   mmp7_vlo , 
                   ccl5_vlo,ncol = 2)


####  go terms

PTTG1_diff_all_upgene <- (PTTG1_diff_IUA %>% dplyr::filter(sig == 'up'))$gene

PTTG1_diff_all_upgene.df <- bitr(PTTG1_diff_all_upgene, fromType = "SYMBOL",
                              toType = c("ENSEMBL", "ENTREZID"),
                              OrgDb = org.Hs.eg.db)

PTTG1_diff_all_upgene.ego <- enrichGO(gene =PTTG1_diff_all_upgene.df$ENTREZID,
                                   OrgDb = org.Hs.eg.db,
                                   ont = "BP",
                                   readable = TRUE,
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05)

PTTG1_diff_all_downgene <- (PTTG1_diff_IUA %>% dplyr::filter(sig == 'down'))$gene

PTTG1_diff_all_downgene.df <- bitr(PTTG1_diff_all_downgene, fromType = "SYMBOL",
                                toType = c("ENSEMBL", "ENTREZID"),
                                OrgDb = org.Hs.eg.db)

PTTG1_diff_all_downgene.ego <- enrichGO(gene =PTTG1_diff_all_downgene.df$ENTREZID,
                                     OrgDb = org.Hs.eg.db,
                                     ont = "BP",
                                     readable = TRUE,
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)



PTTG1_diff_all_upgene.plot <- PTTG1_diff_all_upgene.ego@result[c('GO:0030199',
                                                                 'GO:0032465',
                                                                 'GO:0000302',
                                                                 'GO:0060485',
                                                                 'GO:0044773'),]

PTTG1_diff_all_downgene.plot <- PTTG1_diff_all_downgene.ego@result[c('GO:0071222',
                                                                     'GO:0030216',
                                                                     'GO:0003382',
                                                                     'GO:0048865',
                                                                     'GO:0072089'),]




PTTG1_diff_all_updowngene.plot <- rbind(PTTG1_diff_all_upgene.plot,
                                     PTTG1_diff_all_downgene.plot)

PTTG1_diff_all_updowngene.plot <- PTTG1_diff_all_updowngene.plot %>% 
  mutate(type=case_when(ID == 'GO:0030199' ~ 'up',
                        ID == 'GO:0032465' ~ 'up',
                        ID == 'GO:0000302' ~ 'up',
                        ID == 'GO:0060485' ~ 'up',
                        ID == 'GO:0044773' ~ 'up',
                        ID == 'GO:0071222'  ~ 'down',
                        ID == 'GO:0030216'  ~ 'down',
                        ID == 'GO:0003382'  ~ 'down',
                        ID == 'GO:0048865'  ~ 'down',
                        ID == 'GO:0072089'  ~ 'down'
))

PTTG1_diff_all_updowngene.plot$Description<-factor(PTTG1_diff_all_updowngene.plot$Description,
                                                levels=c("mitotic DNA damage checkpoint signaling",
                                                         "mesenchyme development",
                                                         "response to reactive oxygen species",
                                                         "regulation of cytokinesis",
                                                         "collagen fibril organization",
                                                         "cellular response to lipopolysaccharide",
                                                 "keratinocyte differentiation" ,
                                                 "epithelial cell morphogenesis",
                                                 "stem cell fate commitment" ,
                                                 "stem cell proliferation"     
                                                         ))


ggplot(data=PTTG1_diff_all_updowngene.plot,mapping=aes(x=Description,y=-log10(p.adjust),fill=type))+
  geom_bar(stat="identity") + coord_flip() + theme_test() + scale_fill_manual(values =  c('skyblue','pink'))


####### go terms


U2_diff_all_upgene <- (MMP7_diff_IUA %>% dplyr::filter(sig == 'up'))$gene

U2_diff_all_upgene.df <- bitr(U2_diff_all_upgene, fromType = "SYMBOL",
                              toType = c("ENSEMBL", "ENTREZID"),
                              OrgDb = org.Hs.eg.db)

U2_diff_all_upgene.ego <- enrichGO(gene =U2_diff_all_upgene.df$ENTREZID,
                                   OrgDb = org.Hs.eg.db,
                                   ont = "BP",
                                   readable = TRUE,
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05)

U2_diff_all_downgene <- (MMP7_diff_IUA %>% dplyr::filter(sig == 'down'))$gene

U2_diff_all_downgene.df <- bitr(U2_diff_all_downgene, fromType = "SYMBOL",
                                toType = c("ENSEMBL", "ENTREZID"),
                                OrgDb = org.Hs.eg.db)

U2_diff_all_downgene.ego <- enrichGO(gene =U2_diff_all_downgene.df$ENTREZID,
                                     OrgDb = org.Hs.eg.db,
                                     ont = "BP",
                                     readable = TRUE,
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)


U7_diff_all_upgene <- (CCL5_diff_IUA %>% dplyr::filter(sig == 'up'))$gene

U7_diff_all_upgene.df <- bitr(U7_diff_all_upgene, fromType = "SYMBOL",
                              toType = c("ENSEMBL", "ENTREZID"),
                              OrgDb = org.Hs.eg.db)

U7_diff_all_upgene.ego <- enrichGO(gene =U7_diff_all_upgene.df$ENTREZID,
                                   OrgDb = org.Hs.eg.db,
                                   ont = "BP",
                                   readable = TRUE,
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05)

U7_diff_all_downgene <- (CCL5_diff_IUA %>% dplyr::filter(sig == 'down'))$gene

U7_diff_all_downgene.df <- bitr(U7_diff_all_downgene, fromType = "SYMBOL",
                                toType = c("ENSEMBL", "ENTREZID"),
                                OrgDb = org.Hs.eg.db)

U7_diff_all_downgene.ego <- enrichGO(gene =U7_diff_all_downgene.df$ENTREZID,
                                     OrgDb = org.Hs.eg.db,
                                     ont = "BP",
                                     readable = TRUE,
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)


U2_diff_all_upgene.plot <- U2_diff_all_upgene.ego@result[c('GO:0030199',
                                                           'GO:0030198',
                                                           'GO:0043062',
                                                           'GO:0048144',
                                                           'GO:0001837'),]

U2_diff_all_downgene.plot <- U2_diff_all_downgene.ego@result[c('GO:0042060',
                                                               'GO:0002064',
                                                               'GO:0007249',
                                                               'GO:0050727',
                                                               'GO:0097193'),]


U7_diff_all_upgene.plot <- U7_diff_all_upgene.ego@result[c('GO:0050852',
                                                           'GO:0001909',
                                                           'GO:0001906',
                                                           'GO:0030101',
                                                           'GO:0001913'),]


U7_diff_all_downgene.plot <- U7_diff_all_downgene.ego@result[c('GO:0002064',
                                                               'GO:0042060',
                                                               'GO:0060562',
                                                               'GO:0048732',
                                                               'GO:0051403'),]


U7_diff_all_updowngene.plot <- rbind(U7_diff_all_upgene.plot,
                                     U7_diff_all_downgene.plot)

U2_diff_all_updowngene.plot <- rbind(U2_diff_all_upgene.plot,
                                     U2_diff_all_downgene.plot)


U2_diff_all_updowngene.plot <- U2_diff_all_updowngene.plot %>% 
  mutate(type=case_when(ID == "GO:0030199"  ~ 'up',
                        ID == "GO:0030198"  ~ 'up',
                        ID == "GO:0043062"  ~ 'up',
                        ID == "GO:0048144"  ~ 'up',
                        ID == "GO:0001837"  ~ 'up',
                        ID == "GO:0042060"  ~ 'down',
                        ID == "GO:0002064"  ~ 'down',
                        ID == "GO:0007249"  ~ 'down',
                        ID == "GO:0050727"  ~ 'down',
                        ID == "GO:0097193"  ~ 'down'
  ))



U7_diff_all_updowngene.plot <- U7_diff_all_updowngene.plot %>% 
  mutate(type=case_when(ID == "GO:0050852"  ~ 'up',
                        ID == "GO:0001909"  ~ 'up',
                        ID == "GO:0001906"  ~ 'up',
                        ID == "GO:0030101"  ~ 'up',
                        ID == "GO:0001913"  ~ 'up',
                        ID == "GO:0002064"  ~ 'down',
                        ID == "GO:0042060"  ~ 'down',
                        ID == "GO:0060562"  ~ 'down',
                        ID == "GO:0048732"  ~ 'down',
                        ID == "GO:0051403"  ~ 'down'))

U7_diff_all_updowngene.plot$Description<-factor(U7_diff_all_updowngene.plot$Description,
                                                levels=c('T cell mediated cytotoxicity',
                                                         'natural killer cell activation',
                                                         'cell killing',
                                                         'leukocyte mediated cytotoxicity',
                                                         'T cell receptor signaling pathway',
                                                         'epithelial cell development',
                                                         'wound healing',
                                                         'epithelial tube morphogenesis',
                                                         'gland development',
                                                         'stress-activated MAPK cascade'))

U2_diff_all_updowngene.plot$Description<-factor(U2_diff_all_updowngene.plot$Description,
                                                levels=c(
                                                  'epithelial to mesenchymal transition',
                                                  'fibroblast proliferation',
                                                  'extracellular structure organization',
                                                  'extracellular matrix organization',
                                                  'collagen fibril organization',
                                                  'wound healing',
                                                  'epithelial cell development',
                                                  'I-kappaB kinase/NF-kappaB signaling',
                                                  'regulation of inflammatory response',
                                                  'intrinsic apoptotic signaling pathway'))


CCL5_barplot <- ggplot(data=U7_diff_all_updowngene.plot,mapping=aes(x=Description,y=-log10(p.adjust),fill=type))+
  geom_bar(stat="identity") + coord_flip() + theme_test() + scale_fill_manual(values =  c('skyblue','pink'))


MMP7_barplot <- ggplot(data=U2_diff_all_updowngene.plot,mapping=aes(x=Description,y=-log10(p.adjust),fill=type))+
  geom_bar(stat="identity") + coord_flip() + theme_test() + scale_fill_manual(values =  c('skyblue','pink'))


PTTG1_barplot <- ggplot(data=PTTG1_diff_all_updowngene.plot,mapping=aes(x=Description,y=-log10(p.adjust),fill=type))+
  geom_bar(stat="identity") + coord_flip() + theme_test() + scale_fill_manual(values =  c('skyblue','pink'))


PTTG1_barplot + MMP7_barplot + CCL5_barplot

  
  pttg1_vlo + mmp7_vlo + ccl5_vlo

cowplot::plot_grid(pttg1_vlo,
                   mmp7_vlo , 
                   ccl5_vlo,ncol = 2)



#########  scores  stem cell

stem <- list(unique(stem_cell_differentiation$V1))

STEM <- AddModuleScore(
  object = adhesion_CILIA_UNCILIA.sct.MMP7,
  features = stem,
  ctrl = 100, #默认值是100
  name = 'STEM'
)
colnames(STEM@meta.data)[19] <- 'stem' 
data<- FetchData(STEM, vars = c("celltype","group","stem"))

data$group <- factor(data$group, 
                     levels=c("Normal",
                              "IUA"))

stem_mmp7 <- ggplot(data, aes(x=group,y=stem), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = stem,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  ) 





###### EMT

emt <- list(unique(EMT$V1))

EMT <- AddModuleScore(
  object = adhesion_CILIA_UNCILIA.sct.CCL5,
  features = emt,
  ctrl = 100, #默认值是100
  name = 'EMT'
)
colnames(EMT@meta.data)[19] <- 'emt' 
data<- FetchData(EMT, vars = c("celltype","group","emt"))

data$group <- factor(data$group, 
                     levels=c("Normal",
                              "IUA"))

emt_ccl5 <- ggplot(data, aes(x=group,y=emt), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = emt,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  ) 









EMT_pttg1 <- ggplot(data, aes(x=group,y=EMT), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = EMT,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  ) 


stem_pttg1




########   wnt

wnt <- list(unique(wnt$V1))

WNT <- AddModuleScore(
  object = adhesion_CILIA_UNCILIA.sct.MMP7,
  features = wnt,
  ctrl = 100, #默认值是100
  name = 'WNT'
)
colnames(WNT@meta.data)[19] <- 'wnt' 
data<- FetchData(WNT, vars = c("celltype","group","wnt"))

data$group <- factor(data$group, 
                     levels=c("Normal",
                              "IUA"))


wnt_mmp7 <- ggplot(data, aes(x=group,y=wnt), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = wnt,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  ) 


########   proflieratio

proliferation <- list(unique(proliferation$V1))

PRO <- AddModuleScore(
  object = adhesion_CILIA_UNCILIA.sct.CCL5,
  features = proliferation,
  ctrl = 100, #默认值是100
  name = 'PRO'
)
colnames(PRO@meta.data)[19] <- 'proliferation' 
data<- FetchData(PRO, vars = c("celltype","group","proliferation"))

data$group <- factor(data$group, 
                     levels=c("Normal",
                              "IUA"))


pro_ccl5 <- ggplot(data, aes(x=group,y=proliferation), fill=group) +
  geom_violin(aes(color=group),position = position_dodge(0.45),width=1)+  
  scale_color_manual(values = c("#3288BD","#D53E4F")) +
  geom_boxplot(aes(color=group),width=.15,position = position_dodge(.45))+
  scale_fill_manual(values = c("#3288BD","#D53E4F")) +
  xlab(NULL) + ylab(NULL) + theme_classic2()+
  stat_compare_means(data = data,aes(x=group,y = proliferation,fill=group),label = 'p.signif')  + NoLegend()+
  theme(
    axis.text = element_text(size = 12), #轴标签大小调整
    axis.title.x = element_text(size = 13), #x轴标题大小调整
    legend.title = element_text(size = 10), #图例标题大小调整
    legend.text = element_text(size = 10), #图例标签大小调整
  ) 













wnt_plot <- wnt_pttg1 + wnt_mmp7+ wnt_ccl5

stem_plot <- stem_pttg1 + stem_mmp7 + stem_ccl5

emt_plot <- emt_PTTG1 + emt_mmp7 + emt_ccl5


wnt_plot + stem_plot + emt_plot

cowplot::plot_grid(wnt_plot,
                   stem_plot , 
                   emt_plot,nrow = 3)


###########  NK cells


DefaultAssay(adhesion_NK.sct) <- 'RNA'

Idents(adhesion_NK.sct) <- 'celltype'
adhesion_NK.sct.NK3 <- subset(adhesion_NK.sct,
                              idents = "NK3")

Idents(adhesion_NK.sct.NK3) <- 'group'
NK3_diff_IUA <- adhesion_NK.sct.NK3@misc$NK3_diff_IUA <- 
  FindMarkers(adhesion_NK.sct.NK3,
              ident.1 = c('IUA'), 
              ident.2 = c('Normal'), 
              verbose = FALSE,
              logfc.threshold = 0)


NK3_diff_IUA <- NK3_diff_IUA %>% 
  mutate(sig=case_when(avg_log2FC <= -0.25 & p_val_adj < 0.05 ~ 'down',
                       avg_log2FC >=  0.25 & p_val_adj < 0.05 ~ 'up',
                       p_val_adj > 0.05 ~ 'non',
                       avg_log2FC > -0.25 | avg_log2FC < 0.25 ~ 'non')) %>%
  mutate(celltype='NK3') %>% 
  mutate(gene=row.names(NK3_diff_IUA)) %>% 
  rownames_to_column("gene1")

NK3_diff_IUA.gene.ego <- NK3_diff_IUA  %>%
  dplyr::filter(sig != 'non')  %>%
  dplyr::pull(gene, name=avg_log2FC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)




NK3_diff_IUA.upgene.ego <- NK3_diff_IUA  %>%
  dplyr::filter(sig == 'up')  %>%
  dplyr::pull(gene, name=avg_log2FC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)


NK3_diff_IUA.downgene.ego <- NK3_diff_IUA  %>%
  dplyr::filter(sig == 'down')  %>%
  dplyr::pull(gene, name=avg_log2FC) %>% 
  bitr(fromType = "SYMBOL",
       toType = c("ENSEMBL", "ENTREZID"),
       OrgDb = org.Hs.eg.db) %>% 
  .$ENTREZID %>%
  enrichGO(gene = .,
           OrgDb = org.Hs.eg.db,
           ont = "BP",
           readable = TRUE,
           pvalueCutoff = 0.05,
           qvalueCutoff = 0.05)

geneList <- NK3_diff_IUA  %>%
  dplyr::filter(sig != 'non')  %>%
  dplyr::pull(gene, name=avg_log2FC)

#########  scores  hypoxia  ros  GO:0072593 GO:0001666


NK3_diff_ros <- NK3_diff_IUA %>% dplyr::filter(sig != 'non') %>%
  dplyr::filter(gene %in% ros[[1]])


NK3_diff_cytoxi <- NK3_diff_IUA %>% dplyr::filter(sig != 'non') %>%
  dplyr::filter(gene %in% cytoxi )



geneList <- NK3_diff_IUA  %>%
  dplyr::filter(sig != 'non')  %>%
  dplyr::pull(gene, name=avg_log2FC)


cnetplot(NK3_diff_IUA.gene.ego, categorySize="pvalue", 
         foldChange=geneList,showCategory = 20,
         color_category = "#E5C494",)
















