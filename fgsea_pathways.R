library(fgsea)
library(tidyverse)


res <- read.table("~/Desktop/projects/data/fgsea/XXX",header = TRUE, sep = "")
head(res)
colnames(res)[1]  <- "SYMBOL"    # change column name for x column


### Further, all you’ll care about later on is the gene symbol and the test statistic.####
##Get just those, and remove the NAs. Finally, if you have multiple test statistics for###
####the same symbol, you’ll want to deal with that in some way. Here I’m just averaging them.####
res2 <- res %>% 
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(SYMBOL) %>% 
  summarize(stat=mean(stat))
res2
ranks <- deframe(res2)
head(ranks, 20)




########################################################################
#########################################################################


fgseaRes_reactome=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c2/c2.cp.reactome.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy <- fgseaRes_reactome %>%
  as_tibble() %>%
  arrange(desc(NES))
readr::write_csv(fgseaRes_reactome, "fgsea_reactome_XXX.csv")

fgseaRes_kegg=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c2/c2.cp.kegg.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy <- fgseaRes_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))


# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

fgseaResTidy=data.frame(fgseaResTidy)

p=ggplot(fgseaResTidy%>%filter(padj <0.01) , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Reactome_Pathways_XXX") + 
  theme_minimal()
pdf("Reactome_Pathways_XXX_only_True.pdf", width=20,height=50)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Reactome_Pathways_XXX") + 
  theme_minimal()
pdf("Reactome_Pathways_XXX_all.pdf", width=20,height=120)
print(p)
dev.off()


########################################################################
#########################################################################


fgseaRes_pid=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c2/c2.cp.pid.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy <- fgseaRes_pid %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
p=ggplot(fgseaResTidy%>%filter(padj <0.01) , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="PID_Pathways_XXX") + 
  theme_minimal()
pdf("PID_Pathways_XXX_only_True.pdf", width=10,height=10)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="PID_Pathways_XXX") + 
  theme_minimal()
pdf("PID_Pathways_XXX_all.pdf", width=10,height=20)
print(p)
dev.off()

########################################################################
#########################################################################


fgseaRes_kegg=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c2/c2.cp.kegg.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy <- fgseaRes_kegg %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
p=ggplot(fgseaResTidy%>%filter(padj <0.01) , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG_Pathways_XXX") + 
  theme_minimal()
pdf("KEGG_Pathways_XXX_only_True.pdf", width=9,height=10)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="KEGG_Pathways_XXX") + 
  theme_minimal()
pdf("KEGG_Pathways_XXX_all.pdf", width=10,height=20)
print(p)
dev.off()


########################################################################
#########################################################################


fgseaRes_biocarta=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c2/c2.cp.biocarta.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy <- fgseaRes_biocarta %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
p=ggplot(fgseaResTidy%>%filter(padj <0.01) , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="biocarta_Pathways_XXX") + 
  theme_minimal()
pdf("biocarta_Pathways_XXX_only_True.pdf", width=10,height=10)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="biocarta_Pathways_XXX") + 
  theme_minimal()
pdf("biocarta_Pathways_XXX_all.pdf", width=10,height=30)
print(p)
dev.off()



########################################################################
###############################GO ANNONTATIONS #########################
########################################################################
fgseaRes_GO_MF=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c5/c5.go.mf.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)

fgseaResTidy <- fgseaRes_GO_MF %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
p= ggplot(fgseaResTidy%>%filter(padj <0.01), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Molecular Function ontology related pathways") + 
  theme_minimal()
pdf("GO_MF_XXX_only_True.pdf", width=20,height=30)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Molecular Function ontology related pathways") + 
  theme_minimal()
pdf("GO_MF_XXX_all.pdf", width=20,height=80)
print(p)
dev.off()

######################################################################


fgseaRes_GO_CC=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c5/c5.go.cc.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)

fgseaResTidy <- fgseaRes_GO_CC %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
p= ggplot(fgseaResTidy%>%filter(padj <0.01), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Cellular Component ontology related pathways") + 
  theme_minimal()
pdf("GO_CC_XXX_only_True.pdf", width=20,height=40)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Cellular Component ontology related pathways") + 
  theme_minimal()
pdf("GO_CC_XXX_all.pdf", width=15,height=90)
print(p)
dev.off()

######################################################################


fgseaRes_GO_bp=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c5/c5.go.bp.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)

fgseaResTidy <- fgseaRes_GO_bp %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
p= ggplot(fgseaResTidy%>%filter(padj <0.01), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Biological Process ontology related pathways") + 
  theme_minimal()
pdf("GO_BP_XXX_only_True.pdf", width=20,height=160)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GO Biological Process ontology related pathways") + 
  theme_minimal()
pdf("GO_BP_XXX_all.pdf", width=20,height=200)
print(p)
dev.off()

######################################################################

#





#########################################################################


fgseaRes_wiki=fgsea(pathways=gmtPathways("~/Desktop/projects/fgsea/c2/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt"), ranks, nperm=1000) %>% 
  as_tibble() %>% 
  arrange(padj)
fgseaResTidy <- fgseaRes_wiki %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
p=ggplot(fgseaResTidy%>%filter(padj <0.01) , aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Wiki_Pathways_XXX") + 
  theme_minimal()
pdf("Wiki_Pathways_XXX_only_True.pdf", width=13,height=20)
print(p)
dev.off()
p=ggplot(fgseaResTidy,aes(reorder(pathway, NES), NES) ) +
  geom_col(aes(fill=padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Wiki_Pathways_XXX") + 
  theme_minimal()
pdf("WIKI_Pathways_XXX_all.pdf", width=13,height=60)
print(p)
dev.off()
pathways.hallmark %>% 
  enframe("pathway", "SYMBOL") %>% 
  unnest() %>% 
  inner_join(res, by="SYMBOL")






















































