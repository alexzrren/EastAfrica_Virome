library(vegan)
library(proxy)
library(readxl)
library(tidyverse)
library(geosphere)
library(ape)
library(reshape2)
library(scales)

sample.df <- read_xlsx('/Users/zirui/BGI_Projects/WIV_EastAfrica/VirSpecPhylo_2Spec/sampleinfo_use.xlsx')
requant.df <- read_xlsx('/Users/zirui/BGI_Projects/WIV_EastAfrica/Requantification_new/EastAfrica_Requantification_0724_rerun_v1_batfecesfiltedspecies_filtered.xlsx') |>
  mutate(vspecname = str_sub(vsname, end=-2))

revised.gps.df <- read_xlsx('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/Tables/1_Metdata_v2.1.xlsx') |>
  select(sn, revised_gps)

sample.df <- sample.df |>
  left_join(revised.gps.df, by='sn')


host_species <- unique((requant.df |> filter(host_species!="Otomops martiensseni"))$host_species)

result <- data.frame()
for( hspec in host_species){
  hspec_fasta <- paste(gsub(' ', '_', hspec), 'fasta', sep='.')
  
  requant.df.hspec <- filter(requant.df, host_species==hspec)
  positive.sn <- requant.df.hspec$sampleindex
  sample.df.hspec <- sample.df |> filter(Adjusted_Host_Species==hspec) |> filter(sn %in% positive.sn)
  
  # Virome Jaccard distance matrix
  sample.virome <- as.matrix(table(requant.df.hspec$sampleindex, requant.df.hspec$vspecname))
  virome.jaccard.matx <- as.matrix(dist(sample.virome, method="Jaccard"))
  
  # GPS matrix
  sample_gps_df <- sample.df.hspec |>
    select(sn, revised_gps) |>
    separate(col = revised_gps, into = c("Lat", "Lon"), sep = ",") # separate gps coord string into longitude and latitude
  
  sample.geodist.matx <- sample_gps_df |>
    select(Lon, Lat) |>
    as.data.frame() |>
    mutate(
      Lon = as.numeric(Lon),
      Lat = as.numeric(Lat)
    ) |>
    select(Lon, Lat) |>
    distm(fun=distVincentyEllipsoid) # calc spherical distance
  
  rownames(sample.geodist.matx) <- sample_gps_df$sn # add row and col names
  colnames(sample.geodist.matx) <- sample_gps_df$sn
  
  sample.geodist.matx <- sample.geodist.matx[rownames(virome.jaccard.matx), colnames(virome.jaccard.matx)] # sub-matrix 
  
  # Host genetic distance matrix
  host.genodist.matx <- dist.dna(read.FASTA(paste('/Users/zirui/BGI_Projects/WIV_EastAfrica/Mitocondrial_Heterozygosity/snp_alignment', hspec_fasta, sep='/')),
                                 as.matrix = T, pairwise.deletion = T)
  
  rownames(host.genodist.matx) <- unlist(lapply(rownames(host.genodist.matx), function(e) str_split(e, '_')[[1]][1])) # sample name converting 
  colnames(host.genodist.matx) <- unlist(lapply(colnames(host.genodist.matx), function(e) str_split(e, '_')[[1]][1])) # eg. R2101013270_1_2021-04-16 -> R2101013270
  
  host.genodist.matx <- host.genodist.matx[rownames(virome.jaccard.matx), colnames(virome.jaccard.matx)] 
  
  vg.h <- mantel.partial(virome.jaccard.matx, sample.geodist.matx, host.genodist.matx, method = "spearman", permutations = 999, 
                         strata = NULL, na.rm = T, parallel = 8)
  vh.g <- mantel.partial(virome.jaccard.matx, host.genodist.matx, sample.geodist.matx, method = "spearman", permutations = 999, 
                         strata = NULL, na.rm = T, parallel = 8)
  vg <- mantel(virome.jaccard.matx, sample.geodist.matx, method = "spearman", permutations = 999, 
               strata = NULL, na.rm = T, parallel = 8)
  vh <- mantel(virome.jaccard.matx, host.genodist.matx, method = "spearman", permutations = 999, 
               strata = NULL, na.rm = T, parallel = 8)
  hg <- mantel(host.genodist.matx, sample.geodist.matx, method = "spearman", permutations = 999, 
               strata = NULL, na.rm = T, parallel = 8)
  #sampled.sn.list <- sample(rownames(virome.jaccard.matx), 50)
  #sampled.virome.matx <- virome.jaccard.matx[sampled.sn.list, sampled.sn.list]
  #sampled.geodist.matx <- sample.geodist.matx[sampled.sn.list, sampled.sn.list]
  
  row <- data.frame(hspec, vg.h$statistic, vg.h$signif,
                    vh.g$statistic, vh.g$signif,
                    vg$statistic, vg$signif,
                    vh$statistic, vh$signif,
                    hg$statistic, hg$signif)
  
  result <- rbind(result, row)
  
  print(hspec)
}
df.fordraw <-  result |>
  mutate(
    vg.h.padj = p.adjust(vg.h.signif, method='BH', n=6),
    vh.g.padj = p.adjust(vh.g.signif, method='BH', n=6),
    vg.padj = p.adjust(vg.signif, method='BH', n=6),
    vh.padj = p.adjust(vh.signif, method='BH', n=6),
    hg.padj = p.adjust(hg.signif, method='BH', n=6)
  ) |>
  melt() |>
  mutate(
    value_type = as.factor(unlist(lapply(strsplit(as.character(variable), '.', fixed = T), function(e) e[length(e)]))),
    compare_group = as.factor(unlist(lapply(strsplit(as.character(variable), '.', fixed = T), function(e) paste(e[1:length(e)-1],collapse='|'))))
  ) |>
  filter(value_type != 'signif') |>
  select(-variable)|>
  dcast(hspec+compare_group~value_type) |>
  mutate(signif = cut(padj, breaks=c(0,0.01,0.05,1), labels = c("P<0.01","P<0.05","NS"))) 


ggplot(df.fordraw)+
  geom_point(aes(x=fct_relevel(compare_group, c('vg','hg','vg|h','vh','vh|g')), y=hspec, fill=statistic, 
                 size=-log10(padj), color=signif, stroke=1), shape=21)+
  scale_fill_gradientn(colours = c("DodgerBlue1","white","DarkOrange1"), 
                       values = rescale(c(-0.2,0,1)),
                       guide = "colorbar", limits=c(-0.2,1), breaks = c(-0.2,0,0.2,0.4,0.6,0.8,1.0))+
  scale_color_manual(values = c("NavyBlue","royalblue1","lightgray"))+
  theme_bw()+
  labs(x = 'Compare variable', y = 'Bat species', title='(Partial) Mantel Test')+
  theme(
    axis.text.x = element_text(angle = 30, hjust=1, vjust = 1, color='black'),
    axis.text.y = element_text(color='black', face='italic', size=10),
    axis.title.x = element_text(color='black', size=12, face='bold'),
    axis.title.y = element_text(color='black', size=12, face='bold')
  )
ggsave('/Users/zirui/BGI_Projects/WIV_EastAfrica_Clean/Plot/ExtendedFig/ExtendedFig9/ExtFig.9a.pdf', height=4.5, width = 4.8)
#mutate(compare = tail(strsplit(as.character(variable), '.', fixed = T), n=1))

