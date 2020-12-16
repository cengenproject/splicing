# Analyze results of SplAdder and generate figure 6D bottom/left half
# 
# SplAdder was previously run on cluster and all pairwise tests performed.
# Here the test results are loaded for further analysis and to generate the 
# figures for the paper.
# Note: updated on 2020-07-15 for using sg_min_edge_count = 3 (vs 10 by default)
# Note: updated on 2020-07-26 to use alignments "AB_bss2" to be consistent with other figure

# Initializations ----
library(tidyverse)

library(wbData) # defines s2i() and i2s() to convert between Wormbase IDs and gene symbols
gids <- wb_load_gene_ids(274)

dpar <- par()




# Gene coordinates (Wormbase)
dfcoords <- wb_load_gene_coords(274) %>%
  select(gene_id, chr, start, end) %>%
  mutate(length = end - start) %>%
  data.frame(row.names = coords$gene_id)




# Read the SplAdder result table for that pair of neurons and that event type
read_spl_results <- function(neurA, neurB, type){
  read_tsv(paste0("data/",
                  "testing_",neurA,"_vs_",neurB,"/",
                  "test_results_C3_",type,".tsv"),
           col_types = cols(
             event_id = col_character(),
             gene = col_character(),
             p_val = col_double(),
             p_val_adj = col_double(),
             mean_event_count_A = col_double(),
             mean_event_count_B = col_double(),
             log2FC_event_count = col_double(),
             mean_gene_exp_A = col_double(),
             mean_gene_exp_B = col_double(),
             log2FC_gene_exp = col_double()
           )
  )
}



# Read and process data ----



# define neuron modalities (sensory, inter or motor)
list_modalities <- read_csv("data/neur_modalities.csv")
list_neurons <- list_modalities$neuron
list_event_types <- c("alt_3prime", "alt_5prime", "exon_skip",
                      "intron_retention", "mult_exon_skip", "mutex_exons")




#~~ Read quantification----

# read results of all SplAdder pairwise tests
all_res <- tibble()
for(i in 1:(length(list_neurons)-1)){
  for(j in (i+1):length(list_neurons)){
    all_res <- bind_rows(all_res,
                         map_dfr(list_event_types,
                                 ~read_spl_results(list_neurons[i], list_neurons[j], .x)) %>%
                           add_column(neurA=list_neurons[i], neurB=list_neurons[j]))
  }
}

# Return the modality for that neuron
neur2modal <- function(neuron)
  return(list_modalities$modality[match(neuron, list_modalities$neuron)])

# Concatenate string pair in alphabetical order
alpha_order <- function(aaa, bbb){
  map2_chr(aaa, bbb, ~ paste(pmin(.x, .y), pmax(.x, .y), sep = "-"))
}


# Annotate data
res <- all_res %>%
  mutate(padj = p.adjust(p_val, method = "BH")) %>%
  filter(padj < 0.1) %>%
  mutate(symbol = i2s(gene, gids),
         row = row_number(),
         neuron_pair = paste(neurA, neurB, sep="-"),
         modalities = alpha_order(neur2modal(neurA), neur2modal(neurB)),
         event_type = str_match(event_id, "^([a-z_35]+)_[0-9]{1,4}$")[,2], .after=gene) %>%
  filter(event_type != "intron_retention")
#> 348 significant events at 0.1 FDR
#> 228 at 0.05 FDR
#> 112 at 0.01 FDR



#~~ Event-level annotation ----

# Format big number (e.g. genomic position)
pfrm <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)

# Using the coordinates of the exons from the gff3 files, return a single string describing position
get_centr_coord <- function(df, ev_type){
  stopifnot(ncol(df) == 3)
  switch(ev_type,
         `alt_3prime+`=,
         `alt_5prime-`={
           stopifnot(nrow(df) == 4)
           return(paste0(paste(pfrm(df[2,1:2]), collapse = "-"),
                         "_vs_",
                         paste(pfrm(df[4,1:2]), collapse = "-")))
         },
         `alt_3prime-`=,
         `alt_5prime+`={
           stopifnot(nrow(df) == 4)
           return(paste0(paste(pfrm(df[1,1:2]), collapse = "-"),
                         "_vs_",
                         paste(pfrm(df[3,1:2]), collapse = "-")))
         },
         `exon_skip+`=,
         `exon_skip-`={
           stopifnot(nrow(df) == 5)
           return(paste(pfrm(df[4,1:2]), collapse = "-"))
         },
         `intron_retention+`=,
         `intron_retention-`={
           stopifnot(nrow(df) == 3)
           return(paste(pfrm(df[3,1:2]), collapse = "-"))
         },
         `mult_exon_skip+`=,
         `mult_exon_skip-`={
           stopifnot(nrow(df) > 5)
           return(paste(apply(df[4:(nrow(df)-1), 1:2], 1,
                              function(l) paste(pfrm(l), collapse = "-")),
                        collapse = "_and_"))
         },
         `mutex_exons+`=,
         `mutex_exons-`={
           stopifnot(nrow(df) == 6) # here always 6 rows, not sure if generalizes to other datasets
           return(paste0(paste(pfrm(df[2,1:2]), collapse = "-"),
                         "_vs_",
                         paste(pfrm(df[5,1:2]), collapse = "-")))
         })
  stop("Unknown error.")
}

# read spladder annotation from GFF3
all_events <- map_dfr(list_event_types,
                      ~ read_tsv(paste0(
                        "data/merge_graphs_",
                        .x, "_C3.confirmed.gff3"),
                        skip=1,
                        col_names = c("chr","event_type", "feature",
                                      "start", "end", "score", "strand",
                                      "frame", "attributes"),
                        col_types = cols(chr = col_character(),
                                         event_type = col_character(),
                                         feature = col_character(),
                                         start = col_double(),
                                         end = col_double(),
                                         score = col_character(),
                                         strand = col_character(),
                                         frame = col_character(),
                                         attributes = col_character())))
# annotate data frame
spl_events_pos <- all_events %>%
  filter(feature == "exon") %>%
  extract(attributes,
          into = c("event_id", "isoform"),
          "^Parent=([a-z_35]+\\.[0-9]{1,4})_(iso[12])$") %>%
  nest(coords = c(start, end, isoform)) %>%
  mutate(event_id = str_replace(event_id,"\\.","_"),
         event_type_strand = paste0(event_type,strand),
         central_coord = map2_chr(coords, event_type_strand, get_centr_coord))




# Plot: Fig 6D ----

pair_cnts <- res %>%
  group_by(neuron_pair, neurA, neurB) %>%
  count() %>%
  ungroup() %>%
  select(-neuron_pair)


# Make sure we have both sides and make matrix
pair_cnts_square <- pair_cnts %>%
  mutate(neurB2 = neurA,
         neurA = neurB,
         neurB = neurB2) %>%
  select(- neurB2) %>%
  full_join(pair_cnts) %>%
  pivot_wider(names_from = neurB, values_from = n)

mat_cnts <- as.matrix(select(pair_cnts_square, - neurA))
rownames(mat_cnts) <- pair_cnts_square$neurA
mat_cnts[is.na(mat_cnts)] <- 0


# set order of neurons
mat_cnts <- mat_cnts[, c("ASG", "AVG", "AWB", "AWA", "AVE", "PVD", "DD", "VD")]
mat_cnts <- mat_cnts[colnames(mat_cnts),]

# Be consistent with Fig 5, use ggplot2 default palette for modality encoding
modal_colors <- hcl(seq(15,365,length.out = 5), l=65,c=100)[1:3]
names(modal_colors) <- c("Interneuron", "Sensory", "Motor")

# Prepare margin annotations for modality
annot <- data.frame(Modality = list_modalities$modality,
                    row.names = list_modalities$neuron)
annot_colors <- list(Modality = modal_colors)


# Figure 6D, bottom/left half
pheatmap::pheatmap(mat_cnts,
                   scale = "none",
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 9,
                                                                     name = "Blues")
                   )(max(mat_cnts)),
                   border_color = NA,
                   annotation_row = annot,
                   annotation_col = annot,
                   annotation_names_row = FALSE,
                   annotation_names_col = FALSE,
                   annotation_colors = annot_colors,
                   width = 5,
                   height = 4,
                   filename = "output/figures/splicing_heatmap_pairs.pdf",
                   fontsize = 6,
                   fontsize_row = 11,
                   fontsize_col = 11)



# Table S12 ----

all_res %>%
  mutate(padj = p.adjust(p_val, method = "BH")) %>%
  filter(padj < 0.05) %>%          #Note: p<0.05 here, vs p<0.1 above for the figure
  mutate(symbol = i2s(gene, gids),
         row = row_number(),
         neuron_pair = paste(neurA, neurB, sep="-"),
         modalities = alpha_order(neur2modal(neurA), neur2modal(neurB)),
         event_type = str_match(event_id, "^([a-z_35]+)_[0-9]{1,4}$")[,2], .after=gene) %>%
  filter(event_type != "intron_retention") %>%
  left_join(dfcoords, by=c("gene" = "gene_id")) %>%
  mutate(spacing = round(0.05*length),
         igv = paste0(chr,":",
                      start-spacing,"-",
                      end+spacing),
         jbrowse = paste0(chr,":",
                          start-spacing,"..",
                          end+spacing)) %>%
  arrange(padj) %>%
  group_by(symbol) %>%
  nest() %>%
  mutate(minp = map_dbl(data, ~min(.x$padj))) %>%
  arrange(minp) %>%
  unnest(everything()) %>%
  left_join(spl_events_pos, by = "event_id")  %>%
  select(event_id, symbol, event_type = event_type.x, gene_id = gene, strand,
         neuronA=neurA, neuronB=neurB, gene_position = igv, event_position = central_coord) %>%
  write_csv("output/tables/table_S12_events.csv")



