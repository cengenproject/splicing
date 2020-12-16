# Load splicing graph to identify novel exons


library(tidyverse)
library(GenomicRanges)
library(wbData)

gids <- wb_load_gene_ids(274)


# Load data ----

#~~ Get Wormbase exon annotation ----
wb_exons <- wb_load_exon_coords(274) %>%
  select(-c(exon_id, transcript_id, exon_number, ends_with("biotype"),position)) %>%
  distinct() %>%
  mutate(symbol = i2s(gene_id, gids),
         full_pos = paste0(chr,":", start,"-", end))


wb_pos <- unique(wb_exons$full_pos)

# 651 exons annotated as being in 2 genes/strands with same start/end
wb_exons %>% group_by(full_pos) %>% summarize(n=n()) %>% pull(n) %>% table()





#~~ Read raw splicing graph vertices ----

# Exported with Python from SplAdder's splicing graph pickle
vertice_list <- readLines("data/200726_vertices.csv")

n_genes <- 46909 #nb of genes found with Python
stopifnot(length(vertice_list)/3 == n_genes)

spl_exons <- tibble(
  gene_id = str_split_fixed(vertice_list[(0:(n_genes - 1))*3+1], " ", 4)[,2],
  chr = str_split_fixed(vertice_list[(0:(n_genes - 1))*3+1], " ", 4)[,3],
  strand = str_split_fixed(vertice_list[(0:(n_genes - 1))*3+1], " ", 4)[,4],
  start = str_split(vertice_list[(0:(n_genes - 1))*3+2], ", "),
  end = str_split(vertice_list[(0:(n_genes - 1))*3+3], ", ")) %>%
  unnest(c(start, end)) %>%
  filter(! (start == "" & end == "")) %>%
  distinct() %>%
  mutate(start = as.integer(start) + 1L,
         end = as.integer(end),
         symbol = i2s(gene_id, gids),
         full_pos = paste0(chr,":", start,"-", end))

spl_pos <- unique(spl_exons$full_pos)
length(spl_pos)
# => note 655 exons are shared by 2 genes (many look like ncRNA on both strands from modEncode)
spl_exons %>% group_by(full_pos) %>% summarize(n=n()) %>% pull(n) %>% table()


# Processing ----


#~~ Novel exons ----

# Compare vertices dumped vs wormbase exons
table(wb_pos %in% spl_pos)
table(! spl_pos %in% wb_pos)
# => 3,860 potential novel exons


# Detect if a new exon is the assembly of existing exons with one or several IR
new_ex_gr <- GRanges(spl_exons)
known_ex_gr <- GRanges(wb_exons)
novel_exons <- spl_exons[which(countOverlaps(new_ex_gr, known_ex_gr, type = "start") *
                  countOverlaps(new_ex_gr, known_ex_gr, type = "end") == 0),] #countOverlaps returns 0 if not found

# => 2,142 exons do not share either start or end with known ones




# Check: difference btw filtering on full exon coords vs start and end separately

not_exact <- spl_exons[which(! spl_exons$full_pos %in% wb_exons$full_pos),]
# -> these exons are not known by exact coords

table(not_exact$full_pos %in% novel_exons$full_pos)
# -> 2138 are not known, as they have an unanot. start or end. 1722 are not exactly known, but share both start and end with annotated exons

not_exact[which(! not_exact$full_pos %in% novel_exons$full_pos), ]
# these exons are not exactly = known, but their start and end known (e.g. IR, often small != in UTR)
# eg rab-11.1, WB has an exon at 108670-109446 and two at 108680-xx + xx-109446, but not 108680-109446

table(! novel_exons$full_pos %in% not_exact$full_pos)
# -> 4 do not have start+end known, but have same full position as known.
# Actually they are on != strands in cases where 2 genes are on opposite strands and could share these exons.




#~~ Remove exons with small differences from known ----

# Filter on overlap size to find exons that are separate from known
novel_ex_gr <- GRanges(novel_exons)
mcols(known_ex_gr) <- NULL

length(countOverlaps(novel_ex_gr, known_ex_gr)) # 2142 novel exons
table(countOverlaps(novel_ex_gr, known_ex_gr)) #42 exons no overlap, 2040 single, 60 two or more overlaps
sum(countOverlaps(novel_ex_gr, known_ex_gr)) # 2255 pairs (from 2100 exons with overlap)


opairs <- mergeByOverlaps(novel_ex_gr, known_ex_gr)
opairs$prop_overlap <- width(pintersect(opairs$novel_ex_gr, opairs$known_ex_gr))/width(opairs$novel_ex_gr)

sep_novel_exons <- opairs %>%
  as_tibble() %>%
  group_by(novel_ex_gr.start, novel_ex_gr.end, novel_ex_gr.seqnames, novel_ex_gr.strand,
           novel_ex_gr.gene_id, novel_ex_gr.symbol) %>%
  summarize(prop_ol = max(prop_overlap)) %>%
  ungroup() %>%
  filter(prop_ol < .9) %>%
  dplyr::select(gene_id = novel_ex_gr.gene_id, chr = novel_ex_gr.seqnames,
                start = novel_ex_gr.start, end = novel_ex_gr.end,
                strand = novel_ex_gr.strand, symbol = novel_ex_gr.symbol) %>%
  mutate(full_pos = paste0(chr, ":", start, "-", end))
# -> only 21 exons with 0 < overlap < 0.9

# add back the 42 exons with no overlap at all
sep_novel_exons <- sep_novel_exons %>%
  add_row(novel_exons[countOverlaps(novel_ex_gr, known_ex_gr) == 0,])


# Export ----

# Save as suppl table S13
sep_novel_exons %>%
  dplyr::rename(position = full_pos) %>%
  arrange(chr, start, end) %>%
  mutate(start = format(start, big.mark = ",", scientific = FALSE, trim = TRUE),
         end = format(end, big.mark = ",", scientific = FALSE, trim = TRUE)) %>%
  write_csv("output/tables/table_S13_novel_exons.csv")
  


