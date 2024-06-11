library(tidyverse)
options(readr.show_col_types=F)


rawtbl = read_tsv("tmp/all_sample_coverages_10k.bed", col_names=c("chrom", "start", "end", sprintf("%s_cover", bams),  sprintf("%s_reads", bams))) %>%
    glimpse()

rawtbl %>%
    filter(chrom==zoom_chrom, start >= zoom_start, end <= zoom_end) %>%
    kview()

comparisons = tribble(
    ~comparison, ~sample, ~testctrl,
    "S30GY_01vSA_GN",  "S30GY_01", "test",
    "S30GY_01vSA_GN",  "SA_GN", "control",
    "S40GY_16vSA_GN",  "S40GY_16", "test",
    "S40GY_16vSA_GN",  "SA_GN", "control",
    "S40GY_7vSA_GN",  "S40GY_7", "test",
    "S40GY_7vSA_GN",  "SA_GN", "control",
    "S40GY_9vSA_GN",  "S40GY_9", "test",
    "S40GY_9vSA_GN",  "SA_GN", "control",
)

ggplot(x, aes(x=start, y=log2_fc_norm)) +
    geom_point(aes(alpha=abs(log2_fc_norm)>1.5, colour=log10(pvalue))) +
    geom_segment(aes(x=start, xend=end, y=0, yend=0), colour="red", data=hits) +
    #scale_alpha_manual(values=c(0.1, 1)) +
    #geom_smooth() +
    scale_color_viridis_c() +
    facet_grid(comparison~chrom, space = "free_x", scales="free_x") +
    ylim(c(-3, 3)) +
    theme_classic()
ggsave("testplot.png", width=65, height=8, dpi=300, limitsize=F)

zoom_chrom = "chr29"
zoom_start = 0e6
zoom_end = 5e6

x %>%
    filter(chrom==zoom_chrom, start >= zoom_start, end <= zoom_end) %>%
ggplot(aes(x=start, y=log2_fc_norm)) +
    #geom_vline(xintercept = 22.06e6, colour="red", alpha=0.5) +
    geom_point(aes(alpha=abs(log2_fc_norm)>1.5, colour=log10(pvalue))) +
    geom_segment(aes(x=start, xend=end, y=0, yend=0), colour="red",
                 data=filter(hits, chrom==zoom_chrom, start >= zoom_start, end <= zoom_end)) +
    scale_color_viridis_c() +
    #scale_alpha_manual(values=c(0.1, 1)) +
    facet_grid(comparison~chrom, space = "free_x", scales="free_x") +
    ylim(c(-3, 3)) +
    theme_classic()
ggsave(sprintf("zoomplot_%s-%d-%d.png", zoom_chrom, zoom_start, zoom_end), width=8, height=8, dpi=300, limitsize=F)

zoom_chrom = "chr32"
zoom_start = 0e6
zoom_end = 2e6

x %>%
    filter(chrom==zoom_chrom, start >= zoom_start, end <= zoom_end) %>%
ggplot(aes(x=start, y=log2_fc_norm)) +
    #geom_vline(xintercept = 22.06e6, colour="red", alpha=0.5) +
    geom_point(aes(alpha=abs(log2_fc_norm)>1.5, colour=log10(pvalue))) +
    geom_segment(aes(x=start, xend=end, y=0, yend=0), colour="red",
                 data=filter(hits, chrom==zoom_chrom, start >= zoom_start, end <= zoom_end)) +
    scale_color_viridis_c() +
    #scale_alpha_manual(values=c(0.1, 1)) +
    facet_grid(comparison~chrom, space = "free_x", scales="free_x") +
    ylim(c(-3, 3)) +
    theme_classic()
ggsave(sprintf("zoomplot_%s-%d-%d.png", zoom_chrom, zoom_start, zoom_end), width=8, height=8, dpi=300, limitsize=F)

zoom_chrom = "chr21"
zoom_start = 3.5e6
zoom_end = 4e6


x %>%
    filter(chrom==zoom_chrom, start >= zoom_start, end <= zoom_end) %>%
ggplot(aes(x=start, y=log2_fc_norm)) +
    #geom_vline(xintercept = 22.06e6, colour="red", alpha=0.5) +
    geom_point(aes(alpha=abs(log2_fc_norm)>1.5, colour=log10(pvalue))) +
    geom_segment(aes(x=start, xend=end, y=0, yend=0), colour="red",
                 data=filter(hits, chrom==zoom_chrom, start >= zoom_start, end <= zoom_end)) +
    scale_color_viridis_c() +
    #scale_alpha_manual(values=c(0.1, 1)) +
    facet_grid(comparison~chrom, space = "free_x", scales="free_x") +
    ylim(c(-3, 3)) +
    theme_classic()
ggsave(sprintf("zoomplot_%s-%d-%d.png", zoom_chrom, zoom_start, zoom_end), width=8, height=8, dpi=300, limitsize=F)

zoom_chrom = "chr06"
zoom_start = 35e6
zoom_end = 60e6



x %>%
    filter(chrom==zoom_chrom, start >= zoom_start, end <= zoom_end) %>%
ggplot(aes(x=start, y=log2_fc_norm)) +
    #geom_vline(xintercept = 22.06e6, colour="red", alpha=0.5) +
    geom_point(aes(alpha=abs(log2_fc_norm)>1.5, colour=log10(pvalue))) +
    geom_segment(aes(x=start, xend=end, y=0, yend=0), colour="red",
                 data=filter(hits, chrom==zoom_chrom, start >= zoom_start, end <= zoom_end)) +
    scale_color_viridis_c() +
    #scale_alpha_manual(values=c(0.1, 1)) +
    facet_grid(comparison~chrom, space = "free_x", scales="free_x") +
    ylim(c(-3, 3)) +
    theme_classic()
ggsave(sprintf("zoomplot_%s-%d-%d.png", zoom_chrom, zoom_start, zoom_end), width=8, height=8, dpi=300, limitsize=F)

zoom_chrom = "chr10"
zoom_start = 33500000
zoom_end = 34500000



x %>%
    filter(chrom==zoom_chrom, start >= zoom_start, end <= zoom_end) %>%
ggplot(aes(x=start, y=log2_fc_norm)) +
    #geom_vline(xintercept = 22.06e6, colour="red", alpha=0.5) +
    geom_point(aes(alpha=abs(log2_fc_norm)>1.5, colour=log10(pvalue))) +
    geom_segment(aes(x=start, xend=end, y=0, yend=0), colour="red",
                 data=filter(hits, chrom==zoom_chrom, start >= zoom_start, end <= zoom_end)) +
    scale_color_viridis_c() +
    #scale_alpha_manual(values=c(0.1, 1)) +
    facet_grid(comparison~chrom, space = "free_x", scales="free_x") +
    ylim(c(-3, 3)) +
    theme_classic()
ggsave(sprintf("zoomplot_%s-%d-%d.png", zoom_chrom, zoom_start, zoom_end), width=8, height=8, dpi=300, limitsize=F)


