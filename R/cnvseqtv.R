#' Convert NAs to 0
#' @param x  numeric vector
#' @return A numeric vector with NA/infinite values replaced with 0
na0 = function(x) {ifelse(is.na(x) | !is.finite(x), 0, x)}

#' Calculate the Gearey-Hinkley transform p-value.
#' @param rat  Ratio between test and control
#' @param xmean  Mean coverage of test
#' @param ymean  Mean coverage of control
#' @return Vector of p-values
ghratiop = function(rat, xmean, ymean) {
    # Gearey-Hinkley transform makes this follow the normal distribution
	gh_zscore = (ymean*rat-xmean)/sqrt(ymean*rat^2+xmean)
    purrr::map2_dbl(gh_zscore, gh_zscore<0, \(x, y) stats::pnorm(x, lower.tail=y))
}

#' Overlap adjacent windows to define regions of interest
#' @param starts Vector of start positions
#' @param ends Vector of end positions
#' @return Vector of group identities
overlap_groups = function(starts, ends) {
    irange = IRanges::IRanges(starts, ends)
    unname(as.matrix((IRanges::findOverlaps(irange, IRanges::reduce(irange))))[,2])
}

#' Calculate CNVs
#' @param count_table  A table of windowed counts read directly from the output of samtools depth. See worked example.
#' @param comparisons   A table of comparisons, must have columns
#'                      `"comparison", "sample", "ctrltest"`. Comparison is a
#'                      name for this comparison, sample is which sample name
#'                      (must match bam name in `count_table`), and ctrltest is
#'                      a string either `"control"` or `"test"`, denoting if
#'                      this row describes a test or control sample. Therefore,
#'                      there should be two rows for each comparison, one for
#'                      the test and one for the control.
#' @export
#' @return A table of windows in which the coverage between the control and test sample(s) are compared.
calculate_cnvs = function(count_table, comparisons) {
    count_table %>%
	dplyr::rename(chrom=`#chrom`) %>%
        dplyr::select(c(chrom, start, end, dplyr::ends_with("_count"))) %>%
        tidyr::pivot_longer(-c(chrom, start, end), names_to = "sample",
                            names_transform = function(x) sub(".*/([^.]+).bam_count", "\\1", x, perl=T), values_to = "nreads") %>%
        dplyr::group_by(sample) %>%
        dplyr::mutate(nreads=nreads+1, rpm = nreads/sum(nreads)*1e6, meanreads=mean(nreads)) %>%
        dplyr::ungroup() %>%
        dplyr::inner_join(comparisons, by="sample", relationship = "many-to-many") %>%
        tidyr::pivot_longer(c(nreads, rpm, meanreads), names_to = "variable", values_to = "value") %>%
        dplyr::select(-sample) %>%
        tidyr::pivot_wider(names_from=c("testctrl", "variable"), values_from="value") %>%
        dplyr::group_by(comparison) %>%
        dplyr::mutate(
            log2_fc=na0(log2(test_rpm/control_rpm)),
            log2_fc_norm=log2_fc-mean(na0(log2(test_rpm/control_rpm))),
            pvalue = ghratiop(test_nreads/control_nreads, test_meanreads, control_meanreads)
        )
}

#' Extract ROIs from CNV calculation
#' @param cnv_table Windowed comparisons between control and test sample(s), output directly from `calculate_cnvs`.
#' @param log2fc_thresh log2 fold-change threshold value, keep windows with abs(log2_fc) > X
#' @param pval_thresh P-value threshold for signficant difference between test & control, windows with p < X are kept.
#' @export
#' @return ROIs of merged genome windows passing these filtering thresholds
cnv_rois = function(cnv_table, log2fc_thresh=1.5, pval_thresh=1e-12) {
    cnv_table %>%
    dplyr::filter(abs(log2_fc_norm) > log2fc_thresh & pvalue < pval_thresh) %>%
    dplyr::group_by(comparison, chrom) %>%
    dplyr::mutate(hit_group = overlap_groups(start, end)) %>%
    dplyr::group_by(comparison, chrom, hit_group) %>%
    dplyr::summarise(
        start=min(start),
        end=max(end),
        n_windows=dplyr::n(),
        dplyr::across(dplyr::ends_with("nreads"), sum),
        dplyr::across(c(dplyr::ends_with("meanreads"), dplyr::ends_with("rpm"), dplyr::starts_with("log2_fc"), pvalue), mean),
        .groups="drop"
    )
}

#' Create a ggplot figure.
#' @param cnv_table Windowed comparisons between control and test sample(s), output directly from `calculate_cnvs`.
#' @param hits Hit ROIs from cnv_rois
#' @param log2fc_thresh log2 fold-change threshold value, keep windows with abs(log2_fc) > X
#' @param pval_thresh P-value threshold for signficant difference between test & control, windows with p < X are kept.
#' @param zoom_chrom  Zoom to chrom
#' @param zoom_start  Start position for zoom-in
#' @param zoom_end  End position for zoom-in
#' @export
#' @return ROIs of merged genome windows passing these filtering thresholds
cnv_plot = function(cnv_table, hits, log2fc_thresh=1.5, pval_thresh=1e-12, zoom_chrom=NULL, zoom_start=NULL, zoom_end=NULL) {
    pdat = cnv_table
    if (!is.null(zoom_chrom)) {
        if (is.null(zoom_start)) zoom_start=0
        if (is.null(zoom_end)) zoom_end=1e100
        pdat = pdat %>% dplyr::filter(chrom==zoom_chrom, start >= zoom_start, end <= zoom_end)
        hits = dplyr::filter(hits, chrom==zoom_chrom, start >= zoom_start, end <= zoom_end)
    }

    ggplot2::ggplot(pdat, ggplot2::aes(x=start, y=log2_fc_norm)) +
        ggplot2::geom_point(ggplot2::aes(alpha=abs(log2_fc_norm)>log2fc_thresh & pvalue < pval_thresh, colour=log10(pvalue))) +
        ggplot2::geom_segment(ggplot2::aes(x=start, xend=end, y=0, yend=0), colour="red", data=hits) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::facet_grid(comparison~chrom, space = "free_x", scales="free_x") +
	ggplot2::labs(alpha="hit") +
        ggplot2::theme_classic()
}
