#' @importFrom assertthat assert_that
#' @importFrom tibble tibble
#' @importFrom ggplot2 aes ggplot scale_color_manual geom_abline geom_point ggtitle facet_wrap theme_bw
#' @importFrom jaspGraphs geom_rangeframe themeJaspRaw

#' @export
pplot <- function(x, y, ...) {
  tibble(x = x, y = y) |>
    ggplot(mapping = ggplot2::aes(x = x, y = y)) +
    geom_abline() +
    geom_point() +
    ggtitle(sprintf("cor: %.3f", cor(x, y))) +
    geom_rangeframe() +
    themeJaspRaw()
}

#' @export
scatterplot_retrieval <- function(tib, mapping = NULL, facets = ~parameter, scales = "free") {

  assert_that(nrow(tib) > 0L)

  if (is.null(mapping)) {
    if ("rater_group" %in% names(tib) && !all(is.na(tib$rater_group)) && !all(tib$rater_group == 1L)) {
      mapping <- ggplot2::aes(x = true_value, y = estimate, color = rater_group)
      no_rater_groups <- max(as.integer(tib$rater_group), na.rm = TRUE)
      color_scale <- scale_color_manual(values = scales::hue_pal()(no_rater_groups), breaks = seq_len(no_rater_groups))
    } else {
      mapping <- ggplot2::aes(x = true_value, y = estimate)
      color_scale <- NULL
    }
  }

  if (length(vars <- all.vars(facets)) == 2L) {
    ncol <- vctrs::vec_unique_count(tib[[vars[2L]]])
  } else {
    ncol <- NULL
  }

  ggplot(data = tib, mapping = mapping) +
    geom_abline() +
    geom_point() +
    facet_wrap(facets = facets, scales = scales, ncol = ncol) +
    color_scale +
    theme_bw()
}
