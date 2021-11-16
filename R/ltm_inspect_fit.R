#' @export
obs_proportions <- function(x, ..., nc = 18) {
  if (vctrs::vec_unique_count(x) != nc)
    tb <- table(c(x, 1:nc), ...) - 1L
  else
    tb <- table(x, ...)

  return(tb / sum(tb))
}

#' @export
inspect_probability_plot <- function(observed_proportions, model_probabilities) {

  nc <- length(observed_proportions)
  diff <- observed_proportions - model_probabilities
  tib <- tibble(
    category    = factor(rep(1:nc, 3)),
    probability = c(observed_proportions, model_probabilities, diff),
    group       = factor(rep(c("observed", "model", "diff"), each = nc)),
    panel       = rep(c("raw", "difference"), c(2*nc, nc))
  )

  tibAbline <- tibble(
    intercept = 0,
    slope     = 0,
    panel     = "difference"
  )

  yBreaks <- jaspGraphs::getPrettyAxisBreaks(c(0, tib$probability))
  ggplot(data = tib, aes(x = category, y = probability, group = group, fill = group)) +
    jaspGraphs::geom_abline2(data = tibAbline, mapping = aes(intercept = intercept, slope = slope)) +
    geom_bar(position="dodge", stat="identity", width = .4) +
    # scale_y_continuous(breaks = yBreaks, limits = range(yBreaks)) +
    labs(fill = NULL, x = "Category", y = "Probability") +
    facet_wrap(~panel, scales = "free_y") +
    jaspGraphs::geom_rangeframe() +
    jaspGraphs::themeJaspRaw(legend.position = "right")
}
