# PLUD: pathogen load upon death, per group
# Input: the per-fly table saved by CFU Drop Counter (Load step), one row per fly, with the log10 load.
# Screens each group for flies that did not die of the infection (Rosner test) and compares groups
# with a linear mixed model of log10 PLUD, as in the protocol ("Analysis Bacterial Load Upon Death").
# Runs as a normal Shiny app, or in the browser with shinylive (no server).

library(shiny)
library(bslib)
library(ggplot2)
library(lme4)
library(EnvStats)
library(patchwork)
if (FALSE) { library(svglite); library(jsonlite) }   # declares them to shinylive, which bundles the packages it finds in library() calls
has_svglite <- requireNamespace("svglite", quietly = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a

APP_VERSION <- "1.0"
EXAMPLE_FILE <- "example_PLUD_RAL821.csv"
PALETTE <- c("#1f5fbf", "#e8850c", "#2a9d6f", "#c0398a", "#7a5cc4", "#8c6d1f", "#3a3a3a", "#d1453b")
# Female orange and male blue in every app of the protocol, whatever the order; other groups take the palette in order
group_colours <- function(lv, custom = NULL) {
  cols <- setNames(rep(NA_character_, length(lv)), lv)
  low <- tolower(lv)
  cols[low %in% c("female", "females", "f")] <- "#e8850c"
  cols[low %in% c("male", "males", "m")] <- "#1f5fbf"
  free <- setdiff(PALETTE, cols)
  cols[is.na(cols)] <- rep(free, length.out = sum(is.na(cols)))
  hit <- intersect(names(custom), lv); cols[hit] <- unlist(custom[hit])      # colours chosen by the user
  cols
}
# a difference in bacteria per fly, as a plotmath label: +5.0 x 10^6
bact_label <- function(v) {
  out <- rep("", length(v)); ok <- is.finite(v) & v != 0
  e <- floor(log10(abs(v[ok]))); m <- abs(v[ok]) / 10^e
  out[ok] <- sprintf('"%s%.1f" %%*%% 10^%d', ifelse(v[ok] < 0, "\u2212", "+"), m, e)
  out
}
fmt_bact <- function(v) ifelse(is.finite(v), formatC(v, format = "e", digits = 1), "")

# ------------------------------------------------------------------ reading
read_any <- function(path) {
  first <- readLines(path, n = 1, warn = FALSE)
  seps <- c(";" = ";", "," = ",", "\t" = "\t")
  counts <- sapply(seps, function(s) lengths(regmatches(first, gregexpr(s, first, fixed = TRUE))))
  sep <- if (max(counts) > 0) seps[[which.max(counts)]] else ""     # "" = any white space
  read.table(path, header = TRUE, sep = sep, quote = "\"'", dec = ".", stringsAsFactors = FALSE,
             check.names = FALSE, comment.char = "", strip.white = TRUE, na.strings = c("NA", ""))
}
guess <- function(cols, candidates) { hit <- candidates[candidates %in% cols]; if (length(hit)) hit[1] else "" }
num <- function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(sub(",", ".", trimws(as.character(x)), fixed = TRUE)))   # accepts decimal commas
chr <- function(x) { x <- trimws(as.character(x)); x[is.na(x) | x == ""] <- "NA"; x }    # missing values become a group "NA"
# default order: alphabetical, with a wild-type or control group first (it then becomes the reference)
default_order <- function(values) {
  lv <- sort(unique(chr(values)))
  first <- grepl("^(wild|wt([^a-z]|$)|control|ctrl)", tolower(lv))
  c(lv[first], lv[!first])
}
order_levels <- function(values, text) {
  all_lv <- sort(unique(chr(values)))
  given <- trimws(strsplit(text %||% "", ",")[[1]])
  c(given[given %in% all_lv], setdiff(all_lv, given))
}

# ------------------------------------------------------------------ screen
# Rosner test (generalised ESD) within each group, at most k values, two-sided.
screen_groups <- function(pf, k, alpha, low_only) {
  pf$Flagged <- FALSE; pf$Removed <- FALSE; pf$Direction <- ""
  cells <- unique(pf$.cell[order(pf$.panel, pf$.x, pf$.col)])
  summary <- lapply(cells, function(g) {
    i <- which(pf$.cell == g)
    res <- data.frame(Group = g, Flies = length(i), Flagged_low = 0L, Flagged_high = 0L, Removed = 0L, Note = "")
    if (length(i) < max(3, k + 2)) { res$Note <- "too few flies to test"; return(res) }
    rt <- tryCatch(EnvStats::rosnerTest(pf$.y[i], k = k, alpha = alpha), error = function(e) NULL)
    if (is.null(rt)) { res$Note <- "test failed"; return(res) }
    st <- rt$all.stats
    idx <- i[st$Obs.Num[st$Outlier]]
    if (length(idx)) {
      low <- pf$.y[idx] < mean(pf$.y[i])
      pf$Flagged[idx] <<- TRUE
      pf$Direction[idx] <<- ifelse(low, "low", "high")
      pf$Removed[idx] <<- if (low_only) low else TRUE
      res$Flagged_low <- sum(low); res$Flagged_high <- sum(!low); res$Removed <- if (low_only) sum(low) else length(idx)
    }
    if (length(i) < 25) res$Note <- "fewer than 25 flies: the test is less reliable"
    res
  })
  list(pf = pf, summary = do.call(rbind, summary))
}

# ------------------------------------------------------------------ model
# One analysis: log10 PLUD ~ compared factor (+ within factor), box and day as random effects.
#   joint = FALSE: .y ~ .cmp, on the flies given (one within-group level at a time)
#   joint = TRUE : .y ~ .w + .w:.cmp, on all within-group levels at once, sharing the residual
#                  variance and the random effects; the comparison is made within each level of .w
fit_analysis <- function(x, has_box, has_day, joint) {
  x$.cmp <- droplevels(x$.cmp)
  if (joint) x$.w <- droplevels(x$.w)
  re <- c(if (has_box && length(unique(x$.box)) > 1) "(1 | .box)", if (has_day && length(unique(x$.day)) > 1) "(1 | .day)")
  mixed <- length(re) > 0
  fx1 <- if (joint) ".w + .w:.cmp" else ".cmp"; fx0 <- if (joint) ".w" else "1"
  fxm <- if (joint) "0 + .w:.cmp" else "0 + .cmp"
  f <- function(fx) as.formula(paste(".y ~", paste(c(fx, re), collapse = " + ")))
  ctrl <- lmerControl(check.conv.singular = "ignore")
  if (mixed) {
    m1 <- lmer(f(fx1), data = x, control = ctrl)
    chisq <- as.numeric(2 * (logLik(lmer(f(fx1), data = x, REML = FALSE, control = ctrl)) -
                               logLik(lmer(f(fx0), data = x, REML = FALSE, control = ctrl))))
    mm <- lmer(f(fxm), data = x, control = ctrl)
    b1 <- fixef(m1); bm <- fixef(mm); vm <- as.matrix(vcov(mm))
    vc <- as.data.frame(VarCorr(m1)); re_sd <- setNames(vc$sdcor[vc$grp != "Residual"], sub("^\\.", "", vc$grp[vc$grp != "Residual"]))
    singular <- isSingular(m1)
  } else {
    m1 <- lm(f(fx1), data = x)
    chisq <- as.numeric(2 * (logLik(m1) - logLik(lm(f(fx0), data = x))))
    mm <- lm(f(fxm), data = x)
    b1 <- coef(m1); bm <- coef(mm); vm <- vcov(mm); re_sd <- numeric(0); singular <- FALSE
  }
  # comparison terms: .cmpX (separate) or .wW:.cmpX (joint); NA terms (empty cells) are dropped
  terms <- names(b1)[grepl(".cmp", names(b1), fixed = TRUE) & !is.na(b1)]
  ci <- if (mixed) tryCatch(confint(m1, parm = terms, method = "profile", quiet = TRUE), error = function(e) NULL)
        else confint(m1)[terms, , drop = FALSE]
  ci_kind <- if (mixed) "profile likelihood" else "t-based"
  if (is.null(ci)) { se <- sqrt(diag(as.matrix(vcov(m1))))[terms]; ci <- cbind(b1[terms] - 1.96 * se, b1[terms] + 1.96 * se); ci_kind <- "Wald (profile failed)" }
  ref <- levels(x$.cmp)[1]
  parse_term <- function(t) {
    if (joint) { p <- strsplit(t, ":", fixed = TRUE)[[1]]; c(sub("^\\.w", "", p[1]), sub("^\\.cmp", "", p[2])) }
    else c(NA, sub("^\\.cmp", "", t))
  }
  tl <- t(vapply(terms, parse_term, c("", "")))
  # group means (fixed effects only)
  mn <- names(bm)[!is.na(bm)]
  ml <- t(vapply(mn, function(t) if (joint) { p <- strsplit(t, ":", fixed = TRUE)[[1]]; c(sub("^\\.w", "", p[1]), sub("^\\.cmp", "", p[2])) }
                                   else c(NA, sub("^\\.cmp", "", t)), c("", "")))
  mse <- sqrt(diag(vm))[mn]
  nfl <- if (joint) as.vector(table(paste(x$.w, x$.cmp))[paste(ml[, 1], ml[, 2])]) else as.vector(table(x$.cmp)[ml[, 2]])
  means <- data.frame(Within = ml[, 1], Group = ml[, 2], Flies = nfl, Mean = unname(bm[mn]),
                      Lower = unname(bm[mn] - 1.96 * mse), Upper = unname(bm[mn] + 1.96 * mse))
  ref_mean <- function(w) means$Mean[means$Group == ref & (is.na(means$Within) | means$Within %in% w)][1]
  diffs <- data.frame(Within = tl[, 1], Group = tl[, 2], Difference = unname(b1[terms]),
                      Lower = unname(ci[, 1]), Upper = unname(ci[, 2]))
  in_bact <- function(d) vapply(seq_len(nrow(diffs)), function(i) 10^ref_mean(diffs$Within[i]) * (10^d[i] - 1), 0)
  diffs$Bacteria_per_fly <- in_bact(diffs$Difference)
  diffs$Bacteria_lower <- in_bact(diffs$Lower); diffs$Bacteria_upper <- in_bact(diffs$Upper)
  cell <- if (joint) paste(x$.w, x$.cmp) else as.character(x$.cmp)
  total_sd <- sqrt(sum(tapply(x$.y, cell, function(v) sum((v - mean(v))^2))) / (nrow(x) - length(unique(cell))))
  n_h <- length(unique(cell)) / sum(1 / as.vector(table(cell)))
  list(means = means, diffs = diffs, ref = ref, joint = joint,
       lrt = data.frame(Chisq = chisq, df = length(terms), p = pchisq(chisq, length(terms), lower.tail = FALSE)),
       re_sd = re_sd, resid_sd = sigma(m1), total_sd = total_sd, n_h = n_h,
       detectable = 2.8 * total_sd * sqrt(2 / n_h), mixed = mixed, re = re, ci_kind = ci_kind, singular = singular)
}

# all analyses: one per panel, and within a panel one per level of the other factor unless joint
analyse <- function(pf, s) {
  out <- list()
  for (pn in levels(pf$.panel)) {
    x <- pf[pf$.panel == pn, ]
    if (!nrow(x)) next
    lab <- function(...) { v <- paste(c(if (nzchar(s$panel)) pn, ...), collapse = " / "); if (nzchar(v)) v else "All" }
    if (!nzchar(s$colour)) {
      out[[lab()]] <- tryCatch(fit_analysis(x, nzchar(s$box), nzchar(s$day), FALSE), error = function(e) e)
    } else {
      x$.w <- if (s$compare == "colour") x$.x else x$.col
      if (s$joint) {
        out[[lab()]] <- tryCatch(fit_analysis(x, nzchar(s$box), nzchar(s$day), TRUE), error = function(e) e)
      } else for (w in levels(droplevels(x$.w))) {
        out[[lab(w)]] <- tryCatch(fit_analysis(x[x$.w == w, ], nzchar(s$box), nzchar(s$day), FALSE), error = function(e) e)
      }
    }
  }
  out
}

# ------------------------------------------------------------------ plots
theme_plud <- theme_classic(base_size = 12) +
  theme(axis.line = element_line(linewidth = 0.6), axis.ticks = element_line(linewidth = 0.6),
        axis.text = element_text(colour = "black"), legend.position = "top", legend.title = element_blank(),
        strip.background = element_rect(fill = NA, colour = "black"))

plot_plud <- function(pf, s) {
  kept <- pf[!pf$Removed, ]
  shapes <- scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 4), guide = "none")
  if (nzchar(s$colour)) {
    cols <- group_colours(levels(pf$.col), s$cols)
    dg <- position_dodge(width = 0.8)
    mm <- aggregate(.y ~ .x + .col + .panel, data = kept, FUN = mean)
    p <- ggplot(pf, aes(.x, .y, fill = .col)) +
      geom_violin(data = kept, position = dg, colour = "grey30", linewidth = 0.3, alpha = 0.55, scale = "width", width = 0.75) +
      geom_point(aes(shape = Removed, group = .col), position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.8, seed = 1),
                 colour = "black", size = 1.4, stroke = 0.4) +
      geom_point(data = mm, aes(group = .col), position = dg, shape = 95, size = 9, colour = "black", show.legend = FALSE) +
      scale_fill_manual(values = cols) + shapes
  } else {
    mm <- aggregate(.y ~ .x + .panel, data = kept, FUN = mean)
    p <- ggplot(pf, aes(.x, .y, fill = .x)) +
      geom_violin(data = kept, colour = "grey30", linewidth = 0.3, alpha = 0.55, scale = "width", width = 0.75, show.legend = FALSE) +
      geom_point(aes(shape = Removed), position = position_jitter(width = 0.12, height = 0, seed = 1),
                 colour = "black", size = 1.4, stroke = 0.4, show.legend = FALSE) +
      geom_point(data = mm, shape = 95, size = 9, colour = "black", show.legend = FALSE) +
      scale_fill_manual(values = group_colours(levels(pf$.x), s$cols), guide = "none") + shapes
  }
  top <- max(pf$.y) + 0.12 * diff(range(pf$.y))
  if (nzchar(s$colour)) {
    ns <- aggregate(.y ~ .x + .col + .panel, data = pf, FUN = length)
    p <- p + geom_text(data = ns, aes(.x, top, label = .y, group = .col), position = position_dodge(width = 0.8), inherit.aes = FALSE, size = 3.1) +
      annotate("text", x = -Inf, y = top, label = "N", hjust = -0.3, size = 3.1, fontface = "italic")
  } else {
    ns <- aggregate(.y ~ .x + .panel, data = pf, FUN = length)
    p <- p + geom_text(data = ns, aes(.x, top, label = paste0("(N=", .y, ")")), inherit.aes = FALSE, size = 3.3)
  }
  p <- p +
    labs(x = s$xlab, y = s$ylab) + theme_plud
  if (nzchar(s$panel)) p + facet_wrap(~.panel, scales = "free_x") else p + scale_x_discrete(drop = FALSE)
}

plot_diffs <- function(res, s) {
  rows <- do.call(rbind, lapply(names(res), function(nm) {
    f <- res[[nm]]; if (inherits(f, "error")) return(NULL)
    d <- f$diffs
    refs <- data.frame(Within = unique(d$Within), Group = f$ref, Difference = 0, Lower = 0, Upper = 0)
    refs$Bacteria_per_fly <- NA
    d <- rbind(refs, d[, c("Within", "Group", "Difference", "Lower", "Upper", "Bacteria_per_fly")])
    d$Panel <- ifelse(is.na(d$Within), nm, paste(nm, d$Within, sep = if (nchar(nm)) " / " else ""))
    d$Panel <- sub("^All / |^ / ", "", d$Panel); d
  }))
  if (is.null(rows) || !nrow(rows)) return(NULL)
  lv <- if (s$compare == "colour") s$col_levels else s$x_levels
  rows$Group <- factor(rows$Group, levels = intersect(lv, unique(rows$Group)))
  p <- ggplot(rows, aes(Group, Difference)) +
    geom_hline(yintercept = 0, colour = "red", linetype = 4) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.25, linewidth = 0.5) +
    geom_point(aes(colour = Group), size = 2.5, show.legend = FALSE) +
    geom_text(aes(y = Upper, label = bact_label(Bacteria_per_fly)), parse = TRUE, vjust = -0.6, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_colour_manual(values = group_colours(levels(rows$Group), s$cols)) +
    labs(x = NULL, y = paste0("Difference in log10 PLUD (95% CI)\nrelative to ", s$ref),
         caption = "Numbers: the difference in bacteria per fly") + theme_plud +
    theme(plot.caption = element_text(size = 9, colour = "grey30"))
  if (sum(nchar(levels(rows$Group))) > 24) p <- p + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  if (length(unique(rows$Panel)) > 1) {
    rows$Panel <- factor(rows$Panel, levels = unique(rows$Panel))
    p <- p %+% rows + facet_wrap(~Panel)
  }
  p
}

# ------------------------------------------------------------------ reports
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "f", digits = 3)))

methods_text <- function(res, s, scr, n_flies) {
  cells <- paste(c(if (nzchar(s$panel)) s$panel, s$x, if (nzchar(s$colour)) s$colour), collapse = " x ")
  filt <- if (nzchar(s$filter_col) && length(s$filter_vals)) sprintf("%s %s", s$filter_col, paste(s$filter_vals, collapse = ", "))
  scr_txt <- if (s$screen) {
    sm <- scr$summary
    sprintf("Within each group (%s), log10 PLUD was screened for flies that did not die of the infection with a Rosner test (generalized extreme studentized deviate; EnvStats), at a %s%% threshold and at most k = %d values per group; %s. The screen removed %d of %d flies (%s). Removed values are interpreted as deaths from other causes, an inference that the test does not establish. ",
            cells, 100 * s$alpha, s$k, if (s$low_only) "only values below the group mean were removed" else "flagged values on either side were removed",
            sum(sm$Removed), sum(sm$Flies), paste(sprintf("%s: %d", sm$Group, sm$Removed), collapse = ", "))
  } else "No screen for outlying values was applied. "
  first <- Filter(function(f) !inherits(f, "error"), res)[[1]]
  re_names <- c(box = s$box, day = s$day)[gsub("\\(1 \\| \\.|\\)", "", first$re %||% character(0))]
  cmp <- if (s$compare == "colour") s$colour else s$x
  within <- if (nzchar(s$colour)) (if (s$compare == "colour") s$x else s$colour) else ""
  model <- sprintf("Log10 PLUD was compared between levels of %s%s with a linear mixed model%s (lme4)%s; the effect was tested by a likelihood ratio test on maximum likelihood fits, and differences are given relative to %s with %s 95%% confidence intervals. ",
                   cmp, if (nzchar(within)) paste0(" within each level of ", within) else "",
                   if (length(re_names)) paste0(" with ", paste(re_names, collapse = " and "), " as random effect", if (length(re_names) > 1) "s" else "") else "",
                   if (nzchar(within)) (if (s$joint) sprintf(", fitted to all levels of %s at once (%s + %s:%s)", within, within, within, cmp) else sprintf(", fitted separately for each level of %s", within)) else "",
                   s$ref, first$ci_kind)
  res_txt <- sapply(names(res), function(nm) {
    f <- res[[nm]]
    if (inherits(f, "error")) return(sprintf("%s: the model could not be fitted (%s).", nm, conditionMessage(f)))
    dd <- f$diffs
    paste0(if (nzchar(nm) && nm != "All") paste0(nm, ". ") else "",
           paste(sprintf("%s%s: %+.2f log10 (95%% CI %.2f to %.2f), that is %s bacteria per fly (%s to %s)", ifelse(is.na(dd$Within), "", paste0(dd$Within, ", ")),
                         dd$Group, dd$Difference, dd$Lower, dd$Upper, fmt_bact(dd$Bacteria_per_fly), fmt_bact(dd$Bacteria_lower), fmt_bact(dd$Bacteria_upper)), collapse = "; "),
           sprintf("; likelihood ratio test, chi2 = %.1f, df = %d, p = %s. Within-group SD of log10 PLUD = %.2f; with %.0f flies per group, a difference of %.2f log10 would be detected with 80%% power at 5%%.",
                   f$lrt$Chisq, f$lrt$df, fmt_p(f$lrt$p), f$total_sd, f$n_h, f$detectable))
  })
  paste0(if (length(filt)) paste0("Flies of ", filt, " were left out of the analysis. ") else "",
         scr_txt, model, "Analysis made with the PLUD app version ", APP_VERSION, " on ", n_flies, " flies.\n\n",
         paste(res_txt, collapse = "\n"))
}

r_code <- function(s) {
  re <- c(if (nzchar(s$box)) "(1 | box)", if (nzchar(s$day)) "(1 | day)")
  rhs <- function(fx) paste(c(fx, re), collapse = " + ")
  cmp <- if (s$compare == "colour") "colour" else "x"
  w <- if (nzchar(s$colour)) (if (s$compare == "colour") "x" else "colour") else ""
  fx1 <- if (nzchar(w) && s$joint) paste0(w, " + ", w, ":cmp") else "cmp"
  fx0 <- if (nzchar(w) && s$joint) w else "1"
  fun <- if (length(re)) "lmer" else "lm"
  paste0(
'# Reproduces the PLUD app from its table of flies (Flies tab, "PLUD_flies.csv")
library(lme4)
library(EnvStats)

pf <- read.csv("PLUD_flies.csv")
pf$x <- factor(pf$x, levels = c(', paste0('"', s$x_levels, '"', collapse = ", "), '))
', if (nzchar(s$colour)) paste0('pf$colour <- factor(pf$colour, levels = c(', paste0('"', s$col_levels, '"', collapse = ", "), '))\n') else '',
'pf$group <- interaction(', paste(c(if (nzchar(s$panel)) "pf$panel", "pf$x", if (nzchar(s$colour)) "pf$colour"), collapse = ", "), ', drop = TRUE)

# Rosner screen within each group, k = ', s$k, ', threshold ', s$alpha, '
pf$removed <- FALSE
for (g in levels(pf$group)) {
  i <- which(pf$group == g)
  if (length(i) < ', max(3, s$k + 2), ') next
  st <- rosnerTest(pf$log10_load[i], k = ', s$k, ', alpha = ', s$alpha, ')$all.stats
  out <- i[st$Obs.Num[st$Outlier]]
  ', if (s$low_only) 'pf$removed[out] <- pf$log10_load[out] < mean(pf$log10_load[i])   # low values only' else 'pf$removed[out] <- TRUE', '
}
table(pf$group, pf$removed)

d <- subset(pf, ', if (s$screen) '!removed' else 'TRUE', ')
d$cmp <- relevel(factor(d$', cmp, '), ref = "', s$ref, '")
',
if (nzchar(s$panel)) '# One analysis per panel; repeat for each:\n# d <- subset(d, panel == "...")\n' else '',
if (nzchar(w) && !s$joint) paste0('# One analysis per level of ', w, '; repeat for each:\n# d <- subset(d, ', w, ' == "...")\n') else '',
if (length(re)) paste0('
m1   <- lmer(log10_load ~ ', rhs(fx1), ', data = d)
m1ml <- lmer(log10_load ~ ', rhs(fx1), ', data = d, REML = FALSE)
m0ml <- lmer(log10_load ~ ', rhs(fx0), ', data = d, REML = FALSE)
anova(m0ml, m1ml)                                    # likelihood ratio test
fixef(m1)                                            # differences relative to ', s$ref, '
confint(m1, parm = "beta_", method = "profile")      # with profile 95% CI
sigma(m1)                                            # residual SD
') else paste0('
m1 <- lm(log10_load ~ ', fx1, ', data = d)
m0 <- lm(log10_load ~ ', fx0, ', data = d)
anova(m0, m1)
confint(m1)
'))
}

# ------------------------------------------------------------------ downloads
# Files are made in R and handed to the browser as a data URL, which works under shinylive.
dl_button <- function(id, label) actionButton(id, label, icon = icon("download"), class = "btn-outline-secondary btn-sm me-1")
DL_JS <- "
Shiny.addCustomMessageHandler('send-file', function(m) {
  var a = document.createElement('a');
  a.href = 'data:' + m.mime + ';base64,' + m.data;
  a.download = m.name;
  document.body.appendChild(a); a.click(); a.remove();
});"
send_file <- function(session, name, mime, write_fun) {
  f <- tempfile(fileext = paste0(".", tools::file_ext(name)))
  write_fun(f)
  data <- jsonlite::base64_enc(readBin(f, "raw", file.info(f)$size))
  session$sendCustomMessage("send-file", list(name = name, mime = mime, data = gsub("\n", "", data)))
}

# ------------------------------------------------------------------ UI
ui <- page_sidebar(
  title = "PLUD: pathogen load upon death",
  tags$head(tags$script(HTML(DL_JS)), tags$style(HTML("select.shiny-input-select:not(.selectized) { appearance: auto; -webkit-appearance: menulist; }"))),
  sidebar = sidebar(width = 350,
    accordion(id = "steps", open = c("Data", "Comparison"),
      accordion_panel("Data", icon = icon("table"),
        fileInput("file", NULL, accept = c(".csv", ".txt", ".tsv"), placeholder = "Per-fly table from CFU Drop Counter"),
        actionLink("use_example", "Try the example (RAL821, 36 females and 36 males)"),
        div(class = "small text-muted mt-1", textOutput("data_info"))),
      accordion_panel("Comparison", icon = icon("scale-balanced"),
        selectInput("x", "Groups to compare (x axis)", NULL),
        selectInput("ref", "Which group is the least susceptible?", NULL),
        div(class = "small text-muted", style = "margin-top:-0.6rem; margin-bottom:0.8rem",
            "The group with the lowest mortality in the survival assay. It is the reference: as in the protocol, mortality and PLUD then read in one direction."),
        textInput("x_order", "Order on the axis (comma-separated)", ""),
        selectInput("box", "Box (random effect)", NULL, selectize = FALSE),
        conditionalPanel("input.box != ''",
          checkboxInput("box_restart", "Box labels restart in each group (box 1 of one line is not box 1 of another)", FALSE)),
        selectInput("day", "Day of injection (random effect, optional)", NULL, selectize = FALSE)),
      accordion_panel("Second factor and panels", icon = icon("layer-group"),
        p(class = "small text-muted", "Optional: for designs with a second factor, such as age within each sex, and a third one shown as panels."),
        selectInput("colour", "Second factor, by colour", NULL, selectize = FALSE),
        conditionalPanel("input.colour != ''",
          textInput("col_order", "Order of the colours (comma-separated)", ""),
          radioButtons("compare", "Compare", c("the groups on the x axis, within each colour" = "x", "the colours, within each x group" = "colour")),
          radioButtons("joint", "Fit", c("one model for all groups (shared residual variance and boxes)" = "joint", "a separate model for each group" = "separate"))),
        selectInput("panel", "Panels (analysed separately)", NULL, selectize = FALSE)),
      accordion_panel("Leave out groups", icon = icon("filter"),
        p(class = "small text-muted", "Optional: leave out whole groups, such as a genotype or a sex. Single flies cannot be removed here; the screen does that, as described in the protocol."),
        selectInput("filter_col", "Column", NULL, selectize = FALSE),
        selectizeInput("filter_vals", "Groups to leave out", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "none", plugins = list("remove_button")))),
      accordion_panel("Colours", icon = icon("palette"),
        p(class = "small text-muted", "Click a colour to change it. Females are orange and males blue by default, as in the survival app."),
        uiOutput("colour_pickers"),
        actionLink("reset_cols", "Back to the default colours")),
      accordion_panel("Columns and labels", icon = icon("table-columns"),
        p(class = "small text-muted", "Found automatically in a table from CFU Drop Counter."),
        selectInput("load_col", "Load per fly (log scale)", NULL),
        radioButtons("load_base", NULL, c("log10" = "10", "log2" = "2"), inline = TRUE),
        selectInput("status_col", "Detection status", NULL, selectize = FALSE),
        p(class = "small text-muted", "CFU Drop Counter marks each fly \"Detected\" when it has a load. Flies with another status (no colony, or every drop too crowded to count) have no load: they are left out and counted in a message."),
        selectInput("id", "Fly identifier (optional, for the tables)", NULL, selectize = FALSE),
        textInput("xlab", "x-axis label", ""),
        textInput("ylab", "y-axis label", "log10 PLUD (bacteria per fly)"))
    )
  ),
  uiOutput("messages"),
  navset_card_underline(
    nav_panel("Plots",
      plotOutput("plot", height = "480px"),
      div(dl_button("dl_pdf", "PDF"), dl_button("dl_png", "PNG"), if (has_svglite) dl_button("dl_svg", "SVG"))),
    nav_panel("Results", uiOutput("results_ui"), dl_button("dl_results", "Results (CSV)")),
    nav_panel("Screen",
      p(class = "small text-muted", "Each group is tested separately. Removed flies are shown as crosses on the plot."),
      tableOutput("screen_table"), h6("Flagged flies"), tableOutput("flagged_table")),
    nav_panel("Flies", tableOutput("pf_table"), dl_button("dl_pf", "Table of flies (CSV)")),
    nav_panel("Methods text", verbatimTextOutput("methods"), dl_button("dl_methods", "Text")),
    nav_panel("R code", p(class = "small text-muted", "Runs on the table downloaded from the Flies tab."),
              verbatimTextOutput("rcode"), dl_button("dl_code", "R script"))
  )
)

# ------------------------------------------------------------------ server
server <- function(input, output, session) {
  rv <- reactiveValues(d = NULL, name = NULL)
  observeEvent(input$file, {
    d <- tryCatch(read_any(input$file$datapath), error = function(e) NULL)
    if (is.null(d)) { showNotification("This file could not be read as a table.", type = "error"); return() }
    rv$d <- d; rv$name <- input$file$name
  })
  observeEvent(input$use_example, { rv$d <- read_any(EXAMPLE_FILE); rv$name <- EXAMPLE_FILE })

  observeEvent(rv$d, {
    cols <- names(rv$d); none <- c("(none)" = "")
    group_cols <- cols[vapply(rv$d, function(v) { n <- length(unique(v)); n >= 2 && n <= 30 }, TRUE)]   # columns of groups, not of flies
    load <- guess(cols, c("log10_CFU_per_fly", "log10_load", "log2_CFU_per_fly", "log10_PLUD", "Bacteria.per.fly_log10", "log2b", "log2_load"))
    updateSelectInput(session, "id", choices = c("(none: row number)" = "", cols), selected = guess(cols, c("ID", "Fly_ID", "Sample_ID", "sample_ID", "sampleID")))
    updateSelectInput(session, "load_col", choices = cols, selected = load)
    updateRadioButtons(session, "load_base", selected = if (grepl("log2", load)) "2" else "10")
    updateSelectInput(session, "status_col", choices = c(none, cols), selected = guess(cols, c("Status")))
    varying <- cols[vapply(rv$d, function(v) length(unique(v)) > 1, TRUE)]     # a column with a single value cannot be the x axis
    updateSelectInput(session, "x", choices = group_cols, selected = guess(varying, c("Genotype", "Line", "Treatment", "Sex", "Age_at_infection")))
    updateSelectInput(session, "colour", choices = c(none, group_cols), selected = "")
    updateSelectInput(session, "panel", choices = c(none, group_cols), selected = "")
    updateSelectInput(session, "box", choices = c(none, cols), selected = guess(cols, c("Box_ID", "Box", "box")))
    updateSelectInput(session, "day", choices = c(none, cols), selected = guess(cols, c("Day_of_injection", "Day", "Experiment", "Experiment_ID")))
    updateSelectInput(session, "filter_col", choices = c(none, group_cols), selected = "")
  })

  dat <- reactive({
    req(rv$d); d <- rv$d
    if (nzchar(input$filter_col %||% "") && length(input$filter_vals))
      d <- d[!chr(d[[input$filter_col]]) %in% input$filter_vals, , drop = FALSE]
    d
  })
  observeEvent(input$filter_col, {
    vals <- if (nzchar(input$filter_col)) sort(unique(chr(rv$d[[input$filter_col]]))) else character(0)
    updateSelectizeInput(session, "filter_vals", choices = vals, selected = character(0))
  }, ignoreInit = TRUE)

  # orders follow the chosen columns; the user's order is kept when the data change
  observeEvent(input$x, { req(nzchar(input$x)); updateTextInput(session, "x_order", value = paste(default_order(rv$d[[input$x]]), collapse = ", ")) })
  observeEvent(input$colour, { updateTextInput(session, "col_order", value = if (nzchar(input$colour)) paste(default_order(rv$d[[input$colour]]), collapse = ", ") else "") })
  observeEvent(list(rv$d, input$x, input$colour, input$box), {
    req(nzchar(input$x), nzchar(input$box %||% ""), input$box %in% names(rv$d))
    grp <- paste(rv$d[[input$x]], if (nzchar(input$colour %||% "")) rv$d[[input$colour]])
    reused <- any(tapply(grp, rv$d[[input$box]], function(v) length(unique(v))) > 1)
    updateCheckboxInput(session, "box_restart", value = reused && is.numeric(rv$d[[input$box]]))
  })
  x_levels <- reactive({ req(nzchar(input$x)); order_levels(dat()[[input$x]], input$x_order) })
  col_levels <- reactive(if (nzchar(input$colour %||% "")) order_levels(dat()[[input$colour]], input$col_order) else character(0))
  # the reference is asked, never preselected: it is the least susceptible group, known from the survival assay
  observeEvent(list(x_levels(), col_levels(), input$compare), {
    lv <- if (input$compare == "colour" && length(col_levels())) col_levels() else x_levels()
    cur <- isolate(input$ref)
    updateSelectInput(session, "ref", choices = c("Choose the reference group" = "", lv),
                      selected = if (!is.null(cur) && cur %in% lv) cur else "")
  })
  has_ref <- reactive(nzchar(input$ref %||% "") && input$ref %in% (if (settings()$compare == "colour") col_levels() else x_levels()))

  # colours chosen by the user, kept by group name
  custom_cols <- reactiveVal(list())
  observeEvent(input$pick_col, { cc <- custom_cols(); cc[[input$pick_col$lv]] <- input$pick_col$col; custom_cols(cc) })
  observeEvent(input$reset_cols, custom_cols(list()))
  output$colour_pickers <- renderUI({
    lv <- if (nzchar(input$colour %||% "")) col_levels() else x_levels()
    cur <- group_colours(lv, custom_cols())
    lapply(seq_along(lv), function(i) div(class = "d-flex align-items-center gap-2 mb-1",
      tags$input(type = "color", value = unname(cur[i]), style = "width:2.5rem; height:1.8rem; border:none; padding:0; background:none;",
                 onchange = sprintf("Shiny.setInputValue('pick_col', {lv: %s, col: this.value}, {priority: 'event'})",
                                    jsonlite::toJSON(lv[i], auto_unbox = TRUE))),
      span(lv[i])))
  })

  settings <- reactive(list(
    cols = custom_cols(),
    id = input$id %||% "", load_col = input$load_col, load_base = as.numeric(input$load_base), status_col = input$status_col %||% "",
    x = input$x, x_levels = x_levels(), colour = input$colour %||% "", col_levels = col_levels(), panel = input$panel %||% "",
    compare = if (nzchar(input$colour %||% "")) input$compare else "x", joint = input$joint == "joint", ref = input$ref %||% "",
    box = input$box %||% "", box_restart = isTRUE(input$box_restart), day = input$day %||% "",
    screen = TRUE, k = 2, alpha = 0.05, low_only = TRUE,      # fixed: the protocol's screen
    xlab = if (nzchar(input$xlab)) input$xlab else input$x, ylab = input$ylab,
    filter_col = input$filter_col %||% "", filter_vals = input$filter_vals %||% character(0)))

  # flies used: detected, with a load, and the analysis columns
  flies <- reactive({
    s <- settings(); d <- dat(); req(nzchar(s$load_col), nzchar(s$x))
    validate(need(nrow(d) > 0, "No fly is left: too many groups are left out."),
             need(s$panel == "" || !(s$panel %in% c(s$x, s$colour)), "Panels must be a different column from the x axis and the colour."),
             need(s$colour == "" || s$colour != s$x, "The colour must be a different column from the x axis."),
             need(!(s$box %in% c(s$x, s$colour, s$panel)) || s$box == "", sprintf("%s is chosen as the box and as a factor to compare. The box is the unit flies share; compare what the flies are (sex, genotype, age), with the box as random effect.", s$box)),
             need(!(s$day %in% c(s$x, s$colour, s$panel)) || s$day == "", sprintf("%s is chosen as the day of injection and as a factor to compare.", s$day)),
             need(s$x != s$id || s$id == "", "The fly identifier cannot be the groups to compare."),
             need(length(unique(d[[s$x]])) >= 2, sprintf("%s has a single value: choose a column with at least two groups to compare.", s$x)),
             need(length(unique(d[[s$x]])) <= max(2, nrow(d) / 3), sprintf("%s has almost one value per fly: choose the column with the groups to compare.", s$x)))
    y <- num(d[[s$load_col]]) * log10(s$load_base)                      # to log10
    status <- if (nzchar(s$status_col)) as.character(d[[s$status_col]]) else ifelse(is.finite(y), "Detected", "No load")
    keep <- status == "Detected" & is.finite(y)
    pf <- data.frame(.row = seq_len(nrow(d)), .id = if (nzchar(s$id)) as.character(d[[s$id]]) else paste0("row ", seq_len(nrow(d))),
                     .y = y, Status = status, stringsAsFactors = FALSE)
    pf$.x <- factor(chr(d[[s$x]]), levels = s$x_levels)
    pf$.col <- if (nzchar(s$colour)) factor(chr(d[[s$colour]]), levels = s$col_levels) else factor("")
    pf$.panel <- if (nzchar(s$panel)) factor(chr(d[[s$panel]])) else factor("All")
    cmp <- if (s$compare == "colour") pf$.col else pf$.x
    pf$.cmp <- relevel(droplevels(cmp), ref = if (s$ref %in% levels(droplevels(cmp))) s$ref else levels(droplevels(cmp))[1])
    if (nzchar(s$box)) {
      parts <- list(pf$.panel, pf$.x, if (nzchar(s$colour)) pf$.col, d[[s$box]], if (nzchar(s$day)) d[[s$day]])
      pf$.box <- if (s$box_restart) do.call(interaction, c(Filter(Negate(is.null), parts), drop = TRUE, sep = "_")) else factor(d[[s$box]])
    }
    if (nzchar(s$day)) pf$.day <- factor(d[[s$day]])
    pf$.cell <- do.call(paste, c(Filter(Negate(is.null), list(if (nzchar(s$panel)) pf$.panel, pf$.x, if (nzchar(s$colour)) pf$.col)), sep = " / "))
    list(all = pf, used = pf[keep, ])
  })

  screened <- reactive({
    pf <- flies()$used; s <- settings()
    if (s$screen) screen_groups(pf, s$k, s$alpha, s$low_only)
    else { pf$Flagged <- FALSE; pf$Removed <- FALSE; pf$Direction <- ""; list(pf = pf, summary = NULL) }
  })
  results <- reactive({ req(has_ref()); s <- settings(); pf <- screened()$pf; analyse(if (s$screen) pf[!pf$Removed, ] else pf, s) })
  results_all <- reactive({ req(has_ref()); s <- settings(); if (s$screen) analyse(screened()$pf, s) else NULL })

  output$data_info <- renderText({
    if (is.null(rv$d)) return("No data loaded.")
    n <- nrow(dat())
    sprintf("%s: %d rows%s, %d columns.", rv$name, nrow(rv$d), if (n < nrow(rv$d)) sprintf(" (%d after leaving out groups)", n) else "", ncol(rv$d))
  })

  output$messages <- renderUI({
    if (is.null(rv$d)) return(div(class = "alert alert-info", "Load the per-fly table saved by CFU Drop Counter, or try the example."))
    fl <- tryCatch(flies(), error = function(e) e)
    if (inherits(fl, "error")) return(div(class = "alert alert-danger", conditionMessage(fl)))
    s <- settings(); w <- character(0)
    st <- table(fl$all$Status[!fl$all$.row %in% fl$used$.row])
    if (length(st)) w <- c(w, paste0("Flies without a load, left out: ", paste(sprintf("%s %d", names(st), st), collapse = ", "), "."))
    if (s$screen) {
      sm <- screened()$summary
      small <- sm$Group[sm$Flies < 25]
      if (length(small)) w <- c(w, sprintf("Fewer than 25 flies in %s: the screen is less reliable there. Pool experimental replicates of the same group where possible.", paste(small, collapse = ", ")))
      if (sum(sm$Flagged_high) > 0) w <- c(w, sprintf("%d high value(s) were flagged; they are %s.", sum(sm$Flagged_high), if (s$low_only) "shown in the Screen tab but kept" else "removed"))
    }
    if (nzchar(s$box)) {
      pf <- fl$used
      if (all(tapply(as.character(pf$.cmp), pf$.box, function(v) length(unique(v))) == 1)) {
        nb <- tapply(as.character(pf$.box), pf$.cell, function(v) length(unique(v)))
        w <- c(w, sprintf("Each box holds a single one of the groups compared, so the groups are compared between boxes: the comparison rests on %s box%s per group, more than on the number of flies.",
                          if (min(nb) == max(nb)) min(nb) else paste(min(nb), "to", max(nb)), if (max(nb) > 1) "es" else ""))
      }
    }
    if (!has_ref()) return(tagList(
      div(class = "alert alert-info", strong("Choose the reference group."),
          " Take the least susceptible group, the one with the lowest mortality in the survival assay. Differences in PLUD are then given relative to it, and read in the same direction as the hazard ratios (see the protocol)."),
      if (length(w)) div(class = "alert alert-warning", tags$ul(lapply(w, tags$li)))))
    for (nm in names(results())) {
      f <- results()[[nm]]
      if (inherits(f, "error")) w <- c(w, sprintf("%s: the model could not be fitted (%s).", nm, conditionMessage(f)))
      else if (isTRUE(f$singular)) w <- c(w, sprintf("%s: a random-effect variance was estimated at zero; the model is then equivalent to one without it.", if (nm == "All") "Model" else nm))
    }
    if (length(w)) div(class = "alert alert-warning", tags$ul(lapply(w, tags$li)))
  })

  the_plot <- reactive({
    s <- settings(); p1 <- plot_plud(screened()$pf, s); p2 <- if (has_ref()) plot_diffs(results(), s)
    if (is.null(p2)) p1 else p1 + p2 + plot_layout(widths = c(2, 1))
  })
  output$plot <- renderPlot(the_plot(), res = 96)

  table_of <- function(f) {
    m <- f$means; d <- f$diffs
    key_m <- paste(m$Within, m$Group); key_d <- paste(d$Within, d$Group)
    dd <- d[match(key_m, key_d), ]
    data.frame(Within = m$Within, Group = m$Group, Flies = m$Flies, `Mean log10 PLUD` = sprintf("%.2f", m$Mean),
               `95% CI` = sprintf("%.2f to %.2f", m$Lower, m$Upper),
               `Difference to reference` = ifelse(m$Group == f$ref, "reference", sprintf("%+.2f", dd$Difference)),
               `95% CI of difference` = ifelse(m$Group == f$ref, "", sprintf("%.2f to %.2f", dd$Lower, dd$Upper)),
               `Difference in bacteria per fly` = ifelse(m$Group == f$ref, "", fmt_bact(dd$Bacteria_per_fly)),
               `95% CI in bacteria per fly` = ifelse(m$Group == f$ref, "", paste(fmt_bact(dd$Bacteria_lower), "to", fmt_bact(dd$Bacteria_upper))),
               check.names = FALSE)
  }
  output$results_ui <- renderUI({
    s <- settings()
    sets <- c(list(main = results()), if (s$screen) list(all = results_all()))
    tagList(lapply(names(sets), function(set) {
      res <- sets[[set]]
      tagList(h5(if (set == "main") (if (s$screen) "After the screen" else "All flies") else "All flies, without the screen"),
        lapply(names(res), function(nm) {
          f <- res[[nm]]; if (inherits(f, "error")) return(p(sprintf("%s: the model could not be fitted.", nm)))
          tab <- table_of(f); if (all(is.na(tab$Within))) tab$Within <- NULL
          tagList(if (nm != "All") h6(nm), renderTable(tab),
            p(class = "small", sprintf("Likelihood ratio test: chi2 = %.2f, df = %d, p = %s.", f$lrt$Chisq, f$lrt$df, fmt_p(f$lrt$p)), br(),
              sprintf("Within-group SD of log10 PLUD = %.2f (residual SD = %.2f%s).", f$total_sd, f$resid_sd,
                      if (length(f$re_sd)) paste0("; ", paste(sprintf("SD of %s = %.2f", names(f$re_sd), f$re_sd), collapse = ", ")) else ""), br(),
              sprintf("With %.0f flies per group, a difference of %.2f log10 is detected with 80%% power at 5%% (2.8 x SD x sqrt(2/n), ignoring the boxes). Confidence intervals: %s.",
                      f$n_h, f$detectable, f$ci_kind)))
        }))
    }))
  })
  results_csv <- reactive({
    s <- settings(); sets <- c(list(screened = results()), if (s$screen) list(all_flies = results_all()))
    do.call(rbind, lapply(names(sets), function(set) do.call(rbind, lapply(names(sets[[set]]), function(nm) {
      f <- sets[[set]][[nm]]; if (inherits(f, "error")) return(NULL)
      cbind(Flies_used = if (set == "screened" && s$screen) "after screen" else "all", Analysis = nm, table_of(f),
            LRT_chisq = round(f$lrt$Chisq, 3), LRT_df = f$lrt$df, LRT_p = signif(f$lrt$p, 3), Within_group_SD = round(f$total_sd, 3))
    }))))
  })

  output$screen_table <- renderTable({ validate(need(settings()$screen, "The screen is switched off.")); screened()$summary })
  output$flagged_table <- renderTable({
    sc <- screened(); f <- sc$pf[sc$pf$Flagged, ]
    validate(need(nrow(f) > 0, "No fly was flagged."))
    data.frame(Fly = f$.id, Group = f$.cell, log10_load = round(f$.y, 2), Side = f$Direction, Removed = ifelse(f$Removed, "yes", "no"))
  })
  pf_out <- reactive({
    fl <- flies()$all; sc <- screened()$pf; s <- settings()
    out <- data.frame(ID = fl$.id, check.names = FALSE)
    if (nzchar(s$panel)) out$panel <- as.character(fl$.panel)
    out$x <- as.character(fl$.x)
    if (nzchar(s$colour)) out$colour <- as.character(fl$.col)
    if (nzchar(s$box)) out$box <- as.character(fl$.box)
    if (nzchar(s$day)) out$day <- as.character(fl$.day)
    out$log10_load <- round(fl$.y, 4); out$Status <- fl$Status
    out$Screen <- ifelse(fl$.row %in% sc$.row[sc$Removed], "removed", ifelse(fl$.row %in% sc$.row[sc$Flagged], "flagged, kept", ""))
    out[out$Status == "Detected" & is.finite(out$log10_load), ]
  })
  output$pf_table <- renderTable(head(pf_out(), 50))

  methods <- reactive(methods_text(results(), settings(), screened(), nrow(screened()$pf)))
  output$methods <- renderText(methods())
  code <- reactive(r_code(settings()))
  output$rcode <- renderText(code())

  observeEvent(input$dl_pdf, send_file(session, "PLUD_plots.pdf", "application/pdf", function(f) ggsave(f, the_plot(), width = 11, height = 5.5)))
  observeEvent(input$dl_png, send_file(session, "PLUD_plots.png", "image/png", function(f) ggsave(f, the_plot(), width = 11, height = 5.5, dpi = 300)))
  observeEvent(input$dl_svg, send_file(session, "PLUD_plots.svg", "image/svg+xml", function(f) ggsave(f, the_plot(), width = 11, height = 5.5, device = svglite::svglite)))
  observeEvent(input$dl_results, send_file(session, "PLUD_results.csv", "text/csv", function(f) write.table(results_csv(), f, sep = ";", row.names = FALSE)))
  observeEvent(input$dl_pf, send_file(session, "PLUD_flies.csv", "text/csv", function(f) write.csv(pf_out(), f, row.names = FALSE)))
  observeEvent(input$dl_methods, send_file(session, "PLUD_methods.txt", "text/plain", function(f) writeLines(methods(), f)))
  observeEvent(input$dl_code, send_file(session, "PLUD_analysis.R", "text/plain", function(f) writeLines(code(), f)))
}

shinyApp(ui, server)
