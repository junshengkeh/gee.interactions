#' Plot 2/3/4-way interactions for GEE
#'
#' @description
#' `gee_plot` wraps around (and combines) `geeglm` and `ggplot` to provide
#'     visualizations for 2-, 3-, and 4-way interactions. For 3- and 4-way
#'     interactions, the output is a faceted plot.
#'
#' @param data
#'     A `data.frame` object
#' @param dv
#'     Dependent variable
#' @param dv_label
#'     Label name for `dv`
#' @param a
#'     Independent variable (variable plotted on the x-axis)
#' @param a_label
#'     Label name for `a`
#' @param b
#'     Independent variable (variable plotted as simple slopes)
#' @param b_label
#'     Label name for `b`
#' @param c
#'     Independent variable (variable plotted along columns)
#' @param c_label
#'     Label name for `c`
#' @param d
#'     Independent variable (variable plotted along rows)
#' @param d_label
#'     Label name for `d`
#' @param ...
#'     Covariates
#' @param corstr
#'     (`geeglm` argument) See documentation for `geeglm`. Defaults to
#'     "independence".
#' @param int
#'     A character string specifying the interaction to be probed. Accepts the
#'     following: "2-way", "3-way", and "4-way".
#' @param line_size
#'     (`ggplot` argument) Line widths for plot. Defaults to 0.8.
#' @param line_color_low
#'     (`ggplot` argument) Line color when level of `b` is low.
#'     Defaults to bright red (#EE4B2B").
#' @param line_color_high
#'     (`ggplot` argument) Line color when level of `b` is high.
#'     Defaults to dark blue ("#00AFBB").
#' @param x_tick_break
#'     (`ggplot` argument) Tick mark breaks for x-axis
#' @param x_tick_label
#'     (`ggplot` argument) Tick mark labels for x-axis. If `x_tick_label` is
#'     specified, the length must be the same as `x_tick_break`.
#' @param x_limit
#'     (`ggplot` argument) Upper and lower limits of x-axis
#'
#' @details
#' Since geeglm only works for complete data, NAs are removed via the following
#'     specification in `geeglm` (`data = na.omit(data)`).
#'
#' @return
#' When `int = "2-way"`, `gee_plot` returns a non-faceted plot; when
#'     `int = "3-way"`, `gee_plot` returns a faceted plot consisting of 3
#'     sub-plots; and when `int = "4-way"`, `gee_plot` returns a faceted plot
#'     consisting of 9 sub-plots.
#'
#' @author Jun Sheng Keh, \email{junshengkeh@gmail.com}
#'
#' @seealso \code{\link{geeglm}}
#'
#' @examples
#' gee_plot(data = df1c,
#'          dv = "MEA_c", dv_label = "Meaning",
#'          a = "POS", a_label = "Positive Affect",
#'          b = "NEG", b_label = "Negative Affect",
#'          "EXT", "NEU", "tMEA", "OPT", "SubHea",
#'          corstr = "exchangeable",
#'          int = "2-way",
#'          x_tick_break = c(2.9, 3.8, 4.7),
#'          x_tick_label = c("low", "mid", "high"),
#'          x_limit = c(1, 5))
#'
#' @export
gee_plot <- function(data, dv, dv_label, a, a_label, b, b_label,
                     c = NULL, c_label = NULL, d = NULL, d_label = NULL,
                     ...,
                     corstr = "independence", int = c("2-way", "3-way", "4-way"),
                     line_color_low = "#EE4B2B", line_color_high = "#00AFBB",
                     line_size = 0.8, x_tick_break = NULL,
                     x_tick_label = NULL, x_limit = NULL) {

  # Setup formula 1
  if (int == "2-way") {
    formula.1 <- as.formula(
      paste0(dv, " ~ ",
             paste(
               paste(paste0("(", a),
                     paste0(b, ")"),
                     sep = "*"),
               ...,
               sep = "+")))
  } else if (int == "3-way") {
    formula.1 <- as.formula(
      paste0(dv, " ~ ",
             paste(
               paste(paste0("(", a),
                     b,
                     paste0(c, ")"),
                     sep = "*"),
               ...,
               sep = "+")))
  } else if (int == "4-way") {
    formula.1 <- as.formula(
      paste0(dv, " ~ ",
             paste(
               paste(paste0("(", a),
                     b,
                     c,
                     paste0(d, ")"),
                     sep = "*"),
               ...,
               sep = "+")))
  }

  # Run GEE model 1
  mod <- geeglm(formula.1,
                data = na.omit(data),
                id = SubID,
                corstr = corstr)

  # FUNCTION: Find simple slope
  simple_intercept <- function(b_level = c("low", "mean", "high"),
                               c_level = c("low", "mean", "high"),
                               d_level = c("low", "mean", "high")) {

    # FUNCTION: Extract coefficients from GEE models
    #     mod       GEE model
    #     k         coefficient (k)
    #     l         coefficient (k*l)
    #     m         coefficient (k*l*m)
    extract_coef <- function(mod, k, l = NULL, m = NULL) {
      if (int == "2-way") {
        mod$coefficients[k]
      } else if (int == "3-way") {
        if (is.null(l)) {
          mod$coefficients[k]
        } else {
          if (is_na(mod$coefficients[paste0(l, ":", k)])) {
            mod$coefficients[paste0(k, ":", l)]          # k*l
          } else {
            mod$coefficients[paste0(l, ":", k)]          # l*k
          }
        }
      } else if (int == "4-way") {
        if (is.null(c(l, m))) {
          mod$coefficients[k]
        } else if (is.null(m)) {
          if (is_na(mod$coefficients[paste0(l, ":", k)])) {
            mod$coefficients[paste0(k, ":", l)]          # k*l
          } else if (is_na(mod$coefficients[paste0(k, ":", l)])) {
            mod$coefficients[paste0(l, ":", k)]          # l*k
          }
        } else {
          if (is_na(mod$coefficients[paste0(m, ":", k, ":", l)])) {
            mod$coefficients[paste0(k, ":", l, ":", m)]  # k*l*m
          } else if (is_na(mod$coefficients[paste0(k, ":", l, ":", m)])) {
            mod$coefficients[paste0(k, ":", m, ":", l)]  # k*m*l
          } else if (is_na(mod$coefficients[paste0(k, ":", m, ":", l)])) {
            mod$coefficients[paste0(l, ":", k, ":", m)]  # l*k*m
          } else if (is_na(mod$coefficients[paste0(l, ":", k, ":", m)])) {
            mod$coefficients[paste0(l, ":", m, ":", k)]  # l*m*k
          } else if (is_na(mod$coefficients[paste0(l, ":", m, ":", k)])) {
            mod$coefficients[paste0(m, ":", l, ":", k)]  # m*l*k
          } else if (is_na(mod$coefficients[paste0(l, ":", m, ":", k)])) {
            mod$coefficients[paste0(m, ":", k, ":", l)]  # m*k*l
          }
        }
      }
    }

    # FUNCTION: Compute multiplier (for slope levels)
    #     k
    #     level     when k is at low, mean, or high level
    slope_level <- function(k, level = c("low", "mean", "high")) {
      if (level == "low") {
        mean(data[, k], na.rm = TRUE) - sd(data[, k], na.rm = TRUE)
      } else if (level == "mean") {
        mean(data[, k], na.rm = TRUE)
      } else if (level == "high") {
        mean(data[, k], na.rm = TRUE) + sd(data[, k], na.rm = TRUE)
      }
    }

    # Compute simple intercept
    #     The simple intercept can be found by finding the value of the dv when
    #     the predictor is 0. Thus, to find the simple intercept, find dv after
    #     substituting (a = 0) into the general formula.
    #
    #     For a 2-way interaction:
    #       (1) dv = b0 + b1(a) + b2(b) + b3(ab)
    #       (2) (When a = 0):
    #           dv = b0 + b2(b)
    #
    #     For a 3-way interaction:
    #       (1) dv = b0 + b1(a) + b2(b) + b3(c) + b4(ab) + b5(ac) + b6(bc) +
    #                b7(abc)
    #       (2) (When a = 0):
    #           dv = b0 + b2(b) + b3(c) + b6(b*c)
    #
    #     For a 4-way interaction:
    #       (1) dv = b0 + b1(a) + b2(b) + b3(c) + b4(d) + b5(ab) + b6(ac) +
    #                b7(ad) + b8(bc) + b9(bd) + b10(cd) + b11(abc) + b12(abd) +
    #                b13(acd) + b14(bcd) + b15(abcd)
    #       (2) (When a = 0):
    #           dv = b0 + b2(b) + b3(c) + b4(d) + b8(bc) + b9(bd) + b10(cd) +
    #                b14(bcd)
    #
    if (int == "2-way") {
      extract_coef(mod, "(Intercept)") +
        (extract_coef(mod, b) * slope_level(b, b_level))
    } else if (int == "3-way") {
      extract_coef(mod, "(Intercept)") +
        (extract_coef(mod, b) * slope_level(b, b_level)) +
        (extract_coef(mod, c) * slope_level(c, c_level)) +
        (extract_coef(mod, b, c) *
           slope_level(b, b_level) * slope_level(c, c_level))
    } else if (int == "4-way") {
      extract_coef(mod, "(Intercept)") +
        (extract_coef(mod, b) * slope_level(b, b_level)) +
        (extract_coef(mod, c) * slope_level(c, c_level)) +
        (extract_coef(mod, d) * slope_level(d, d_level)) +
        (extract_coef(mod, b, c) *
           slope_level(b, b_level) * slope_level(c, c_level)) +
        (extract_coef(mod, b, d) *
           slope_level(b, b_level) * slope_level(d, d_level)) +
        (extract_coef(mod, c, d) *
           slope_level(c, c_level) * slope_level(d, d_level)) +
        (extract_coef(mod, b, c, d) *
           slope_level(b, b_level) * slope_level(c, c_level) *
           slope_level(d, d_level))

    }
  }

  # FUNCTION: Find simple slope
  #     The simple slope is the value of the first predictor (a) when GEE is
  #     run with the corresponding (e.g. low or high) values of the other
  #     predictors (i.e. b, c, and/or d).
  simple_slope <- function(b_level = c("low", "mean", "high"),
                           c_level = c("low", "mean", "high"),
                           d_level = c("low", "mean", "high")) {

    a_ensym <- ensym(a)
    a_low  <- as.name(paste0(ensym(a),"_low"))
    a_mean <- as.name(paste0(ensym(a),"_mean"))
    a_high <- as.name(paste0(ensym(a),"_high"))

    b_ensym <- ensym(b)
    b_low  <- as.name(paste0(ensym(b),"_low"))
    b_mean <- as.name(paste0(ensym(b),"_mean"))
    b_high <- as.name(paste0(ensym(b),"_high"))

    if (int == "2-way") {
      data2 <- data %>%
        mutate({{a_low}}  :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE) -
                                  sd({{a_ensym}}, na.rm = TRUE)),
               {{a_mean}} :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE)),
               {{a_high}} :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE) +
                                  sd({{a_ensym}}, na.rm = TRUE)),
               {{b_low}}  :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE) -
                                  sd({{b_ensym}}, na.rm = TRUE)),
               {{b_mean}} :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE)),
               {{b_high}} :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE) +
                                  sd({{b_ensym}}, na.rm = TRUE)))

      formula.2 <- as.formula(
        paste0(dv, " ~ ",
               paste(
                 paste(paste0("(", a),
                       paste0(b, "_", b_level, ")"),
                       sep = "*"),
                 ...,
                 sep = "+")))

    } else if (int == "3-way") {
      c_ensym <- ensym(c)
      c_low  <- as.name(paste0(ensym(c),"_low"))
      c_mean <- as.name(paste0(ensym(c),"_mean"))
      c_high <- as.name(paste0(ensym(c),"_high"))

      data2 <- data %>%
        mutate({{a_low}}  :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE) -
                                  sd({{a_ensym}}, na.rm = TRUE)),
               {{a_mean}} :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE)),
               {{a_high}} :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE) +
                                  sd({{a_ensym}}, na.rm = TRUE)),
               {{b_low}}  :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE) -
                                  sd({{b_ensym}}, na.rm = TRUE)),
               {{b_mean}} :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE)),
               {{b_high}} :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE) +
                                  sd({{b_ensym}}, na.rm = TRUE)),
               {{c_low}}  :=
                 {{c_ensym}} - (mean({{c_ensym}}, na.rm = TRUE) -
                                  sd({{c_ensym}}, na.rm = TRUE)),
               {{c_mean}} :=
                 {{c_ensym}} - (mean({{c_ensym}}, na.rm = TRUE)),
               {{c_high}} :=
                 {{c_ensym}} - (mean({{c_ensym}}, na.rm = TRUE) +
                                  sd({{c_ensym}}, na.rm = TRUE)))

      formula.2 <- as.formula(
        paste0(dv, " ~ ",
               paste(
                 paste(paste0("(", a),
                       paste0(b, "_", b_level),
                       paste0(c, "_", c_level, ")"),
                       sep = "*"),
                 ...,
                 sep = "+")))

    } else if (int == "4-way") {
      c_ensym <- ensym(c)
      c_low  <- as.name(paste0(ensym(c),"_low"))
      c_mean <- as.name(paste0(ensym(c),"_mean"))
      c_high <- as.name(paste0(ensym(c),"_high"))

      d_ensym <- ensym(d)
      d_low  <- as.name(paste0(ensym(d),"_low"))
      d_mean <- as.name(paste0(ensym(d),"_mean"))
      d_high <- as.name(paste0(ensym(d),"_high"))

      data2 <- data %>%
        mutate({{a_low}}  :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE) -
                                  sd({{a_ensym}}, na.rm = TRUE)),
               {{a_mean}} :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE)),
               {{a_high}} :=
                 {{a_ensym}} - (mean({{a_ensym}}, na.rm = TRUE) +
                                  sd({{a_ensym}}, na.rm = TRUE)),
               {{b_low}}  :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE) -
                                  sd({{b_ensym}}, na.rm = TRUE)),
               {{b_mean}} :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE)),
               {{b_high}} :=
                 {{b_ensym}} - (mean({{b_ensym}}, na.rm = TRUE) +
                                  sd({{b_ensym}}, na.rm = TRUE)),
               {{c_low}}  :=
                 {{c_ensym}} - (mean({{c_ensym}}, na.rm = TRUE) -
                                  sd({{c_ensym}}, na.rm = TRUE)),
               {{c_mean}} :=
                 {{c_ensym}} - (mean({{c_ensym}}, na.rm = TRUE)),
               {{c_high}} :=
                 {{c_ensym}} - (mean({{c_ensym}}, na.rm = TRUE) +
                                  sd({{c_ensym}}, na.rm = TRUE)),
               {{d_low}}  :=
                 {{d_ensym}} - (mean({{d_ensym}}, na.rm = TRUE) -
                                  sd({{d_ensym}}, na.rm = TRUE)),
               {{d_mean}} :=
                 {{d_ensym}} - (mean({{d_ensym}}, na.rm = TRUE)),
               {{d_high}} :=
                 {{d_ensym}} - (mean({{d_ensym}}, na.rm = TRUE) +
                                  sd({{d_ensym}}, na.rm = TRUE)))

      formula.2 <- as.formula(
        paste0(dv, " ~ ",
               paste(
                 paste(paste0("(", a),
                       paste0(b, "_", b_level),
                       paste0(c, "_", c_level),
                       paste0(d, "_", d_level, ")"),
                       sep = "*"),
                 ...,
                 sep = "+")))

    }

    # Run GEE model 2
    mod <- geeglm(formula.2,
                  data = na.omit(data2),
                  id = SubID,
                  corstr = corstr)

    # extract value of simple slope
    mod$coefficients[a]
  }

  # FUNCTION: Creates plot showing predicted values of one predictor (a) at
  #           specific levels of other predictors (b, c, and/or d)
  #     int_low         value of simple intercept
  #     int_high
  #     slope_low       value of simple slope
  #     slope_high
  plot_expand <- function(int_low, slope_low,
                          int_high, slope_high) {

    # Convert to symbols
    x_var <- ensym(a)
    y_var <- ensym(dv)

    # plot settings
    ggplot(data, aes(x = {{x_var}}, y = {{y_var}})) +
      geom_blank() +
      geom_abline(aes(intercept = int_low,
                      slope = slope_low,
                      color = paste("Low", b)),
                  size = line_size) +
      geom_abline(aes(slope = slope_high,
                      intercept = int_high,
                      color = paste("High", b)),
                  size = line_size) +
      # general aesthetics
      theme_classic() +
      # axis aesthetics
      scale_x_continuous(breaks = x_tick_break,
                         labels = x_tick_label,
                         limits = x_limit) +
      theme(axis.title = element_blank()) +
      # legend aesthetics
      labs(colour = b_label) +
      scale_colour_manual(values = c(line_color_low, line_color_high)) +
      theme(legend.position = "none")

  }

  if (int == "2-way") {
    # left plot
    p1 <- plot_expand(int_low = simple_intercept("low"),
                      slope_low = simple_slope("low"),
                      int_high = simple_intercept("high"),
                      slope_high = simple_slope("high"))

    # Plot settings
    #   - Axis titles/subtitles size
    size_title    <- 11
    size_subtitle <- 8.5
    #   - Axis titles/subtitles
    left    <- text_grob(dv_label, rot = 90, size = size_title)
    bottom  <- text_grob(a_label, size = size_title)
    #   - Create common legend
    legend <- as_ggplot(
      get_legend(
        plot_expand(int_low = simple_intercept("low", "low"),
                    slope_low = simple_slope("low", "low"),
                    int_high = simple_intercept("high", "low"),
                    slope_high = simple_slope("high", "low")),
        position = "bottom"))
    #   - List of grobs
    grobs <- list(p1,
                  left, bottom,
                  legend)
    #   - Create layout
    lay <- rbind(c(2 , 1 ),
                 c(NA, 3 ),
                 c(NA, 4 ))

    # Create plot
    grid.arrange(
      arrangeGrob(grobs = grobs,
                  layout_matrix = lay,
                  widths = c(4, 96),
                  heights = c(88, 4, 8)))
  } else if (int == "3-way") {
    # left plot
    p1 <- plot_expand(int_low = simple_intercept("low", "low"),
                      slope_low = simple_slope("low", "low"),
                      int_high = simple_intercept("high", "low"),
                      slope_high = simple_slope("high", "low"))
    # middle plot
    p2 <- plot_expand(int_low = simple_intercept("low", "mean"),
                      slope_low = simple_slope("low", "mean"),
                      int_high = simple_intercept("high", "mean"),
                      slope_high = simple_slope("high", "mean")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # right plot
    p3 <- plot_expand(int_low = simple_intercept("low", "high"),
                      slope_low = simple_slope("low", "high"),
                      int_high = simple_intercept("high", "high"),
                      slope_high = simple_slope("high", "high")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

    # Plot settings
    #   - Axis titles/subtitles size
    size_title    <- 11
    size_subtitle <- 8.5
    #   - Axis titles/subtitles
    left    <- text_grob(dv_label, rot = 90, size = size_title)
    bottom  <- text_grob(a_label, size = size_title)
    top     <- text_grob(c_label, size = size_title)
    top_1   <- text_grob("(-1 SD)", size = size_subtitle)
    top_2   <- text_grob("(Mean)", size = size_subtitle)
    top_3   <- text_grob("(+1 SD)", size = size_subtitle)
    #   - Create common legend
    legend <- as_ggplot(
      get_legend(
        plot_expand(int_low = simple_intercept("low", "low"),
                    slope_low = simple_slope("low", "low"),
                    int_high = simple_intercept("high", "low"),
                    slope_high = simple_slope("high", "low")),
        position = "bottom"))
    #   - List of grobs
    grobs <- list(p1, p2, p3,
                  left, bottom, top, top_1, top_2, top_3,
                  legend)
    #   - Create layout
    lay <- rbind(c(NA, NA, 6 , NA),
                 c(NA, 7 , 8 , 9 ),
                 c(4 , 1 , 2 , 3 ),
                 c(NA, NA, 5 , NA),
                 c(NA, 10, 10, 10))

    # Create plot
    grid.arrange(
      arrangeGrob(grobs = grobs,
                  layout_matrix = lay,
                  widths = c(4, 32, 32, 32),
                  heights = c(6, 2, 80, 4, 8)))

  } else if (int == "4-way") {
    # top-left plot
    p1 <- plot_expand(int_low = simple_intercept("low", "low", "high"),
                      slope_low = simple_slope("low", "low", "high"),
                      int_high = simple_intercept("high", "low", "high"),
                      slope_high = simple_slope("high", "low", "high"))
    # top-middle plot
    p2 <- plot_expand(int_low = simple_intercept("low", "mean", "high"),
                      slope_low = simple_slope("low", "mean", "high"),
                      int_high = simple_intercept("high", "mean", "high"),
                      slope_high = simple_slope("high", "mean", "high")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # top-right plot
    p3 <- plot_expand(int_low = simple_intercept("low", "high", "high"),
                      slope_low = simple_slope("low", "high", "high"),
                      int_high = simple_intercept("high", "high", "high"),
                      slope_high = simple_slope("high", "high", "high")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # middle-left plot
    p4 <- plot_expand(int_low = simple_intercept("low", "low", "mean"),
                      slope_low = simple_slope("low", "low", "mean"),
                      int_high = simple_intercept("high", "low", "mean"),
                      slope_high = simple_slope("high", "low", "mean"))
    # middle-middle plot
    p5 <- plot_expand(int_low = simple_intercept("low", "mean", "mean"),
                      slope_low = simple_slope("low", "mean", "mean"),
                      int_high = simple_intercept("high", "mean", "mean"),
                      slope_high = simple_slope("high", "mean", "mean")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # middle-right plot
    p6 <- plot_expand(int_low = simple_intercept("low", "high", "mean"),
                      slope_low = simple_slope("low", "high", "mean"),
                      int_high = simple_intercept("high", "high", "mean"),
                      slope_high = simple_slope("high", "high", "mean")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # bottom-left plot
    p7 <- plot_expand(int_low = simple_intercept("low", "low", "low"),
                      slope_low = simple_slope("low", "low", "low"),
                      int_high = simple_intercept("high", "low", "low"),
                      slope_high = simple_slope("high", "low", "low"))
    # bottom-middle plot
    p8 <- plot_expand(int_low = simple_intercept("low", "mean", "low"),
                      slope_low = simple_slope("low", "mean", "low"),
                      int_high = simple_intercept("high", "mean", "low"),
                      slope_high = simple_slope("high", "mean", "low")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # bottom-right plot
    p9 <- plot_expand(int_low = simple_intercept("low", "high", "low"),
                      slope_low = simple_slope("low", "high", "low"),
                      int_high = simple_intercept("high", "high", "low"),
                      slope_high = simple_slope("high", "high", "low")) +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())

    # Plot settings
    #   - Axis titles/subtitles size
    size_title    <- 11
    size_subtitle <- 8.5
    #   - Axis titles/subtitles
    left    <- text_grob(dv_label, rot = 90, size = size_title)
    bottom  <- text_grob(a_label, size = size_title)
    top     <- text_grob(c_label, size = size_title)
    top_1   <- text_grob("(-1 SD)", size = size_subtitle)
    top_2   <- text_grob("(Mean)", size = size_subtitle)
    top_3   <- text_grob("(+1 SD)", size = size_subtitle)
    right   <- text_grob(d_label, rot = 270, size = size_title)
    right_1 <- text_grob("(+1 SD)", rot = 270, size = size_subtitle)
    right_2 <- text_grob("(Mean)", rot = 270, size = size_subtitle)
    right_3 <- text_grob("(-1 SD)", rot = 270, size = size_subtitle)
    #   - Create common legend
    legend <- as_ggplot(
      get_legend(
        plot_expand(int_low = simple_intercept("low", "low", "high"),
                    slope_low = simple_slope("low", "low", "high"),
                    int_high = simple_intercept("high", "low", "high"),
                    slope_high = simple_slope("high", "low", "high")),
        position = "bottom"))
    #   - List of grobs
    grobs <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9,
                  left, bottom, top, top_1, top_2, top_3, right,
                  right_1, right_2, right_3,
                  legend)
    #   - Create layout
    lay <- rbind(c(NA, NA, 12, NA, NA, NA),
                 c(NA, 13, 14, 15, NA, NA),
                 c(NA, 1 , 2 , 3 , 17, NA),
                 c(10, 4 , 5 , 6 , 18, 16),
                 c(NA, 7 , 8 , 9 , 19, NA),
                 c(NA, NA, 11, NA, NA, NA),
                 c(NA, 20, 20, 20, NA, NA))

    # Create plot
    grid.arrange(
      arrangeGrob(grobs = grobs,
                  layout_matrix = lay,
                  widths = c(4, (86/3), (86/3), (86/3), 4, 6),
                  heights = c(6, 2, (80/3), (80/3), (80/3), 4, 8)))

  }
}
