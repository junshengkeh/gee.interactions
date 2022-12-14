% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gee_plot.R
\name{gee_plot}
\alias{gee_plot}
\title{Plot 2/3/4-way interactions for GEE}
\usage{
gee_plot(
  data,
  dv,
  dv_label,
  a,
  a_label,
  b,
  b_label,
  c = NULL,
  c_label = NULL,
  d = NULL,
  d_label = NULL,
  ...,
  corstr = "independence",
  int = c("2-way", "3-way", "4-way"),
  line_color_low = "#EE4B2B",
  line_color_high = "#00AFBB",
  line_size = 0.8,
  x_tick_break = NULL,
  x_tick_label = NULL,
  x_limit = NULL
)
}
\arguments{
\item{data}{A \code{data.frame} object}

\item{dv}{Dependent variable}

\item{dv_label}{Label name for \code{dv}}

\item{a}{Independent variable (variable plotted on the x-axis)}

\item{a_label}{Label name for \code{a}}

\item{b}{Independent variable (variable plotted as simple slopes)}

\item{b_label}{Label name for \code{b}}

\item{c}{Independent variable (variable plotted along columns)}

\item{c_label}{Label name for \code{c}}

\item{d}{Independent variable (variable plotted along rows)}

\item{d_label}{Label name for \code{d}}

\item{...}{Covariates}

\item{corstr}{(\code{geeglm} argument) See documentation for \code{geeglm}. Defaults to
"independence".}

\item{int}{A character string specifying the interaction to be probed. Accepts the
following: "2-way", "3-way", and "4-way".}

\item{line_color_low}{(\code{ggplot} argument) Line color when level of \code{b} is low.
Defaults to bright red (#EE4B2B").}

\item{line_color_high}{(\code{ggplot} argument) Line color when level of \code{b} is high.
Defaults to dark blue ("#00AFBB").}

\item{line_size}{(\code{ggplot} argument) Line widths for plot. Defaults to 0.8.}

\item{x_tick_break}{(\code{ggplot} argument) Tick mark breaks for x-axis}

\item{x_tick_label}{(\code{ggplot} argument) Tick mark labels for x-axis. If \code{x_tick_label} is
specified, the length must be the same as \code{x_tick_break}.}

\item{x_limit}{(\code{ggplot} argument) Upper and lower limits of x-axis}
}
\value{
When \code{int = "2-way"}, \code{gee_plot} returns a non-faceted plot; when
\code{int = "3-way"}, \code{gee_plot} returns a faceted plot consisting of 3
sub-plots; and when \code{int = "4-way"}, \code{gee_plot} returns a faceted plot
consisting of 9 sub-plots.
}
\description{
\code{gee_plot} wraps around (and combines) \code{geeglm} and \code{ggplot} to provide
visualizations for 2-, 3-, and 4-way interactions. For 3- and 4-way
interactions, the output is a faceted plot.
}
\details{
Since geeglm only works for complete data, NAs are removed via the following
specification in \code{geeglm} (\code{data = na.omit(data)}).
}
\examples{
gee_plot(data = df1c,
         dv = "MEA_c", dv_label = "Meaning",
         a = "POS", a_label = "Positive Affect",
         b = "NEG", b_label = "Negative Affect",
         "EXT", "NEU", "tMEA", "OPT", "SubHea",
         corstr = "exchangeable",
         int = "2-way",
         x_tick_break = c(2.9, 3.8, 4.7),
         x_tick_label = c("low", "mid", "high"),
         x_limit = c(1, 5))

}
\seealso{
\code{\link{geeglm}}
}
\author{
Jun Sheng Keh, \email{junshengkeh@gmail.com}
}
