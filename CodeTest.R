
devtools::document()

# Test functions
usethis::use_testthat(3)

#################################################################.
##                     Data initialization                   ####
#################################################################.

## Load packages
library(gee.interactions)

## Global settings
options(scipen = 100)
options(digits = 5)

##---------------------------------------------------------------.
##                       Data processing                     ----
##---------------------------------------------------------------.
## Load raw data ----
load(file = "C:\\Users\\junsh\\Dropbox\\R\\01b_3way.RData")

# ------------------------------------------------------------------.
#__ Step 1: Compute fixed slope values ----
df1a <- df %>%
  select(SubID, MEA, POS, NEG,
         EXT, NEU, tMEA, OPT, SubHea) %>%
  mutate(MEA_c = MEA - mean(MEA, na.rm = TRUE))

#__ Step 2: Test for interactions ----
mod <- geeglm(MEA_c ~ (POS * NEG) +
                EXT + NEU + tMEA + OPT + SubHea,
              data = na.omit(df1a),
              id = SubID,
              corstr = "exchangeable"); summary(mod)

#__ Step 3: Plot simple slope plots ----
gee_plot(data = df1a,
         dv = "MEA_c", dv_label = "Meaning",
         a = "POS", a_label = "Positive Affect",
         b = "NEG", b_label = "Negative Affect",
         "EXT", "NEU", "tMEA", "OPT", "SubHea",
         corstr = "exchangeable",
         int = "2-way",
         # line_color_low = "#EE4B2B",
         # line_color_high = "#00AFBB",
         # line_size = 0.8,
         x_tick_break = c(2.9, 3.8, 4.7),
         x_tick_label = c("low", "mid", "high"),
         x_limit = c(1, 5))

# ------------------------------------------------------------------.
#__ Step 1: Compute fixed slope values ----
df1b <- df %>%
  select(SubID, MEA, POS, NEG, AGF,
         EXT, NEU, tMEA, OPT, SubHea) %>%
  mutate(MEA_c = MEA - mean(MEA, na.rm = TRUE))

#__ Step 2: Test for interactions ----
mod <- geeglm(MEA_c ~ (POS * NEG * AGF) +
                EXT + NEU + tMEA + OPT + SubHea,
              data = na.omit(df1b),
              id = SubID,
              corstr = "exchangeable"); summary(mod)

#__ Step 3: Plot simple slope plots ----
gee_plot(data = df1b,
         dv = "MEA_c", dv_label = "Meaning",
         a = "POS", a_label = "Positive Affect",
         b = "AGF", b_label = "Activity Goal Focus",
         c = "NEG", c_label = "Negative Affect",
         "EXT", "NEU", "tMEA", "OPT", "SubHea",
         corstr = "exchangeable",
         int = "3-way",
         line_color_low = "#00AFBB",
         line_color_high = "#EE4B2B"
         # line_size = 0.8,
         # x_tick_break = c(2.9, 3.8, 4.7),
         # x_tick_label = c("low", "mid", "high"),
         # x_limit = c(1, 5)
         )

# ------------------------------------------------------------------.
#__ Step 1: Compute fixed slope values ----
df1c <- df %>%
  select(SubID, MEA, POS, NEG, AGF, SOC,
         EXT, NEU, tMEA, OPT, SubHea) %>%
  mutate(MEA_c = MEA - mean(MEA, na.rm = TRUE))

#__ Step 2: Test for interactions ----
mod <- geeglm(MEA_c ~ (POS * NEG * AGF * SOC) +
                EXT + NEU + tMEA + OPT + SubHea,
              data = na.omit(df1c),
              id = SubID,
              corstr = "exchangeable"); summary(mod)

#__ Step 3: Plot simple slope plots ----
gee_plot(data = df1c,
         dv = "MEA_c", dv_label = "Meaning",
         a = "POS", a_label = "Positive Affect",
         b = "AGF", b_label = "Activity Goal Focus",
         c = "NEG", c_label = "Negative Affect",
         d = "SOC", d_label = "Social Interaction",
         "EXT", "NEU", "tMEA", "OPT", "SubHea",
         corstr = "exchangeable",
         int = "4-way",
         # line_color_low = "#EE4B2B",
         # line_color_high = "#00AFBB",
         # line_size = 0.8,
         x_tick_break = c(2.9, 3.8, 4.7),
         x_tick_label = c("low", "mid", "high"),
         x_limit = c(1, 5)
         )
