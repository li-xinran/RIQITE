
# Reconstruct individual-level data for the birth-month / school-entry-age
# instrumental-variable illustration (see demo_code/RIllustration20221028.Rmd
# and the main paper). The data is NOT a set of real individual records: it
# is built to match the aggregate proportions reported in Black et al.
# (2011), Table 1 (month of birth x school entry timing) and Table 3 (total
# size of the "discontinuity subsample" of subjects born in January or
# December: 104,023, rounded here to 104,000, split evenly between the two
# birth months following Fiorini and Stevens (2014)).

library(dplyr)

x <- matrix(0, 2, 3, dimnames = list(Month = c("Dec", "Jan"),
                                      Entry = c("Early", "On-Time", "Late")))
n <- 104000
x["Dec", "On-Time"] <- .85
x["Dec", "Late"]     <- .15
x["Jan", "Early"]    <- .10
x["Jan", "On-Time"]  <- .90

birth_school_age <- data.frame(
    BirthMonth  = rep(c("Dec", "Jan"), each = n / 2),
    EntryTiming = c(rep("On-Time", n / 2 * .85), rep("Late", n / 2 * .15),
                    rep("Early", n / 2 * .10), rep("On-Time", n / 2 * .90))
)

birth_school_age <- birth_school_age %>%
    mutate(EntryAge = NA,
           EntryAge = ifelse(BirthMonth == "Dec" & EntryTiming == "On-Time", 6.75, EntryAge),
           EntryAge = ifelse(BirthMonth == "Dec" & EntryTiming == "Late", 7.75, EntryAge),
           EntryAge = ifelse(BirthMonth == "Jan" & EntryTiming == "On-Time", 23 / 3, EntryAge),
           EntryAge = ifelse(BirthMonth == "Jan" & EntryTiming == "Early", 20 / 3, EntryAge))

table(birth_school_age$BirthMonth, birth_school_age$EntryTiming)

usethis::use_data(birth_school_age, overwrite = TRUE)
