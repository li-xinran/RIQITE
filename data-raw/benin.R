
# Code to clean and write original csv file.
if ( FALSE ) {
    library( tidyverse )
    benin <- read.csv("data/Benin.csv", stringsAsFactors = FALSE )
    benin <- benin %>% rename(VoteShare = vote_pop ) %>%
        mutate( Treatment = factor(Treatment,
                                   labels = c("Control", "Clientelist", "Policy")),
                Treatment = reorder(Treatment, -VoteShare))
    head( benin )

    save( benin, file = here::here( "data/benin.RData"  ) )
}
