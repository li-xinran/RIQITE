# Documentation of package and data



#' Benin Data
#'
#' These are the small Benin dataset, taken from Table III.H of Watchekon (2002).
#'
#' @format A data frame with 24 rows and 7 variables:
#' \describe{
#'   \item{Village}{Village ID}
#'   \item{District}{District of village}libra
#'   \item{Candidate}{Candidate in village}
#'   \item{reg_voters}{Number of registered voters}
#'   \item{VoteShare}{}
#'   \item{vote_male}{}
#'   \item{vote_male96}{}
#'   \item{vote_female}{}
#'   \item{vote_female96}{}
#'   \item{Treatment}{Is Control, Clientelist, or Policy}
#'   \item{Note}{Notes on our entering data from original paper.}
#' }
"benin"


# Code to clean and write original csv file.
if ( FALSE ) {
    library( tidyverse )
    benin <- read.csv("Examples/ExampleData/Benin.csv", stringsAsFactors = FALSE )
    benin <- benin %>% rename(VoteShare = vote_pop ) %>%
        mutate( Treatment = factor(Treatment,
                              labels = c("Control", "Clientelist", "Policy")),
           Treatment = reorder(Treatment, -VoteShare))
    head( benin )

    save( benin, file = "perminf/data/benin.RData" )
}