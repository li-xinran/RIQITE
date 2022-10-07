
library( tidyverse )
dat = read_csv( here::here( "demo_code/cleanTeacher.csv" ), na = "." )
dat

dat = dplyr::select( dat, T.ID, Site, Tx, Per.1, Per.2, gain )
dat = mutate( dat, TxAny = ifelse( Tx == "D", 0, 1 ) )
# remove units with missing outcomes and units in Site 1 that has no control outcomes
dat = filter( dat, !is.na( gain ) & Site != "S1" )
n = nrow(dat)
table(dat$TxAny)

#dat = sample_n( dat, nrow(dat) )

dat

electric_teachers = dat

usethis::use_data(dat, overwrite = TRUE)
