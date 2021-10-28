##
## This script should be finished to look at how long it takes to do testing
## (e.g., run the exact test) when we have very large datasets.
##
## Script is unfinished
##


benin_ci_precise.steph <- find_ci_stephenson( benin_sub,
                                              tr_var = "Treatment",
                                              tr_level = "Policy",
                                              y_var = "VoteShare",
                                              strat_var = "District",
                                              alternative = "less",
                                              verbosity = 0,
                                              distribution = approximate(100000),
                                              subset_size = 6 )

benin_ci_precise.steph$ci


benin_ci_precise.steph <- find_ci_stephenson( benin_sub,
                                              tr_var = "Treatment",
                                              y_var = "VoteShare",
                                              strat_var = "District",
                                              alternative = "less",
                                              verbosity = 0,
                                              distribution = approximate(100000),
                                              subset_size = 6 )

benin_ci_precise.steph$ci


benin_ci_precise.steph <- find_ci_stephenson( benin_sub,
                                              tr_var = "Treatment",
                                              y_var = "VoteShare",
                                              alternative = "less",
                                              verbosity = 0,
                                              distribution = approximate(100000),
                                              subset_size = 6 )

benin_ci_precise.steph$ci


benin_ci_precise.steph <- find_ci_stephenson( benin_sub,
                                              tr_var = "Treatment",
                                              y_var = "VoteShare",
                                              alternative = "less",
                                              verbosity = 0,
                                              distribution = "exact",
                                              subset_size = 6 )

benin_ci_precise.steph$ci


benin_sub$Z = as.numeric( benin_sub$Treatment == "Policy" )

benin_ci_precise.steph <- find_ci_stephenson( benin_sub,
                                              tr_var = "Z",
                                              y_var = "VoteShare",
                                              alternative = "less",
                                              verbosity = 0,
                                              distribution = "exact",
                                              subset_size = 6 )

benin_ci_precise.steph$ci



dat = make.obs.data( 50, p=0.5, p.extreme=0.05, tau=0.4, DGP=make.data.mix )
dat$blk = sample( LETTERS[1:3], replace=TRUE, nrow( dat ))
table( dat$blk, dat$Z )

benin_ci_precise.steph <- find_ci_stephenson( dat,
                                              tr_var = "Z",
                                              y_var = "Yobs",
                                              alternative = "less",
                                              verbosity = 0,
                                              distribution = "exact",
                                              subset_size = 6 )

benin_ci_precise.steph$ci
benin_ci_precise.steph$pvalues





dat = make.obs.data( 100, p=0.5, p.extreme=0.05, tau=0.4, DGP=make.data.mix )
dat$blk = sample( LETTERS[1:3], replace=TRUE, nrow( dat ))
table( dat$blk, dat$Z )

benin_ci_precise.steph <- find_ci_stephenson( dat,
                                              tr_var = "Z",
                                              y_var = "Yobs",
                                              alternative = "less",
                                              verbosity = 2,
                                              distribution =  approximate(100000),
                                              subset_size = 6 )

benin_ci_precise.steph$ci
benin_ci_precise.steph$pvalues


