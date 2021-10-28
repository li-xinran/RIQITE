



#' Make fake dataset for testing and illustration
#'
#' These are all fixed set-size and -outcome datasets.  The pretend
#' randomization is perfectly balanced, so we have two units of each type with
#' one in treatment and one in control.
#'
#' The scenarios are:
#'
#' A: Long right tail, but also some very low treatment units. Clear evidence of
#' harm to some units and strong positive effects for others (the tx units are
#' in the tails, co units are interior)
#'
#' B: Treatment has lowest values, control has highest.  It could be a simple
#' shift down (constant negative tx effect visually plausible).
#'
#' C: Treatment helps in general, and helps some a lot.
#'
#' D: Treatment helps those at the low end, but hurts those at the high end. So
#' testing for positive tx (alt = greater) on left tail should reject. Testing
#' for negative tx on right tail(alt = less) should also reject.
#'
#' @param scenario   "A", "B", "C", or "D"; different kinds of data structure.
#'
#' @return Dataset with Y0, Y1, Yobs, and Z (treatment assignment).
#' @export
make_fake_data = function( scenario = "D" ) {


  if ( scenario == "A" ) {

    dat_co = data.frame( Z = "C",
                         Yobs = c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 8.5, 9, 10, 11, 12 ),
                         stringsAsFactors = FALSE)
    dat_tx = data.frame( Z = "T",
                         Yobs = c(1.5, 1, 2, 3, 4, 5, 6, 7, 8, 19, 30, 31, 32, 33) - 5.5,
                         stringsAsFactors = FALSE)

  } else if ( scenario == "B" ) {

    dat_co = data.frame( Z = "C",
                         Yobs = c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 19, 20, 21, 22 ),
                         stringsAsFactors = FALSE)
    dat_tx = data.frame( Z = "T",
                         Yobs = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 22) - 5.5,
                         stringsAsFactors = FALSE)

  } else if ( scenario == "C" ) {
    dat_co = data.frame( Z = "C",
                         Yobs = c( 0, 1, 2, 3, 4, 5, 8, 8.5, 9, 10, 11, 12 ),
                         stringsAsFactors = FALSE)
    dat_tx = data.frame( Z = "T",
                         Yobs = c(1.5, 1, 2, 3, 4, 5, 8, 19, 30, 31, 32,33),
                         stringsAsFactors = FALSE)

  } else {  # scenario D
    dat_co = data.frame( Z = "C",
                         Yobs = c( 0.5, 1, 1, 1.5,  8, 9, 10, 15, 17 ),
                         stringsAsFactors = FALSE)
    dat_tx = data.frame( Z = "T",
                         Yobs = c( 2, 4, 6, 6, 7, 4, 5, 6, 5 ),
                         stringsAsFactors = FALSE)
  }
  nrow( dat_co )
  nrow( dat_tx )
  stopifnot( nrow( dat_co ) == nrow( dat_tx ) )

  dat_co = arrange( dat_co, Yobs )
  dat_tx = arrange( dat_tx, Yobs )
  dat_co = dat_co %>% mutate( Y0 = Yobs,
                              Y1 = dat_tx$Yobs )
  dat_tx = dat_tx %>% mutate( Y0 = dat_co$Yobs,
                              Y1 = Yobs )

  dat = bind_rows( dat_co, dat_tx )
  dat = arrange( dat, Y0 )

  if ( FALSE ) {
    plot( dat$Y0, dat$Y1, asp=1 )
    abline( 0, 1 )
  }



  dat = mutate( dat,
                Z.f = factor( dat$Z, levels=c("T","C") ),
                Z = ifelse( Z.f == "T", 1, 0 ),
                Y0 = 2 * Y0,
                Y1 = 2 * Y1,
                Yobs = 2 * Yobs )

  dat = mutate( dat,
                tau = Y1 - Y0 ) %>%
    relocate( Y0, Y1, tau, Z, Z.f, Yobs)

  dat
}
