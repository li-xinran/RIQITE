


# Regression tests for li-xinran/RIQITE#1 (reported by jwbowers):
#
# min_stat() and null_dist() sum the same rank scores in different orders.
# Floating point addition is not associative, so two sums that are equal
# in exact arithmetic can land on different doubles (differing in the
# last bit or two). A strict `stat.null >= stat.min` comparison then fails
# to count a null draw that is mathematically tied with the observed
# statistic, undercounting the p-value.

test_that("pval_from_null treats summation-order floating point ties as equal", {

    # Construct two sums over the *same multiset* of rank scores, added in
    # different orders, to reproduce the exact failure mode from the bug
    # report: deterministic (no RNG), portable IEEE-754 arithmetic.
    n = 50
    score = (1:n) / n
    idx = c(5, 12, 23, 31, 42, 7, 18, 29, 36, 44, 2, 14, 25, 33, 48,
            9, 20, 27, 39, 46, 3, 16, 22, 35, 50)

    stat_a = sum( score[idx] )
    stat_b = sum( score[ sort(idx) ] )

    # If a platform ever makes these identical, the bug can't manifest and
    # this test is moot rather than meaningful.
    skip_if( identical( stat_a, stat_b ),
             "platform did not reproduce the summation-order floating point tie" )

    # Make stat_a the larger representation (plays the role of stat.min /
    # stat.obs) and stat_b the smaller one that must still be counted as tied.
    if ( stat_a < stat_b ) {
        tmp = stat_a; stat_a = stat_b; stat_b = tmp
    }

    # A naive strict comparison undercounts the tie.
    expect_false( stat_b >= stat_a )

    # pval_from_null() must count it.
    stat.null = c( stat_b, 0 )
    expect_equal( RIQITE:::pval_from_null( stat.null, stat_a ), 0.5 )
})


test_that("pval_quantile reproduces the worked example from li-xinran/RIQITE#1", {

    set.seed(54321)
    n = 50
    m = 25
    Z = sample( c( rep(1, m), rep(0, n - m) ) )
    Y = rnorm(n)

    set.seed(98765)
    Z.perm = assign_CRE(n, m, 1000)

    ml = list(name = "Wilcoxon")
    nul = null_dist(n, m, nperm = 1000, Z.perm = Z.perm, method.list = ml)

    pv = pval_quantile(Z, Y, k = 45, c = 0, method.list = ml,
                        stat.null = nul, switch = FALSE)

    expect_equal( pv, 0.772 )
})


test_that("conf_quant_larger (used by ci_quantile) is not corrupted by summation-order ties", {

    # Same worked example as li-xinran/RIQITE#1. conf_quant_larger() does its
    # own `stat.min(c) <= thres`-style comparisons during root-finding, built
    # from the same kind of rank-score sums as pval_H_k_c(), so it was
    # exposed to the same floating point tie mechanism -- just not the one
    # named in the original report. (Before the fix, 8 of the 50 confidence
    # limits in this example differed by up to 0.05.)
    set.seed(54321)
    n = 50
    m = 25
    Z = sample( c( rep(1, m), rep(0, n - m) ) )
    Y = rnorm(n)

    set.seed(98765)
    Z.perm = assign_CRE(n, m, 1000)

    ml = list(name = "Wilcoxon")

    # Integer scores are mathematically equivalent to the normalized scores
    # rank_score() actually produces (just scaled by n), but sums of small
    # integers are exact in double precision -- no tie risk at all. This
    # gives a trustworthy reference to compare the real (normalized, tie-
    # prone) scores against; the confidence limits should be identical.
    score_exact = as.numeric(1:n)
    score_float = score_exact / n

    c_exact = RIQITE:::conf_quant_larger(Z, Y, method.list = ml, score = score_exact,
                                          Z.perm = Z.perm, alpha = 0.05, tol = 1e-3)
    c_float = RIQITE:::conf_quant_larger(Z, Y, method.list = ml, score = score_float,
                                          Z.perm = Z.perm, alpha = 0.05, tol = 1e-3)

    expect_equal( c_exact, c_float )
})
