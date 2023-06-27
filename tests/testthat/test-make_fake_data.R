

test_that("make fake data works", {

    dd = make_fake_data()
    ddA = make_fake_data(scenario = "A")

    expect_true( nrow( dd ) == 18 )
    expect_true( mean( ddA$Z ) == 0.5 )

})
