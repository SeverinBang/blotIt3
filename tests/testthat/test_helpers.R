
# get_symbols() -----------------------------------------------------------

test_that("read_symbols()", {
    expect_equal(
        get_symbols("left1 ~ right11 + right12"),
        c("left1", "right11", "right12")
    )
    expect_equal(
        get_symbols(
            "left2 ~ (right21 * right22 / (right23 * right24)) - right25"
        ),
        c("left2", "right21", "right22", "right23", "right24", "right25")
    )
})



# identify_parameters() ---------------------------------------------------

test_that("identify_parameters()", {
    expect_error(
        identify_effects(
            distinguish = left ~ right,
            scaling = "as ~ string",
            error = left ~ right
        ),
        "Do not pass distinguish, scaling or error as string."
    )


    expect_error(
        identify_effects(
            scaling = left ~ right,
            error = left ~ right
        ),
        "Left and right-hand side of formula 'distinguish' is needed"
    )

    expect_error(
        identify_effects(
            distinguish = left ~ right + right1,
            error = left ~ right
        ),
        "Left and right-hand side of formula 'scaling' is needed"
    )


    expect_error(
        identify_effects(
            distinguish = left ~ right + right1,
            scaling = one ~ two
        ),
        "Left and right-hand side of formula 'error' is needed"
    )
})
