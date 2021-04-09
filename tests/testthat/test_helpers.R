
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



# replace_symbols() -------------------------------------------------------

test_that("replace_symbols()", {
    expect_equal(
        replace_symbols(
            what = "before",
            by = "after",
            x = "this ~ before"
        ),
        "this~after"
    )

    expect_equal(
        replace_symbols(
            what = c("before1", "before2", "before3"),
            by = c("after1", "after2", "after3"),
            x = "this ~ before1*before2/before3"
        ),
        "this~after1*after2/after3"
    )
})
