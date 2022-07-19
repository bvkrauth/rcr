test_that("install_rcrpy produces expected results", {
  skip("not yet tested")
})

test_that("rcr works with default options", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  true_coef <- c(
    12.310599093,
    8.169709965,
    28.935489170,
    5.135043765,
    5.201502574
  )
  true_coef_names <- c(
    "rcInf",
    "effectInf",
    "rc0",
    "effectL",
    "effectH"
  )
  true_cov <- c(
    4.40273105084336, 1.6809105662234, 14.8603396952841,
    0.0262163549131941, 0.0148105699433406, 1.68091056622339,
    936.816073751491, -3305.54493690222, -20.8604783714285,
    0.0945995702080957, 14.8603396952841, -3305.54493690222,
    11776.4762834738, 76.3213527521248, 2.09329547536395,
    0.0262163549131933, -20.8604783714285, 76.3213527521247,
    0.915729395515209, 0.438565220757964, 0.0148105699433403,
    0.094599570208096, 2.09329547536393, 0.438565220757964,
    0.430902711418667
  )
  true_model_dim <- c(5839, 8)
  result <- rcr(f1,
    data = testdata
  )
  # class
  expect_s3_class(result,
    "rcr",
    exact = TRUE
  )
  # result$coefficients
  expect_equal(
    as.vector(result$coefficients),
    true_coef
  )
  expect_equal(
    as.vector(result$coefficients),
    as.vector(result$pyobj$params)
  )
  expect_equal(
    names(result$coefficients),
    true_coef_names
  )
  # result$cov.unscaled
  expect_equal(
    colnames(result$cov.unscaled),
    true_coef_names
  )
  expect_equal(
    rownames(result$cov.unscaled),
    true_coef_names
  )
  expect_equal(
    as.vector(result$cov.unscaled),
    true_cov
  )
  expect_equal(
    as.vector(result$cov.unscaled),
    as.vector(result$pyobj$cov_params)
  )
  # result$xlevels
  # result$call
  expect_type(
    result$call,
    "language"
  )
  # result$terms
  expect_s3_class(
    result$terms,
    "terms"
  )
  expect_equal(
    dim(result$model),
    true_model_dim
  )
  expect_true(is.null(result$weights))
  expect_true(is.null(result$cluster))
})

test_that("rcr allows functions in formulas", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- abs(SAT) ~ abs(Small_Class) | White_Asian * Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  true_coef <- c(
    13.8370293027089,
    8.16829969322228,
    38.8414928511089,
    5.49447815019125,
    5.54093792567284
  )
  result <- rcr(f1,
    data = testdata
  )
  # result$coefficients
  expect_equal(
    as.vector(result$coefficients),
    true_coef
  )
})

test_that("rcr does not require data argument", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  true_coef <- coef(rcr(f1,
    data = testdata
  ))
  result <- rcr(testdata$SAT ~ testdata$Small_Class | testdata$White_Asian +
    testdata$Girl + testdata$Free_Lunch +
    testdata$White_Teacher + testdata$Teacher_Experience +
    testdata$Masters_Degree)
  # result$coefficients
  expect_equal(
    as.vector(result$coefficients),
    as.vector(true_coef)
  )
})


test_that("rcr works with subset option", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  true_result <- rcr(f1, testdata[1:1000, ])
  result <- rcr(f1, testdata, subset = 1:1000)
  expect_equal(
    coef(result),
    coef(true_result)
  )
  expect_equal(
    nrow(result$model),
    1000
  )
})

test_that("rcr works with model and pyobj options", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  result <- rcr(f1, testdata, model = FALSE, pyobj = FALSE)
  expect_true(is.null(result$model))
  expect_true(is.null(result$pyobj))
})


test_that("rcr works with weights option", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  wt <- testdata$TCHID %% 2
  true_result <- rcr(f1, testdata[wt > 0.5, ])
  result <- rcr(f1, testdata, weights = wt)
  expect_equal(
    coef(result),
    coef(true_result)
  )
  expect_equal(
    result$weights,
    wt
  )
  expect_equal(
    result$model$`(weights)`,
    wt
  )
})

test_that("rcr works with na.action option", {
  testdata <- readRDS(test_path("testdata.rds"))
  altdata <- testdata
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  obs <- (testdata$TCHID %% 2) > 0.5
  true_result <- rcr(
    f1,
    testdata[obs, ]
  )
  altdata$SAT[!obs] <- NA
  result <- rcr(f1,
    altdata,
    na.action = na.omit
  )
  expect_equal(
    coef(result),
    coef(true_result)
  )
  expect_error(
    rcr(f1,
      altdata,
      na.action = na.fail
    ),
    "missing values in object"
  )
})

test_that("rcr works with cluster option", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  true_cov <- c(
    69.168147160016, 126.630806365959, -121.5489537968,
    -1.94406588069719, 0.122940162326544, 126.630806365959,
    1904.95751929006, -6125.90888676137, -37.5141026937444,
    3.6638148573893, -121.548953796798, -6125.90888676138,
    21076.4577513673, 129.096256503117, -6.60769825544483,
    -1.94406588069717, -37.5141026937447, 129.096256503117,
    1.84604808673233, 1.00506873788518, 0.122940162326546,
    3.66381485738799, -6.60769825544113, 1.00506873788517,
    1.0621294617475
  )
  base_result <- rcr(f1, testdata)
  result <- rcr(f1, testdata,
    cluster = testdata$TCHID
  )
  expect_equal(
    coef(result),
    coef(base_result)
  )
  expect_equal(
    as.vector(vcov(result)),
    true_cov
  )
  expect_equal(
    result$cluster,
    testdata$TCHID
  )
  expect_equal(
    result$model$`(cluster)`,
    testdata$TCHID
  )
})

test_that("rcr works with rc_range option", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  true_coef <- c(
    12.310599093,
    8.169709965,
    28.935489170,
    5.135043765,
    5.201502574
  )
  result <- rcr(f1,
    data = testdata,
    rc_range = c(0, 0)
  )
  expect_equal(
    as.vector(result$coefficients[4:5]),
    rep(true_coef[5], 2)
  )
  result <- rcr(f1,
    data = testdata,
    rc_range = c(0, 15)
  )
  expect_equal(
    as.vector(result$coefficients[4:5]),
    c(-Inf, Inf)
  )
})

test_that("rcr works with vceadj option", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  true_cov <- c(
    4.40273105084336, 1.6809105662234, 14.8603396952841,
    0.0262163549131941, 0.0148105699433406, 1.68091056622339,
    936.816073751491, -3305.54493690222, -20.8604783714285,
    0.0945995702080957, 14.8603396952841, -3305.54493690222,
    11776.4762834738, 76.3213527521248, 2.09329547536395,
    0.0262163549131933, -20.8604783714285, 76.3213527521247,
    0.915729395515209, 0.438565220757964, 0.0148105699433403,
    0.094599570208096, 2.09329547536393, 0.438565220757964,
    0.430902711418667
  )
  result <- rcr(f1,
    data = testdata,
    vceadj = 0.85
  )
  expect_equal(
    as.vector(vcov(result)),
    0.85 * as.vector(true_cov)
  )
})

test_that("print.rcr works with default options", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_output <- c(
    "",
    "Call:",
    "rcr(formula = f1, data = testdata)",
    "",
    "Coefficients:",
    "    rcInf effectInf       rc0   effectL   effectH ",
    "    12.31      8.17     28.94      5.14      5.20 ",
    ""
  )
  result <- capture.output(print(rcr1, digits = 3))
  expect_equal(
    result,
    true_output
  )
})


test_that("vcov method works with default options", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_cov <- rcr1$cov.unscaled
  result <- vcov(rcr1)
  expect_equal(
    result,
    true_cov
  )
})

test_that("confint method works with default options", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_ci <- t(cbind(
    rcr1$pyobj$params_ci(),
    rcr1$pyobj$effect_ci()
  ))
  result <- confint(rcr1)
  expect_equal(
    rownames(result),
    c(names(rcr1$coefficients), "effect")
  )
  expect_equal(
    colnames(result),
    c("2.5  %", "97.5  %")
  )
  expect_equal(
    as.vector(result),
    as.vector(true_ci)
  )
})


test_that("confint method works with parm option", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_rc0 <- rcr1$pyobj$params_ci()[, 3]
  true_effect <- rcr1$pyobj$effect_ci()
  parm <- c("effect", "rc0", "junk", "effect")
  result <- confint(rcr1,
    parm = parm
  )
  expect_equal(
    rownames(result),
    parm
  )
  expect_equal(
    colnames(result),
    c("2.5  %", "97.5  %")
  )
  expect_equal(
    as.vector(result[1, ]),
    as.vector(true_effect)
  )
  expect_equal(
    as.vector(result[2, ]),
    as.vector(true_rc0)
  )
  expect_true(all(is.na(result[3, ])))
  expect_equal(
    as.vector(result[4, ]),
    as.vector(true_effect)
  )
})

test_that("confint method works with level option", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_ci <- t(cbind(
    rcr1$pyobj$params_ci(cilevel = 90),
    rcr1$pyobj$effect_ci(cilevel = 90)
  ))
  result <- confint(rcr1, level = 0.9)
  expect_equal(
    rownames(result),
    c(names(rcr1$coefficients), "effect")
  )
  expect_equal(
    colnames(result),
    c("5  %", "95  %")
  )
  expect_equal(
    as.vector(result),
    as.vector(true_ci)
  )
})


test_that("confint method works with citype lower", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_effect_ci <- rcr1$pyobj$effect_ci_lower()
  true_ci <- c(
    -Inf, -Inf, -Inf, -Inf, -Inf, -Inf,
    15.76194378392, 58.5144872118226, 207.434139924289,
    6.7090658968995, 6.28123680488214, 6.28123680488214
  )
  result <- confint(rcr1,
    citype = "lower"
  )
  expect_equal(
    as.vector(result),
    as.vector(true_ci)
  )
  expect_equal(
    as.vector(result[6, ]),
    as.vector(true_effect_ci)
  )
})

test_that("confint method works with citype upper", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_effect_ci <- rcr1$pyobj$effect_ci_upper()
  true_ci <- c(
    8.85925440230886, -42.1750672817734, -149.563161584069,
    3.56102163306907, 4.12176834229165, 3.56102163306907,
    Inf, Inf, Inf, Inf, Inf, Inf
  )
  result <- confint(rcr1,
    citype = "upper"
  )
  expect_equal(
    as.vector(result),
    as.vector(true_ci)
  )
  expect_equal(
    as.vector(result[6, ]),
    as.vector(true_effect_ci)
  )
})

test_that("confint method works with citype Imbens-Manksi", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_effect_ci <- rcr1$pyobj$effect_ci_imbensmanski()
  true_ci <- c(
    8.19806823846957, -51.8197921989432, -183.758771908664,
    3.25948071251956, 3.91491988225289, 3.29158005546998,
    16.4231299477593, 68.1592121289925, 241.629750248883,
    7.01060681744901, 6.4880852649209, 6.46606603235554
  )
  result <- confint(rcr1,
    citype = "Imbens-Manski"
  )
  expect_equal(
    as.vector(result),
    as.vector(true_ci)
  )
  expect_equal(
    as.vector(result[6, ]),
    as.vector(true_effect_ci)
  )
  # any other value should return error
  expect_error(
    confint(rcr1,
      citype = "nonexistent"
    ),
    "invalid citype nonexistent"
  )
})

test_that("effect_test method works", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_p <- rcr1$pyobj$test_effect()
  result <- effect_test(rcr1)
  expect_equal(result,
    true_p,
    tolerance = 0.001
  )
  true_p <- rcr1$pyobj$test_effect(3.29158)
  result <- effect_test(rcr1, h0 = 3.29158)
  expect_equal(result,
    true_p,
    tolerance = 0.001
  )
  true_p <- rcr1$pyobj$test_effect(5.15)
  result <- effect_test(rcr1, h0 = 5.15)
  expect_equal(result,
    true_p,
    tolerance = 0.001
  )
})


test_that("summary method works with default options", {
  testdata <- readRDS(test_path("testdata.rds"))
  f1 <- SAT ~ Small_Class | White_Asian + Girl +
    Free_Lunch + White_Teacher +
    Teacher_Experience + Masters_Degree
  rcr1 <- rcr(f1,
    data = testdata
  )
  true_ci <- rcr1$pyobj$effect_ci()
  true_coef <- coef(rcr1)
  true_ses <- sqrt(diag(vcov(rcr1)))
  true_tstat <- true_coef / true_ses
  true_pvalue <- 2 * pnorm(-abs(true_tstat))
  result <- summary(rcr1)
  expect_equal(
    result$call,
    rcr1$call
  )
  expect_equal(
    result$terms,
    rcr1$terms
  )
  # coefficients
  expect_equal(
    rownames(result$coefficients),
    names(coefficients(rcr1))
  )
  expect_equal(
    colnames(result$coefficients),
    c(
      "Estimate", "Std. Error",
      "t value", "Pr(>|t|)"
    )
  )
  expect_equal(
    as.vector(result$coefficients[, 1]),
    as.vector(true_coef)
  )
  expect_equal(
    as.vector(result$coefficients[, 2]),
    as.vector(true_ses)
  )
  expect_equal(
    as.vector(result$coefficients[, 3]),
    as.vector(true_tstat)
  )
  expect_equal(
    as.vector(result$coefficients[, 4]),
    as.vector(true_pvalue)
  )
  # cov.unscaled
  expect_equal(
    result$cov.unscaled,
    rcr1$cov.unscaled
  )
  # effect_ci
  expect_equal(
    rownames(result$effect_ci),
    "effect"
  )
  expect_equal(
    colnames(result$effect_ci),
    c("2.5  %", "97.5  %")
  )
  expect_equal(
    as.vector(result$effect_ci),
    as.vector(true_ci)
  )
  # citype
  expect_equal(result$citype, "conservative")
  # level
  expect_equal(result$level, 0.95)
})
