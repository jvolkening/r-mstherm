library(mstherm)
context("TPP modeling")

# Find demo project files
control <- system.file("extdata", "demo_project/control.tsv", package="mstherm")
annots  <- system.file("extdata", "demo_project/annots.tsv",  package="mstherm")

test_that("demo files are available", {
    expect_match(basename(control), "control.tsv")
    expect_match(basename(annots),  "annots.tsv")
})

# Build experiment
expt <- MSThermExperiment(control, annotations=annots)

test_that("MSThermExperiment creation", {
    expect_is(expt, "MSThermExperiment")
    expect_equal(expt$samples$Control$replicates$C1$data$TMT10.126[6],
        4568753, tolerance=1)
})

# Perform typical normalization
norm <- normalize_to_std(expt, "cRAP_ALBU_BOVIN")

test_that("MSThermExperiment normalization", {
    expect_is(norm, "MSThermExperiment")
    expect_equal(norm$samples$Control$replicates$C1$data$TMT10.126[6],
        4568561, tolerance=1)
})

# Perform typical modeling
res1 <- model_experiment(norm, bootstrap=T, smooth=T, min_rep_psm=0, np=2)
res2 <- model_experiment(norm, bootstrap=T, smooth=T, min_rep_psm=3, np=2)
res3 <- model_experiment(norm, bootstrap=T, smooth=F, min_rep_psm=3, np=2)
sgl1 <- res1$P38707
sgl2 <- res1$cRAP_ALBU_BOVIN


test_that("MSThermExperiment modeling", {
    expect_is(res1, "MSThermResultSet")
    expect_is(res2, "MSThermResultSet")
    expect_is(res3, "MSThermResultSet")
    expect_is(sgl1,  "MSThermResult")
    expect_is(sgl2,  "MSThermResult")

    expect_equal(length(res1), 6)
    expect_equal(length(res2), 5)
    expect_equal(length(res3), 5)

    # protein that should have modeled well
    expect_match(sgl1$annotation, "Asparagine--tRNA ligase, cytoplasmic")
    expect_true(sgl1$series$C1$is.fitted)

    expect_equal(sgl1$series$C1$tm,    48.0,   tolerance=0.1)
    expect_equal(sgl1$series$C1$k,    842,     tolerance=2)
    expect_equal(sgl1$series$C1$plat,   0.06,  tolerance=0.01)
    expect_equal(sgl1$series$C1$slope, -0.090, tolerance=0.001)
    expect_equal(sgl1$series$C1$r2,     0.99,  tolerance=0.01)
    expect_equal(sgl1$series$C1$rmsd,   0.03,  tolerance=0.01)
    expect_equal(sgl1$series$C1$inf,    0.16,  tolerance=0.01)
    expect_equal(sgl1$series$C1$psm,   46,     tolerance=0.01)

    expect_equal(sgl1$series$T1$tm,    52.1,   tolerance=0.1)
    expect_equal(sgl1$series$T1$k,    988,     tolerance=2)
    expect_equal(sgl1$series$T1$plat,   0.02,  tolerance=0.01)
    expect_equal(sgl1$series$T1$slope, -0.090, tolerance=0.001)
    expect_equal(sgl1$series$T1$r2,     0.99,  tolerance=0.01)
    expect_equal(sgl1$series$T1$rmsd,   0.04,  tolerance=0.01)
    expect_equal(sgl1$series$T1$inf,    0.20,  tolerance=0.01)
    expect_equal(sgl1$series$T1$psm,   92,     tolerance=0.01)

    # protein that should not have modeled at all
    expect_false(sgl2$series$C1$is.fitted)
    expect_equal(sgl2$series$C1$tm, NA)

})
