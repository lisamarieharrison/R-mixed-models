library('RUnit')

pathnames <- list.files(pattern="[.]R$", path="C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-functions-southern-ocean/", full.names=TRUE)
invisible(sapply(pathnames, FUN = source))

test.suite <- defineTestSuite("southernOcean",
                              dirs = file.path("C:/Users/Lisa/Documents/phd/southern ocean/Mixed models/R code/R-mixed-models/"),
                              testFileRegexp = 'test_southern_ocean')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)
