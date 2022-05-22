outer <- 0
galname <- "NGC 5846 (E0)"
extent <- 13.3
extent_err <- 1.0
distmod <- 31.98
distmod_err <- 0.2
label1 <- 'Mosaic'
label2 <- 'HST'
datapath <- "~/gcc/n5846/gcfinder3/"
outfile <- "n5846_radprof.pdf"
source('radprof.R',print.eval=TRUE)

galname <- "NGC 4649 (E2)"
datapath <- "~/gcc/n4649/gcfinder3.n4649/"
outfile <- "n4649_radprof.pdf"
extent <- 19.24
extent_err <- 1.0
distmod <- 31.13
distmod_err <- 0.4
label1 <- 'Mosaic'
label2 <- 'HST'
source('radprof.R',print.eval=TRUE)

galname <- "NGC 4621 (E5)"
datapath <- "~/gcc/n4649/gcfinder3.n4621/"
outfile <- "n4621_radprof.pdf"
extent <- 14.2
extent_err <- 0.6
distmod <- 31.31
distmod_err <- 0.2
label1 <- 'Mosaic'
label2 <- 'HST'
source('radprof.R',print.eval=TRUE)

galname <- "NGC 4382 (S0)"
datapath <- "~/gcc/n4382/gcfinder3/"
outfile <- "n4382_radprof.pdf"
extent <- 20.5
extent_err <- 0.8
distmod <- 31.33
distmod_err <- 0.14
label1 <- 'Mosaic'
label2 <- 'HST'
source('radprof.R',print.eval=TRUE)

galname <- "NGC 1023 (SB0)"
datapath <- "~/gcc/n1023/"
outfile <- "n1023_radprof.pdf"
extent <- 6.27
extent_err <- 0.4
distmod <- 30.29
distmod_err <- 0.16
label1 <- 'WIYN'
label2 <- 'HST'
source('radprof.R',print.eval=TRUE)

galname <- "NGC 7332 (S0pec)"
datapath <- "~/gcc/n7332/"
outfile <- "n7332_radprof.pdf"
extent <- 1.95
extent_err <- 0.25
distmod <- 31.81
distmod_err <- 0.20
outer <- 5.0
label1 <- 'WIYN'
label2 <- 'Keck'
source('radprof.R',print.eval=TRUE)


galname <- "NGC 7339 (Sbc)"
datapath <- "~/gcc/n7339/"
outfile <- "n7339_radprof.pdf"
extent <- 1.46
extent_err <- 0.15
distmod <- 31.75
distmod_err <- 0.37
label1 <- 'WIYN'
outer <- 3.5
source('radprof_single.R',print.eval=TRUE)


galname <- "NGC 4013 (Sb)"
datapath <- "~/gcc/n4013_gcfinder/"
outfile <- "n4013_radprof.pdf"
extent <- 2.78
extent_err <- 0.5
distmod <- 30.90
distmod_err <- 0.15
label1 <- 'WIYN'
label2 <- 'HST'
source('radprof.R',print.eval=TRUE)

