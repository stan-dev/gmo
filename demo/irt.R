#!/usr/bin/env Rscript
# The data consists of the following.
#+ econ_attitudes.txt
#  + Responses to 10 econ attitude questions. Lots of missingness:
#    these correspond to times when the questions were not asked (it was a
#    multi-waved surves). All questions are signed so that + is the
#    conservative direction and - is the liberal direction.
#+ soc_attitudes.txt
#  + Same deal, 13 social attitude questions.
#+ econsoc.subset.txt
#  + It's a subset of the survey questions, not a subset of respondents. Columns:
#  + e.index: a variable (I don't remember how we created it)
#    representing econ ideology, somehow it's a combination of those 10
#    variables. We can ignore it.
#  + s.index: same thing, social ideology, we also can ignore.
#  + bush: 1 for Bush supporters, 0 for Gore supporters
#  + st.num: index for the state of residence
#  + CST: state abbreviation
#  + pid: Party id (1-7 scale with 1=strong Democrat and 7=strong Republican)
#  + ideol: ideology (1-5 scale with 1=strong liberal and 5=strong conservative)
#  + relatt: religious attendance (1-5 scale with 1=never and 5=more than once/week)
econ_attitudes <- read.table("data/econ_attitudes.txt", header=T)
soc_attitudes <- read.table("data/soc_attitudes.txt", header=T)
econsoc <- read.table("data/econsoc.subset.txt", header=T)

# only use 'bush' and 'st.num'
econsoc <- econsoc[, c("bush", "st.num")]
# 'bush' is the response in {0, 1}

form_data <- function(econ_attitudes, soc_attitudes, econsoc) {
  data <- list(
    bush=econsoc$bush,
    st.num=econsoc$st.num
  )
  return(data)
}

# TODO write stan programs
data <- form_data(econ_attitudes, soc_attitudes, econsoc)
fit.gmo <- gmo("models/ideal_point.stan", "models/ideal_point_local.stan", data=data)
fit.nuts <- stan("models/ideal_point_original.stan", data=data)

# Fit to a subset of data, where NUTS will converge.
N <- 1000
econ_attitudes <- econ_attitudes[1:N, ]
soc_attitudes <- soc_attitudes[1:N, ]
econsoc <- econsoc[1:N, ]

data <- form_data(econ_attitudes, soc_attitudes, econsoc)
fit.gmo <- gmo("models/ideal_point.stan", "models/ideal_point_local.stan", data=data)
fit.nuts <- stan("models/ideal_point_original.stan", data=data)
