\name{crimCV}
\alias{crimCV}

\title{
  Fits finite mixtures of Zero-inflated Poisson models
}
\description{
  This software fits finite mixtures of ZIP models to longitudinal count data. 
}
\usage{
crimCV(Dat,ng,dpolyp=3,dpolyl=3,model="ZIPt",rcv=FALSE,init=20,Risk=NULL)
}

\arguments{
  \item{Dat}{
    A matrix of the number of criminal offenses.  Each row contains the
    number of criminal offenses committed by a particular individual for
    all time intervals observed (columns).  A negative number is interpreted as a
    missing value.
}
  \item{ng}{
    The number of groups to use in the analysis.
}
  \item{dpolyp}{
    The degree of the polynomial in the log-linear component.  Defaults
    to 3.
}
  \item{dpolyl}{
    The degree of the polynomial in the logistic component (ZIP model
    only).  Defaults to 3.
}
  \item{model}{
    Either "ZIP" for the full ZIP  or "ZIPt" for the
    ZIP(tau) sub-group model.  Defaults to "ZIPt".
}
  \item{rcv}{
    Set to TRUE to compute the CVE.  Defaults to FALSE.
}
  \item{init}{
    This controls how rigorously the initialization procedure searches
    for good starting values of the parameters.  Larger values increase the odds of
    finding the true global solution but can dramatically increase the
    time required to fit the model.  The default is set to 20 and from
    experimentation this seems to be a reasonable value for ng<=5.
    Note: As the number of groups (ng) gets larger finding the optimal
    global maximum of the likelihood becomes progressively more
    challenging.
}
  \item{Risk}{
    A matrix of the same dimension as Dat with the time-at-risk
    correction.  Defaults to 1 for all elements (i.e. at risk for the
    full time period).
}
}
\value{
  A list of class "dmZIP" or "dmZIPt".
}
\references{
J.D. Nielsen, J.S. Rosenthal, Y. Sun, D.M. Day, I. Bevc, and
  T. Duchesne (2011). Group-based Criminal Trajectory Analysis using
  Cross-Validation Criteria.  A draft of the manuscript is available at
  \emph{http://www.probability.ca/jeff/research.html}.
}
\author{
  Jason D. Nielsen
}

\examples{
# Loads crimCV into the interpreter
library(crimCV)

# Load the "divide-and-round" TO1 dataset
data(TO1adj)
# Fit a 2 component ZIP(tau) model of degree 2.  Here the CVE is not
# calculated and only ~1/4 of the data is used so that the code will
# run quickly enough to satisfy CRAN's package policies.  To compute
# CVE run as:
# out1<-crimCV(TO1adj,2,dpolyp=2,rcv=TRUE)
subTO1adj<-TO1adj[1:100,]
out1<-crimCV(subTO1adj,2,dpolyp=2,init=5)
# Plot the component trajectories
plot(out1)
# Print out some useful output
summary(out1)
}

\keyword{ Zero-inflated Poisson model}
\keyword{ Finite mixture model }


