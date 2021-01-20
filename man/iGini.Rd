\name{iGini}
\alias{iGini}
\title{Individual decomposition of the Gini Index}
\description{The function computes individual components of the Gini index
and their group-based decompositions. It takes as input an outcome variable, a grouping variable, and an optional sampling weight. It returns a data matrix of three columns containg individual contributions and their between- and within-group components. Because iGini is computational insensive, parallel processing is recommended, and the number of cores can be specified.}
\usage{
iGini(x, g, w=rep(1,length(x)),core=1)
}
\arguments{
\item{x}{Input continuous variable such as income.
}

\item{g}{A grouping variable containing integers, such gender coded 1 & 2.
}

\item{w}{An optional sampling weight variable.
}

\item{core}{An optional input for specifying the number of processing cores in your computer. When specified, you will need to have the \pkg{doParallel} package and the \pkg{foreach} package installed for conducting parallel processing to speed up the computation.
}
}

\value{The function outputs three variables, \var{g.i}, \var{g.ikb},
and \var{g.ikw}.

\item{g.i}{This variable gives the individual contribtions to the overall Gini index.
}
\item{g.ikb}{This variable provides for each individual component of the Gini its between-group subcomponent.
}
\item{g.ikw}{This variable provides for each individual component of the Gini its within-group subcomponent. The g.ikb and g.ikw sum up to g.i for each i observation.
}
}

\examples{
data(ChickWeight)
attach(ChickWeight)
iGini.result <- iGini(weight,Diet,core=1)
}

\references{
Tim F. Liao. 2019. "Individual Components of Three Inequality Measures for Analyzing Shapes of Inequality.” \emph{Sociological Methods & Research} Advance online publication. doi:10.1177/0049124119875961
}
