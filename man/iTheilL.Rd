\name{iTheilL}
\alias{iTheilL}
\title{Individual decomposition of Theil's L Index}
\description{The function computes individual components of Theil's L index (or Theil's second measure) and their group-based decompositions. It takes as input an outcome variable, a grouping variable, and an optional sampling weight. It returns a data matrix of three columns containg individual contributions and their between- and within-group components.}
\usage{
iTheilL(x, g, w=rep(1,length(x)))
}
\arguments{
\item{x}{Input continuous variable such as income.
}

\item{g}{A grouping variable containing integers, such gender coded 1 & 2.
}

\item{w}{An optional sampling weight variable.
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
iTheilL.result <- iTheilL(weight,Diet)
}

\references{
Tim F. Liao. Forthcoming. "Individual Components of Three Inequality Measures for Analyzing Shapes of Inequality.” \emph{Sociological Methods & Research} 50:xxx–xxx. doi: https://doi.org/10.1177/0049124119875961
}
