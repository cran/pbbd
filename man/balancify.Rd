\name{balancify}
\alias{balancify}
\title{Position balanced and nearly position balanced block design}
\description{This function generates a position balanced or nearly position balanced block design from a given equireplicate and proper block design}
\usage{balancify(d1)}
\arguments{
	\item{d1}{Block design specified in the form of a b x k matrix with elements labelled as 1 to v where b is number of blocks, k is block size and v is number of treatments}	
	}
\value{
	\item{design}{(Nearly) position balanced block design}
	\item{P}{Treatment by Position incidence matrix}
  }
\examples{
d1 = matrix(c(3,  4,    6,
5,    6,    7,
1,   4,    5,
2,    4,    7,
1,    3,    7,
1,    2,    6,
2,    3,    5), ncol = 3, byrow = TRUE)
balancify(d1)

d1 = matrix(c(7	,	8	,	9	,
              1	,	6	,	8	,
              1	,	3	,	9	,
              4	,	6	,	9	,
              5	,	6	,	7	,
              1	,	4	,	5	,
              3	,	5	,	8	,
              3	,	4	,	7	,
              2	,	5	,	9	,
              2	,	4	,	8	,
              1	,	2	,	7	,
              2	,	3	,	6), ncol = 3, byrow = TRUE)
balancify(d1)
}

\note{Input design should be equireplicate that is, each treatment should have equal replications. Block sizes should be same for each block. For any issue, kindly report to author.}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{incidence matrix}
\keyword{position balanced}
\keyword{nearly position balanced}
\keyword{order balanced}