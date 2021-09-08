\name{pbbd}
\alias{pbbd}
\title{Position balanced and nearly position balanced block design}
\description{This function generates a position balanced or nearly position balanced block design with given parameters. User needs to specify number of treatments (v), number of blocks (b) and block size (k)}
\usage{pbbd(v, b, k)}
\arguments{
	\item{v}{Number of treatments}	
	\item{b}{Number of blocks}	
	\item{k}{Block size}	
	}
\value{
	\item{parameters}{Parameters v,b,r,k. Here r is number of replications of each treatment}
	\item{efficiencies}{A- and D-efficiency of the design generated}
	\item{design}{Position balanced block design}
	\item{P}{Treatment verus position incidence matrix} 
  }
\examples{
pbbd(7,7,3)

pbbd(9,12,3)
}

\note{This function works for generating a position balanced block design for upto 30 treatments and block size 10. For getting design with larger number of treatments and/or block size, it is better to use balancify() function with a design supplied by user to make the design position balanced.}
\author{B N Mandal <mandal.stat@gmail.com>}
\keyword{incomplete block design}
\keyword{incidence matrix}
\keyword{position balanced}
\keyword{nearly position balanced}
\keyword{order balanced}