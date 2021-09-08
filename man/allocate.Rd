\name{allocate}
\alias{allocate}
\title{Allocate treatments in a column}
\description{This function alloctes treatments in a column of a block design}
\usage{allocate(j,v,b,k,mvec,x1,x2)}
\arguments{
	\item{j}{An integer less than or equal to block size}
	\item{v}{Number of treatments}
	\item{b}{Number of blocks}
	\item{k}{Block size}
	\item{mvec}{Vector of desired frequencies of treatments in jth column}
	\item{x1}{A b x k matrix}
	\item{x2}{A b x k matrix}
	}
\value{
\item{x1}{A b x k matrix}
\item{x2}{A b x k matrix}
}
\author{Baidya Nath Mandal <mandal.stat@gmail.com>}
\keyword{internal}