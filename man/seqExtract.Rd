% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seqExtract.R
\name{seqExtract}
\alias{seqExtract}
\title{Sequence extraction}
\usage{
seqExtract(positions, genome, seq_range = c(-5000, 5000))
}
\arguments{
\item{positions}{genpomic coordinates}

\item{genome}{genome sequence (from BSgenome data)}

\item{seq_range}{sequence extraction target range from position}
}
\value{
a DNAStringSet format sequences
}
\description{
Sequence extraction from the given range from position.
}
\keyword{genome}
\keyword{sequence}
