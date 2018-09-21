\name{MotScr}
\alias{}
\title{MotScr}
\usage{
MotScr()
}
\description{
Motif overrepresentation analysis for differetially methylated regions.
}
\examples{
MotScr(infile="sel_processed_Mval.txt", motifDBList = motifDBList, cutoff = 2, p.cutoff = 0.001, outname="screening_result", ControlColnum=c(1,2), TreatmentColnum=c(3,4), MethylDemethyl="Demethyl", version = "850")
}
