% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primer3Functions.R
\name{get.primer}
\alias{get.primer}
\title{Get Primer Pairs Given A Gtf And A Transcript ID}
\usage{
get.primer(seq_id, transcriptome, gtf, type = "overlap",
  product_size_range = c(90, 120), Tm = c(60, 62.5, 65),
  primer_size_range = c(15, 18, 21))
}
\arguments{
\item{seq_id}{Sequence to look for.}

\item{transcriptome}{Transcriptome fasta file.}

\item{gtf}{Gtf to look in for sequence.}

\item{type}{One of either "overlap" or "target"}

\item{product_size_range}{default: '90-120'}

\item{Tm}{melting temprature parameters default: c(60, 62.5, 65)}

\item{primer_size_range}{default: c(15, 18, 21)}
}
\description{
`get.primer` returns one list of primer pairs per exon-exon junction. One of
 the primers in each pair will be spanning the junction.
}
