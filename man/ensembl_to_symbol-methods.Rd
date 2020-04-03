% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{ensembl_to_symbol}
\alias{ensembl_to_symbol}
\title{Add transcript symbol column from ensembl id for human and mouse data}
\usage{
ensembl_to_symbol(.data, .ensembl, action = "add")
}
\arguments{
\item{.data}{A `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> |}

\item{.ensembl}{A character string. The column that is represents ensembl gene id}

\item{action}{A character string. Whether to join the new information to the input tbl (add), or just get the non-redundant tbl with the new information (get).}
}
\value{
A `tbl` object including additional columns for transcript symbol
}
\description{
ensembl_to_symbol() takes as imput a `tbl` formatted as | <SAMPLE> | <ENSEMBL_ID> | <COUNT> | <...> | and returns a `tbl` with the additional transcript symbol column
}
\details{
\lifecycle{maturing}

This is useful since different resources use ensembl IDs while others use gene symbol IDs. At the moment this work for human (genes and transcripts) and mouse (genes) data.
}
\examples{


	ensembl_to_symbol(tidybulk::counts_ensembl, ens)


}