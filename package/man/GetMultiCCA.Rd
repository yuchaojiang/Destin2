% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Destin2.R
\name{GetMultiCCA}
\alias{GetMultiCCA}
\title{MultiCCA of Seurat}
\usage{
GetMultiCCA(
  obj,
  reduction.list,
  dims.list,
  assay = NULL,
  reduction.name = NULL,
  reduction.key = NULL,
  n.cca = NULL
)
}
\arguments{
\item{obj}{Seurat object}

\item{reduction.list}{A list of of dimension reduction objects within obj}

\item{dims.list}{A list of which dimensions to use for each combined dimension reduction object}

\item{assay}{A character string specifying the assay used to calculate this dimensional reduction}

\item{reduction.name}{A character string providing the name of the MCCA dimension reduction}

\item{reduction.key}{A character string that acts as a prefix for the MCCA dimension reduction}

\item{n.consensus.pc}{Number of principal components (PC's) to store}
}
\value{
A Seurat object containing MCCA joint dimension reduction in slot "multicca.pca"
}
\description{
This function takes a Seurat object and conducts multiple canonical correlation analysis (MCCA) as joint dimension reduction
}
