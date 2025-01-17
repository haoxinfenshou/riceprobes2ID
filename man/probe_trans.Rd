% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/probe_trans.R
\name{probe_trans}
\alias{probe_trans}
\title{probe_trans}
\usage{
probe_trans(GPL, expr_matrix, merge_by = "mean")
}
\arguments{
\item{GPL}{needs a GPL ID. 2021-12-27, 3 rice microarray platforms are available: GPL6864, GPL8852 and GPL2025}

\item{expr_matrix}{needs a dataframe, rownames should be microarray ID and colnames should be sample GSM ID}

\item{merge_by}{needs a method to merge multi probes and now only mean value is available}
}
\value{
a dataframe
}
\description{
Convert rice microarray ID to RAP-ID and merge multi probes for one genes
}
\details{
Convert rice microarray ID to RAP-ID and merge multi probes for one genes
}
