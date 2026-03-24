\name{simulate_lineage}
\alias{simulate_lineage}
\title{Simulate lineage under dual barcode models}
\usage{
simulate_lineage(Q.list, Q.list1, N, t_max = 500, cell.num = 1,
  gRNA.num = 200, death_rate_per_time = 1/240, p_quiesce_enter = 0.2,
  quiesce_mean = 24, quiesce_sd = 6, p_senesce = 0.01, div_mean = 36,
  div_sd = 8, verbose = TRUE)
}
\arguments{
  \item{Q.list, Q.list1}{Lists of transition matrices.}
  \item{N}{Number of states per gRNA.}
  \item{t_max}{Maximum time.}
  \item{cell.num}{Initial cells.}
  \item{gRNA.num}{Number of gRNAs.}
  \item{death_rate_per_time}{Death hazard.}
  \item{p_quiesce_enter}{Quiescence probability at division.}
  \item{quiesce_mean, quiesce_sd}{Quiescence duration parameters.}
  \item{p_senesce}{Senescence probability.}
  \item{div_mean, div_sd}{Inter-division time parameters.}
  \item{verbose}{Print progress.}
}
\value{
List with \code{edges}, \code{terminal}, and \code{history}.
}
\description{
Main branching-process simulation for hg and sg models.
}
