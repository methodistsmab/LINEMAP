\name{prepare_Q_lists}
\alias{prepare_Q_lists}
\title{Prepare per-gRNA transition matrix lists}
\usage{
prepare_Q_lists(Q, gRNA.num = 200, mutation.rate = seq(0.05, 0.25, 0.05),
  recycle_blocks = 40)
}
\arguments{
  \item{Q}{Square transition matrix.}
  \item{gRNA.num}{Number of gRNA positions.}
  \item{mutation.rate}{Mutation rates for first-row patterns.}
  \item{recycle_blocks}{Block count; \code{gRNA.num = length(mutation.rate) * recycle_blocks}.}
}
\value{
List with \code{Q}, \code{Q.list}, \code{Q.list1}, and \code{N}.
}
\description{
Builds per-gRNA transition lists from a base matrix \code{Q}.
}
