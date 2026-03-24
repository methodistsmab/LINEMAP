\name{sample_subtrees}
\alias{sample_subtrees}
\title{Sample matched subtrees from two phylogenies}
\usage{
sample_subtrees(phy, phy1, n_keep = 30, seed = NULL)
}
\arguments{
  \item{phy, phy1}{Matched \code{phylo} objects.}
  \item{n_keep}{Number of tips to keep.}
  \item{seed}{Optional random seed.}
}
\value{
List with \code{keep_tips}, \code{phy_sub}, \code{phy_sub1}.
}
\description{
Randomly subsamples tips and returns subtrees from both trees.
}
