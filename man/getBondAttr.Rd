\name{getBondAttr}
\alias{getBondAttr}
\title{
	getBondAttr
}
\description{
	On V3000 formatted compounds, returns the value of the given tag on the given bond number.
}
\usage{
getBondAttr(x,bondId,tag)
}
\arguments{
  \item{x}{
		An \code{SDFset} of \code{ExtSDF} objects. \code{ExtSDF} objects are created with
		\code{read.SDFset} with \code{extendedAttributes=TRUE} when reading V3000 sdf files.
	}
  \item{bondId}{
	  The index of the bond to fetch the tag value from.
	}
	\item{tag}{
		The name of the tag to fetch the value of on the given bond.
	}
}
\value{
	The value of the given tag on the given bond.
}
\author{
	Kevin Horan
}
\examples{
	\dontrun{
		getBondAttr(v3Sdfs,10,"CFG")
	}
}
