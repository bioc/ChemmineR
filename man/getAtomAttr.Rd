\name{getAtomAttr}
\alias{getAtomAttr}
\title{
	getAtomAttr
}
\description{
	On V3000 formatted compounds, returns the value of the given tag on the given atom number.
}
\usage{
getAtomAttr(x,atomId,tag)
}
\arguments{
  \item{x}{
		An \code{SDFset} of \code{ExtSDF} objects. \code{ExtSDF} objects are created with
		\code{read.SDFset} with \code{extendedAttributes=TRUE} when reading V3000 sdf files.
	}
  \item{atomId}{
	  The index of the atom to fetch the tag value from.
	}
	\item{tag}{
		The name of the tag to fetch the value of on the given atom.
	}
}
\value{
	The value of the given tag on the given atom.
}
\author{
	Kevin Horan
}
\examples{
	\dontrun{
		getAtomAttr(v3Sdfs,10,"CHG")
	}
}
