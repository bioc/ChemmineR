\name{read.SDFset}
\alias{read.SDFset}
\title{
SD file to \code{SDFset} 
}
\description{
Imports one or many molecules from an SD/MOL file and stores it in an \code{SDFset} container. 
Supports both the V2000 and V3000 formats.
}
\usage{
read.SDFset(sdfstr = sdfstr,skipErrors=FALSE, ...)
}
\arguments{
  \item{sdfstr}{
		path/name to an SD file; alternatively an \code{SDFstr} object can be provided
  }
  \item{skipErrors}{
	  If true, molecules which fail to parse will be removed from the output. Otherwise
	  and error will be thrown and processing of the input will stop. 
  }
  \item{\dots}{
		option to pass on additional arguments. Possible arguments are: 

		datablock: true or false, whether to include the data fields or not.
		Defaults to TRUE.

		tail2vec: true or false, whether to return data feilds as a vector or not. 
		Defaults to TRUE.

		extendedAttributes: true or false, whether to parse the extended attributes available on the
		V3000 format. Defaults to FALSE. When set to TRUE, the resulting objects will be of type
		ExtSDF, which is a sub-class of SDF. However, some functions, such as plot, may not work with
		this type right now. 

	}
}
\details{
...
}
\value{
\item{SDFset}{ for details see ?"SDFset-class"}
}
\references{
SDF format defintion: http://www.symyx.com/downloads/public/ctfile/ctfile.jsp
}
\author{
Thomas Girke
}

\seealso{
Functions: \code{read.SDFstr}
}
\examples{
## Write instance of SDFset class to SD file
data(sdfsample); sdfset <- sdfsample
# write.SDF(sdfset[1:4], file="sub.sdf")

## Import SD file 
# read.SDFset("sub.sdf")

## Pass on SDFstr object
sdfstr <- as(sdfset, "SDFstr")
read.SDFset(sdfstr) 
}
\keyword{ utilities }
