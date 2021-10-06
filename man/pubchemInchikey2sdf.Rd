\name{pubchemInchikey2sdf}
\alias{pubchemInchikey2sdf}
\title{
Query pubchem by InChIKeys sttrings and get SDF back
}
\description{
Use PubChem API to get CIDs by InChIKeys
}
\usage{
pubchemInchikey2sdf(inchikeys)
}

\arguments{
  \item{inchikeys}{
		Character vector, InChIKey strings.
}
}
\value{
a list of 2 items. the first item "sdf_set" is a `SDFset` object. It 
contains all queried and successful SDF infomation. The second item "sdf_index"
is a named numeric vector. It records whether all input InChIKeys have successful 
returns in the  `SDFset` object. If so, a non-zero value is returned as the index of
where it exists in the `SDFset` object, if not, 0 is returned.
}
\references{
PubChem PUG REST: https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST_Tutorial.html

}
\author{
	Le Zhang
}

\examples{
	\dontrun{
	 ## fetch 2 compounds from PubChem
	 inchikeys <- c(
	 "ZFUYDSOHVJVQNB-FZERPYLPSA-N", 
	 "KONGRWVLXLWGDV-BYGOPZEFSA-N", 
	 "AANKDJLVHZQCFG-WLIQWNBFSA-N",
	 "SNFRINMTRPQQLE-JQWAAABSSA-N"
	 )
	 pubchemInchikey2sdf(inchikeys)

	}
}

\keyword{ utilities }
