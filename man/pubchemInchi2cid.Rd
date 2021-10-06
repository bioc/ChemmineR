\name{pubchemInchi2cid}
\alias{pubchemInchi2cid}
\title{
Query pubchem by InChI sttrings and return CIDs
}
\description{
Use PubChem API to get CIDs by InChI sttrings. This function sends 
 one request per InChI. For courtesy, it is not recommended to parellelize this function. 

}
\usage{
pubchemInchi2cid(inchis, verbose = TRUE)
}

\arguments{
  \item{inchis}{
	  Character vector of InChI strings
  }
  \item{verbose}{
	  Logical, show verbose information?
  }
}
\value{
a numeric vector of CIDs with names. Successful requests will have empty 
 names, requests with invalid InChI strings will have name "invalid" and requests 
 with valid InChI but not found in PubChem will have name "not_found"

}
\references{
PubChem PUG REST: https://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST_Tutorial.html

}
\author{
	Le Zhang
}

\examples{
	\dontrun{
	  inchis <- c(
         "InChI=1S/C15H26O/c1-9(2)11-6-5-10(3)15-8-7-14(4,16)13(15)12(11)15/h9-13,16H,5-8H2,1-4H3/t10-,11+,12-,13+,14+,15-/m1/s1", 
         "InChI=1S/C3H8/c1-3-2/h3H2,1-2H3", 
         "InChI=1S/C15H20Br2O2/c1-2-12(17)13-7-3-4-8-14-15(19-13)10-11(18-14)6-5-9-16/h3-4,6,9,11-15H,2,7-8,10H2,1H3/t5-,11-,12+,13+,14-,15-/m1/s1",
         "InChI=abc"
     )
     pubchemInchi2cid(inchis)

	}
}

\keyword{ utilities }
