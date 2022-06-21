\name{Sieve-package}
\alias{Sieve-package}
\alias{Sieve}
\docType{package}
\title{
\packageTitle{Sieve}
}
\description{
\packageDescription{Sieve}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{Sieve}
\packageIndices{Sieve}
~~ An overview of how to use the ~~
~~ package, including the most ~~
~~ important functions ~~
}
\author{
\packageAuthor{Sieve}

Maintainer: \packageMaintainer{Sieve}
}
\references{
Tianyu Zhang and Noah Simon (2022) <arXiv:2206.02994>
}
\keyword{sieve estimation; estimation using orthogonal series}
\examples{

xdim <- 5
basisN <- 1000
type <- 'cosine'

#non-linear additive truth. Half of the features are truly associated with the outcome
TrainData <- GenSamples(s.size = 300, xdim = xdim, 
            frho = 'additive', frho.para = xdim/2)

#noise-free testing samples
TestData <- GenSamples(s.size = 1e3, xdim = xdim, noise.para = 0, 
            frho = 'additive', frho.para = xdim/2)

sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
            basisN = basisN, type = type, interaction_order = 2)

sieve.model <- sieve_solver(sieve.model, TrainData$Y, l1 = TRUE)

sieve_model_prediction <- sieve_predict(testX = TestData[,2:(xdim+1)], 
                testY = TestData$Y, sieve.model)

}
