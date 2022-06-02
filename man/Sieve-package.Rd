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
~~ Literature or other references for background information ~~
}
\keyword{sieve estimation}
\examples{

xdim <- 5
basisN <- 1000
type <- 'cosine'

TrainData <- Sieve:::GenSamples(s.size = 300, xdim = xdim, 
            frho = 'additive', frho.para = xdim/2)

TestData <- Sieve:::GenSamples(s.size = 1e3, xdim = xdim, noise.para = 0, 
            frho = 'additive', frho.para = xdim/2)

sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
            basisN = basisN, type = type, interaction_order = 2)

sieve.model <- sieve_solver(sieve.model, TrainData$Y, l1 = TRUE)

sieve_model_prediction <- sieve_predict(testX = TestData[,2:(xdim+1)], 
                testY = TestData$Y, sieve.model)

}
