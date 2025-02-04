# Example 1.3: Classifying categorical observations into two states

PrA <- 0.7 # Probability of items being produced by the reliable firm
PrAC <- 1 - PrA # Probability of items being produced by the less reliable firm

PrFA <- 0.01 # Failure rate of the reliable firm
PrFAC <- 0.05 # Failure rate of the less reliable firm

N <- 100
y <- 0:6

PosteriorTheta1Unormalized <- PrFA ^ y * (1 - PrFA) ^ (N - y) * PrA
PosteriorTheta0Unormalized <- PrFAC ^ y * (1 - PrFAC) ^ (N - y) * PrAC
NormalizingConstant <- PosteriorTheta0Unormalized + PosteriorTheta1Unormalized

round(PosteriorTheta1Unormalized / NormalizingConstant, 3)
round(PosteriorTheta0Unormalized / NormalizingConstant, 3)
