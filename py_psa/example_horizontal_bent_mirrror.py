from psa_functions import *

SourceI = 1e25
SigmaXSource = 1e-1
SigmaYSource = 1e-2
SigmaYPSource = 10**-4
SigmaXPSource = 3e-5
GammaSource = 0
SigmaSLambda = 1e-3
bMonoX = 1
bMonoY = 1
CoefMonoX = 1
CoefMonoY = 1
CoefAtten = 1


z0 = 10000
z1 = 20000
# # Mirror data
Sigma0 = 0
Delta0 = 1.000000*10**-7
IncAng0 = 0
Lambda0 = 1.127140*10**-7
Fm0 = 8.333333*10**3
S = 1
Mirror = [np.eye(3), matrixMirror(IncAng0, Sigma0, Lambda0, Delta0, Fm0, S)]

ListObject = [Mirror]
MatTabX = [matrixFlight(z0), Mirror[-1], matrixFlight(z1)]
MatTabY = [matrixFlight(z0), Mirror[-2], matrixFlight(z1)]

SigmaXY = beamGeoSize(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)
SigmaXPYP = beamAngularSize(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)
SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, CoefAtten, CoefMonoX, CoefMonoY, bMonoX, bMonoY)

print('SigmaX: ',SigmaXY[0],' SigmaY: ',SigmaXY[1])
print('SigmaXP:', SigmaXPYP[0], 'SigmaYP:', SigmaXPYP[1])
print('SigmaLambda:', SigmaLambdaFlux[0], 'Flux:', SigmaLambdaFlux[2])

#data horizontal bent mirror mathematica
fluxM = 2.96873 * 10 ** 12
sigmaxyM = [2.28035259e-01, 3.00001693]
sigmaxpypM = [1.34163968e+01, 9.99999981e+01]
sigmalambdaM = 9.99999981e-04

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)