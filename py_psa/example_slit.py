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
z1 = 12000
AperturSlitX0 = 20
AperturSlitY0 = 30
Slit = ['Slit',AperturSlitX0, AperturSlitY0, np.eye(3), 0]
ListObject = [Slit]
MatTabX = [matrixFlight(z0), Slit[-2], matrixFlight(z1)]
MatTabY = [matrixFlight(z0), Slit[-2], matrixFlight(z1)]

SigmaXY = beamGeoSize(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)
SigmaXPYP = beamAngularSize(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)
SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, CoefAtten, CoefMonoX, CoefMonoY, bMonoX, bMonoY)

print('SigmaX: ',SigmaXY[0],' SigmaY: ',SigmaXY[1])
print('SigmaXP:', SigmaXPYP[0], 'SigmaYP:', SigmaXPYP[1])
print('SigmaLambda:', SigmaLambdaFlux[0], 'Flux:', SigmaLambdaFlux[2])

# # Results mathematica single slit
fluxM = 4.24833 * 10 ** 12
sigmaxyM = [6.66564735e-01, 2.18549439]
sigmaxpypM = [2.99596142e+01, 9.93400750e+01]
sigmalambdaM = 9.99999981e-04

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=3)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=3)
nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=3)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=3)