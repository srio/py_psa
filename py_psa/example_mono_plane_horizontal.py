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

#distances
z0 = 10000
z1 = 20000

# Mono data
Alpha0  = 0
Wd0     = 1.000000*10**-6
WidthX0 = 2.000000*10**1
WidthY0 = 2.000000*10**1
RMono0  = 9.600000*10**-1
Rint0   = 1.000000*10**-5
thetaB0 = 1.140288*10**1 * pi /180
bMono0  = np.sin(thetaB0+Alpha0)/np.sin(thetaB0-Alpha0)
Mono = ['MonoPlaneHorizontal', thetaB0, Wd0, RMono0, Rint0, Alpha0, WidthX0, WidthY0,bMono0, matrixMonoPlane(bMono0, thetaB0), np.eye(3)]

ListObject = [Mono]
MatTabX = [matrixFlight(z0), Mono[-2], matrixFlight(z1)]
MatTabY = [matrixFlight(z0), Mono[-1], matrixFlight(z1)]

SigmaXY = beamGeoSize(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)
SigmaXPYP = beamAngularSize(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)
SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, SourceI, CoefAtten, CoefMonoX, CoefMonoY, bMonoX, bMonoY)

print('SigmaX: ',SigmaXY[0],' SigmaY: ',SigmaXY[1])
print('SigmaXP:', SigmaXPYP[0], 'SigmaYP:', SigmaXPYP[1])
print('SigmaLambda:', SigmaLambdaFlux[0], 'Flux:', SigmaLambdaFlux[2])

# # Results mathematica
# Result monochromator vertical
fluxM = 5.95495 * 10 ** 10
sigmaxyM = [6.40040415e-01, 3.00001693e+00]
sigmaxpypM = [2.10726506e+01, 9.99999981e+01]
sigmalambdaM = 1.15334503e-04

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=3)