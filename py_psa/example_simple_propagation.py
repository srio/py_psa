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

z0 = 20000
MatTabX = [matrixFlight(z0)]
MatTabY = [matrixFlight(z0)]
ListObject = []

print('Integral', (2*pi)**2.5 * SigmaXSource * SigmaYSource * SigmaYPSource * SigmaXPSource * SigmaSLambda * SourceI)

IXXP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
ISigma = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[2]

def IYint(y, yp, dl):
    return IYYP(y, yp, dl) * ISigma(dl)

print(calculateLimits(IXXP, IYYP, ISigma))
plotXXP(IXXP, 2.68435456, 0.00016384)
#[2.68435456, 0.00016384, 0.16777216, 1.024e-05, 8.192e-05, 0.01048576]
plotYYP(IYYP, 0.01, 0.0001)
# plotAnything(IYint, 0.04, 0.032, 0.032, 0, 12)
# plotAnything(IYint, 0.04, 0.032, 0.032, 0, 23)
# plotAnything(IYint, 0.04, 0.032, 0.032, 0, 13)


SigmaXY = beamGeoSize(IXXP,IYYP,ISigma, SigmaXPSource, SigmaYPSource, SigmaSLambda)
SigmaXPYP = beamAngularSize(IXXP, IYYP, ISigma, SigmaXSource, SigmaYSource, SigmaSLambda)
SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(IXXP, IYYP, ISigma, SigmaXPSource, SigmaYPSource, SigmaXSource, SigmaYSource, SigmaSLambda, CoefAtten, CoefMonoX, CoefMonoY)
print('SigmaX: ',SigmaXY[0],' SigmaY: ',SigmaXY[1])
print('SigmaXP:', SigmaXPYP[0], 'SigmaYP:', SigmaXPYP[1])
print('SigmaLambda:', SigmaLambdaFlux[0], 'Flux:', SigmaLambdaFlux[2])

# Results mathematica simple propagation
fluxM = 2.96873 * 10 ** 12
sigmaxyM = [6.08275780e-01, 2.00002544]
sigmaxpypM = [2.99999994e+01, 9.99999981e+01]
sigmalambdaM = 9.99999981e-04

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)

# x, y, xp, yp and dl
# SigmaSource = [0.1, 0.01, 3e-5, 1e-4, 1e-3]
# [2.68435456, 0.16777216, 0.00065536, 0.00262144, 0.1048576]
# [26.8435456, 16.777216, 21.845333333333333, 26.214399999999998, 104.85759999999]
# fluxM = 2.96873 * 10 ** 12
# sigmaxyM = [0.1, 0.01]
# sigmaxpypM = [2.99999994e+01, 9.99999981e+01]
# sigmalambdaM = 9.99999981e-04
#
#
# [2.68435456, 0.16777216, 0.00016384, 2.048e-05, 0.1048576]
# [26.8435456, 16.777216, 5.461333333333333, 0.20479999999999998, 104.85759999999999]
# fluxM = 1.9816353937496942 * 10**11
# sigmaxyM = [0.316227762923, 1.00004999875]
# sigmaxpypM = [2.99999994e+01, 9.99999981e+01]
# sigmalambdaM = 9.99999981e-04