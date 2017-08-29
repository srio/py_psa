from psa_functions_numeric import *
from psa_functions_symbolic import *

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
ListDistance = [z0]

# Object definition
ListObject = []

# MatTab construction
[MatTabX, MatTabY] = buildMatTab(ListObject, ListDistance)

print('Integral', (2*pi)**2.5 * SigmaXSource * SigmaYSource * SigmaYPSource * SigmaXPSource * SigmaSLambda * SourceI)

IXXP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
ISigma = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[2]

# Symbolic expressions
IXXPSymb = sourceFinaleSymbolic(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYPSymb = sourceFinaleSymbolic(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
print("The symbolic expressions of IXXP is :", IXXPSymb,'and of IYYP :', IYYPSymb)

# def IXint(x, xp, dl):
#     return IXXP(x, xp, dl) * ISigma(dl)
# def IYint(x, xp, dl):
#     return IYYP(x, xp, dl) * ISigma(dl)

# limit calculation
[IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl] = calculateLimits(IXXP, IYYP, ISigma)
print('The integrations boundaries are :', [IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl])

# Sigma calculations
print('Beginning of geometric integration')
SigmaXY = beamGeoSize(IXXP,IYYP,ISigma)
print('Beginning of angular integration')
SigmaXPYP = beamAngularSize(IXXP, IYYP, ISigma)
print('Beginning of flux integration')
# SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(IXXP, IYYP, ISigma, CoefAtten, CoefMonoX, CoefMonoY)
print('SigmaX:%g'%(SigmaXY[0]),' SigmaY:%g'%(SigmaXY[1]))
print('SigmaXP:%g'%(SigmaXPYP[0]), 'SigmaYP:%g'%(SigmaXPYP[1]))
# print('SigmaLambda:%g'%(SigmaLambdaFlux[0]), 'Flux:%g'%(SigmaLambdaFlux[2]))

# Results mathematica simple propagation for z0 = 20000
fluxM = 2.96873 * 10 ** 12
sigmaxyM = [6.08275780e-01, 2.00002544]
sigmaxpypM = [2.99999994e+01, 9.99999981e+01]
sigmalambdaM = 9.99999981e-04

#plotting section
# plotXXP(IXXP, 2, 0.0001, 1000)
# plotYYP(IYYP, 11, 0.0005, 1000)
# plotAnything(IYint, 0.16, 0.0005, 0, 0, 1000)
# plotAnything(IYint, 0, 0.0005, 0.004, 0, 1000)
# plotAnything(IYint, 2, 0, 0.004, 0, 1000)

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=1)

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