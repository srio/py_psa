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

# distances
z0 = 10000
z1 = 20000
ListDistance = [z0, z1]

# Mono data
Alpha0  = 0
Wd0     = 1.000000*10**-6
WidthX0 = 2.000000*10**1
WidthY0 = 2.000000*10**1
RMono0  = 9.600000*10**-1
Rint0   = 1.000000*10**-5
thetaB0 = 1.140288*10**1 * pi /180
bMono0  = np.sin(thetaB0+Alpha0)/np.sin(thetaB0-Alpha0)

# object definition
Mono = ['MonoPlaneHorizontal', thetaB0, Wd0, RMono0, Rint0, Alpha0, WidthX0, WidthY0,bMono0, np.eye(3), matrixMonoPlane(bMono0, thetaB0)]
ListObject = [Mono]

# MatTab construction
[MatTabX, MatTabY] = buildMatTab(ListObject, ListDistance)

# function definition
IXXP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
ISigma = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[2]

# Symbolic expressions
IXXPSymb = sourceFinaleSymbolic(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYPSymb = sourceFinaleSymbolic(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
print("The symbolic expressions of IXXP is :", IXXPSymb,'and of IYYP :', IYYPSymb)

# limit calculation
[IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl] = calculateLimits(IXXP, IYYP, ISigma)
print('The integrations boundaries are :', [IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl])

# plotting section
def IXint(x, xp, dl):
    return IXXP(x, xp, dl) * ISigma(dl)
# def IYint(x, xp, dl):
#     return IYYP(x, xp, dl) * ISigma(dl)
# plotXXP(IXXP, 0.1, 5*10**-6, 500)
# plotYYP(IYYP, 10, 0.0005, 1000)
# plotXXP(IXXP, IotaX, IotaXp, 500)
# plotYYP(IYYP, IotaY, IotaYp, 500)
# plotAnything(IXint, 0.1, 5*10**-6, 0, 0, 500)
# plotAnything(IXint, 0, 2.62144e-05, 5.24288e-05, 0.01, 500)
# plotAnything(IXint, 0, 2.62144e-05, 0.0001048576, 0, 500)

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

# # Results mathematica
# Result monochromator horizontal
fluxM = 5.95495 * 10 ** 10
sigmaxyM = [6.40040415e-01, 3.00001693e+00]
sigmaxpypM = [2.10726506e+01, 9.99999981e+01]
sigmalambdaM = 1.15334503e-04

# Testing section
nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)