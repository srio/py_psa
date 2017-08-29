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

# mosaic monochromator data
eta0 = 5.200000*10**-3
WidthX0 = 1.000000*10**1
WidthY0 = 1.100000*10**-2
RMono0 = 1.000000*10**-2
Rint0 = 1.100000*10**-2
thetaB0 = 9.136369*10**-2*pi/180

# object definition
MonoMosaicHorizontal = ['MonoMosaicHorizontal', thetaB0, eta0, Rint0, np.eye(3,3), matrixMonoMosaic(thetaB0),]
ListObject = [MonoMosaicHorizontal]

# MatTab construction
[MatTabX, MatTabY] = buildMatTab(ListObject, ListDistance)

#function definition
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

# Sigma calculations
print('Beginning of geometric integration')
SigmaXY = beamGeoSize(IXXP,IYYP,ISigma,)
print('Beginning of angular integration')
SigmaXPYP = beamAngularSize(IXXP, IYYP, ISigma)
print('Beginning of flux integration')
# SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(IXXP, IYYP, ISigma, SigmaXPSource, SigmaYPSource, SigmaXSource, SigmaYSource)
print('SigmaX:%g'%(SigmaXY[0]),' SigmaY:%g'%(SigmaXY[1]))
print('SigmaXP:%g'%(SigmaXPYP[0]), 'SigmaYP:%g'%(SigmaXPYP[1]))
# print('SigmaLambda:%g'%(SigmaLambdaFlux[0]), 'Flux:%g'%(SigmaLambdaFlux[2]))

#plotting section
# def IXint(x, xp, dl):
#     return IXXP(x, xp, dl) * ISigma(dl)
# def IYint(x, xp, dl):
#     return IYYP(x, xp, dl) * ISigma(dl)
# plotXXP(IXXP, Iota1, Iota2)
# plotYYP(IYYP, Iota3, Iota4)
# plotAnything(IYint, 0, 0, 0, 0, 1000)

# # Results mathematica
# Result monochromator mosaic horizontal
fluxM = -2.56758 * 10 ** 15
sigmaxyM = [3.22591529e-01, 6.17151193e-01]
sigmaxpypM = [3.01685473e+01, 9.99999981e+01]
sigmalambdaM = 9.99999981e-04

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)