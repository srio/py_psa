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

# mono bent data
Alpha0 = 0
Wd0 = 1.000000*10**-6
RowlandRadius0 = 5.000000*10**3
WidthX0 = 1.000000*10**1
WidthY0 = 1.000000*10**1
RMono0  = 9.000000*10**-1
Rint0   = 1.300000*10**-5
thetaB0 = 1.035423*10**1*pi/180
Fc0     = 4.493334*10**3
DfS0    = 1.000000*10**4
bMono0  = np.sin(thetaB0+Alpha0)/np.sin(thetaB0-Alpha0)

# object definition
MonoBentVertical = ['MonoBentVertical', Alpha0, thetaB0, Wd0, RowlandRadius0, RMono0, Rint0, DfS0, matrixMonoBent(bMono0, Fc0, thetaB0), np.eye(3, 3)]
ListObject = [MonoBentVertical]

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
[IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl]=calculateLimits(IXXP, IYYP, ISigma)
print([IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl])

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

#plotting section
# def IXint(x, xp, dl):
#     return IXXP(x, xp, dl) * ISigma(dl)
def IYint(x, xp, dl):
    return IYYP(x, xp, dl) * ISigma(dl)
# plotXXP(IXXP, 10*IotaX, 10*IotaXp, 500)
# plotYYP(IYYP, Iota3, Iota4, 500)
# plotAnything(IYint, 0, 6.55e-6, 1.31072e-5, 0, 500)
# plotAnything(IYint, 0, 2.62e-5, 0.00029, 0, 500)

# # Results mathematica
# Result monochromator vertical
fluxM = 1.07708 * 10 ** 8
sigmaxyM = [9.05538410e-01, 4.89927034e-02]
sigmaxpypM = [2.99999994e+01, 2.51906723e+00]
sigmalambdaM = 5.69657118e-06

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)