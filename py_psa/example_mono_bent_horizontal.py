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
RowlandRadius0 = 5.000000*10**4
WidthX0 = 1.000000*10**1
WidthY0 = 1.000000*10**1
RMono0  = 9.000000*10**-1
Rint0   = 1.300000*10**-5
thetaB0 = 1.035423*10**1*pi/180
Fc0     = 4.493334*10**3
DfS0    = 1.000000*10**4
bMono0  = np.sin(thetaB0+Alpha0)/np.sin(thetaB0-Alpha0)

# object definition
MonoBentHorizontal = ['MonoBentHorizontal', Alpha0, thetaB0, Wd0, RowlandRadius0, RMono0, Rint0, DfS0, np.eye(3, 3), matrixMonoBent(bMono0, Fc0, thetaB0)]
ListObject = [MonoBentHorizontal]

# MatTab construction
[MatTabX, MatTabY] = buildMatTab(ListObject, ListDistance)
# MatTabX = [matrixFlight(z0), matrixMonoBent(bMono0, Fc0, thetaB0), matrixFlight(z1)]
# MatTabY = [matrixFlight(z0), np.eye(3, 3), matrixFlight(z1)]


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

#plotting section
def IXint(x, xp, dl):
    return IXXP(x, xp, dl) * ISigma(dl)
# def IYint(x, xp, dl):
#     return IYYP(x, xp, dl) * ISigma(dl)
# plotXXP(IXXP, 0.42, 5.24e-5, 500)
# plotYYP(IYYP, 0.215, 0.0002, 500)
# plotAnything(IXint, 0, 2.62144e-05, 0.0002097152, 10**-2, 500)
# plotAnything(IXint, 0, 2.62144e-05, 0.0002097152, 0, 500)
plotAnything(IXint, 0.429496, 0, 0.0002097152, 10**-6, 500)
plotAnything(IXint, 0.429496 , 0, 0.0002097, 0, 500)


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
fluxM = 3.23434*10**9
sigmaxyM = [4.34228510e-01 , 3.00001693e+00]
sigmaxpypM = [3.54741332e+01, 9.99999981e+01]
sigmalambdaM = 4.14803671e-05

# comparing results
nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)