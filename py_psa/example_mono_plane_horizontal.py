from psa_functions_numeric import *

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

IXXP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
ISigma = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[2]

# def IXint(x, xp, dl):
#     return IXXP(x, xp, dl) * ISigma(dl)
# def IYint(x, xp, dl):
#     return IYYP(x, xp, dl) * ISigma(dl)

[IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl] = calculateLimits(IXXP, IYYP, ISigma)
print('The integrations boundaries are :', [IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl])

# plotting section
# plotXXP(IXXP, 0.1, 5*10**-6, 500)
# plotYYP(IYYP, 10, 0.0005, 1000)
# plotAnything(IXint, 0.1, 5*10**-6, 0, 0, 500)
# plotAnything(IXint, 0, 10**-5, 10**-4, 0, 500)
# plotAnything(IXint, 0.1, 0, 5*10**-5, 0, 500)


SigmaXY = beamGeoSize(IXXP,IYYP,ISigma, SigmaXPSource, SigmaYPSource, SigmaSLambda)
print('Beginning of angular integration')
SigmaXPYP = beamAngularSize(IXXP, IYYP, ISigma, SigmaXSource, SigmaYSource, SigmaSLambda)
print('Beginning of flux integration')
SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(IXXP, IYYP, ISigma, SigmaXPSource, SigmaYPSource, SigmaXSource, SigmaYSource, SigmaSLambda, CoefAtten, CoefMonoX, CoefMonoY)
print('SigmaX: ',SigmaXY[0],' SigmaY: ',SigmaXY[1])
print('SigmaXP:', SigmaXPYP[0], 'SigmaYP:', SigmaXPYP[1])
print('SigmaLambda:', SigmaLambdaFlux[0], 'Flux:%g'%(SigmaLambdaFlux[2]))


# # Results mathematica
# Result monochromator vertical
fluxM = 5.95495 * 10 ** 10
sigmaxyM = [6.40040415e-01, 3.00001693e+00]
sigmaxpypM = [2.10726506e+01, 9.99999981e+01]
sigmalambdaM = 1.15334503e-04

nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)