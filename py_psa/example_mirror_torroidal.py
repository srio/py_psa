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

# toroidal mirror data
SigmaX0 = 4.000000*10**-7
SigmaY0 = 4.000000*10**-7
DeltaX0 = 1.000000*10**-6
DeltaY0 = 1.000000*10**-5
IncAng0 = 1.000000*10**-1*pi/180
Lambda0 = 1.127140*10**-7
FmX0 = 1.000000*10**4
FmY0 = 8.000000*10**3
SX = 1
SY = 1

# object definition
MirrorToroidal= ['Mirror', IncAng0, SigmaX0, SigmaY0, Lambda0, DeltaX0, DeltaY0, FmX0, FmY0, SX, SY, matrixMirror(IncAng0, SigmaY0, Lambda0, DeltaY0, FmY0, SY), matrixMirror(IncAng0, SigmaX0, Lambda0, DeltaX0, FmX0, SX)]
ListObject = [MirrorToroidal]

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

#plotting section
# def IXint(x, xp, dl):
#     return IXXP(x, xp, dl) * ISigma(dl)
# def IYint(x, xp, dl):
#     return IYYP(x, xp, dl) * ISigma(dl)
# plotXXP(IXXP, IotaX, IotaXp, 500)
# plotYYP(IYYP, IotaY, IotaYp, 500)
# plotAnything(IYYP, 0, IotaYp, IotaYdl, 0, 500)

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
fluxM = 3.04155*10**12
sigmaxyM = [3.18121549e-01, 3.03529776e-01]
sigmaxpypM = [1.00613208e+01, 3.51776350e+01]
sigmalambdaM = 9.99999981e-04

# comparing results
nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, sigmaxpypM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)

