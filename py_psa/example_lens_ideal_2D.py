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

# lens data
F0 = 1.0/(1.0/z0 + 1.0/z1) # 7.352941*10**3

# object definition
Lens = ['LensIdeal2D', matrixCompLensPara(F0), matrixCompLensPara(F0)]
ListObject = [Lens]

# MatTab construction
[MatTabX, MatTabY] = buildMatTab(ListObject, ListDistance)

#function definition
IXXP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYP = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
ISigma = sourceFinale(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[2]

# Symbolic expressions
IXXPSymb = sourceFinaleSymbolic(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[0]
IYYPSymb = sourceFinaleSymbolic(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY)[1]
print("The symbolic expressions of IXXP is :", IXXPSymb,'\n and of IYYP :', IYYPSymb)


# limit calculation
[IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl] = calculateLimits(IXXP, IYYP, ISigma)
print('The integrations boundaries are :', [IotaX, IotaXp, IotaY, IotaYp, IotaXdl, IotaYdl])


#plotting section
# def IXint(x, xp, dl):
#     return IXXP(x, xp, dl) * ISigma(dl)
# def IYint(x, xp, dl):
    # return IYYP(x, xp, dl) * ISigma(dl)

# plotXXP(IXXP, IotaX, IotaXp, 500)
# plotYYP(IYYP, IotaY, IotaYp, 500)
plotAB(IXXP, IotaX, IotaXp, 500, title="XXP", xtitle="X", ytitle="XP")
plotAB(IYYP, IotaY, IotaYp, 500, title="YYP", xtitle="Y", ytitle="YP")
# plotAnything(IYint, 0, IotaYp, IotaYdl, 0, 500)


# Sigma calculations
print('Beginning of geometric integration')
SigmaXY = beamGeoSize(IXXP,IYYP,ISigma)
print('Beginning of angular integration')
SigmaXPYP = beamAngularSize(IXXP, IYYP, ISigma)
print('SigmaX:%g'%(SigmaXY[0]),' SigmaY:%g'%(SigmaXY[1]))
print('SigmaXP:%g'%(SigmaXPYP[0]), 'SigmaYP:%g'%(SigmaXPYP[1]))

# print('Beginning of flux integration')
# SigmaLambdaFlux = sigma1_MaxFluxL_FluxPhi(IXXP, IYYP, ISigma, CoefAtten, CoefMonoX, CoefMonoY)
# print('SigmaLambda:%g'%(SigmaLambdaFlux[0]), 'Flux:%g'%(SigmaLambdaFlux[2]))



# Results theory
fluxM = SourceI
sigmaxyM   = np.array([SigmaXSource,SigmaYSource])   * z1 / z0 # [1.90707991e-01 , 9.90504563e-02]
sigmaxpypM = np.array([SigmaXPSource,SigmaYPSource]) * z0 / z1 # [1.45663785e+01, 1.25251185e+01]
sigmalambdaM = 9.99999981e-04


print("Magnification (from sizes) Horizontal: %f (theory: %f)"%(SigmaXY[0]/SigmaXSource,z1/z0))
print("Magnification (from sizes) Vertical:   %f (theory: %f)"%(SigmaXY[1]/SigmaYSource,z1/z0))
print("Magnification (from divergences) Horizontal: %g (theory: %f)"%(1e6*SigmaXPSource/SigmaXPYP[0],z1/z0))
print("Magnification (from divergences) Vertical:   %g (theory: %f)"%(1e6*SigmaYPSource/SigmaXPYP[1],z1/z0))

# Testing section
nut.assert_array_almost_equal(SigmaXY, sigmaxyM, decimal=2)
nut.assert_array_almost_equal(SigmaXPYP, 1e6*sigmaxpypM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[0], sigmalambdaM, decimal=2)
# nut.assert_almost_equal(SigmaLambdaFlux[2], fluxM , decimal=2)
