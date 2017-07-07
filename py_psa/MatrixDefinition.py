import numpy as np
from numpy import exp
from numpy import pi
from numpy import sqrt
from numpy import log
import numpy.testing as nut
import sympy as sp
from sympy.integrals.quadrature import gauss_hermite


x, xp, y, yp, dl = sp.symbols('x xp y yp dl')

#Definition of mathematical objects


    #Definition of the transformation matrices

def MatFlight(L):                                                   # For a flight path of length
    return np.array([[1,-L,0],[0,1,0],[0,0,1]])

def MatFlatMono(b, ThetaB):                                         # For a perfect flat crystal monochromator
    return np.array([[b,0,0],[0,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatCurvMono(b, Fc, ThetaB):                                     # For a perfect curved crystal monochromator (meridionally and sagitally focusing)
    return np.array([[b,0,0],[1/Fc,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatMosMono(ThetaB):                                             # For a mosaic monochromator
    return np.array([[1,0,0],[0,-1,2*np.tan(ThetaB)],[0,0,1]])

def MatFlatMirr(IncAng, Sigma, Lambda, Delta):                      # For the plane mirror
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[Delta,1,0],[0,0,1]])

def MatMirr(IncAng, Sigma, Lambda, Delta, Fm, S):                      # For the bent and toroidal mirrors
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[(1+S*Fm*Delta)/Fm,1,0],[0,0,1]])

def MatParaLens(F):                                                 # For the compound refractive lenses with parabolic holes
    return np.array([[1,0,0],[1/F,1,0],[0,0,1]])

def MatCircLens(Coef, F):                                           # For the compound refractive lenses with circular holes
    return np.array([[1,0,0],[Coef/F,1,0],[0,0,1]])

def MatMulti(Fmu):                                                  # For the multilayer
    return np.array([[1,0,0],[1/Fmu,1,0],[0,0,1]])

    #Definition of the beam source

def SourceX(x, xp, SigmaSourceX=1e-6, SigmaSourceXp=2e-6):
    return sp.exp(-( (x/SigmaSourceX)**2 + (xp/SigmaSourceXp)**2) / 2 )
def SourceY(y, yp, SigmaSourceY=1e-6, SigmaSourceYp=2e-6, GammaSource=5e-5):
    return sp.exp( -( (y/SigmaSourceY)**2 + ((yp-GammaSource*y)/SigmaSourceYp)**2)/2 )
def SourceLambda(dl, SigmaSLambda=1e-6):
    return sp.exp(-(dl)**2/2/SigmaSLambda**2)

    #Definition of the acceptance windows of the optical parts

def cotan(x):
    return 1/np.tan(x)

def AccSlit(y, Aperture, calctype):                                                           #Slit acceptance
    if calctype==0:
        return sqrt(6/pi) / sqrt(6*log(2)/pi) * sp.exp( -(y/Aperture)**2/2*12)
    if calctype==1:
        return 1/sqrt(6*log(2)/pi)*sp.exp(-y**2/(2*Aperture**2/2/pi))
    else:
        return 'the value for calctype is 0 or 1 you pipsqueak !'

def AccPinh(y, Diameter):                                                            #Pinhole acceptance
    return sqrt(8/pi) * exp ( -(y/Diameter)**2/2*16 )

def AccMonoPlanAng(RMono, RInt, Wd, yp, DeltaLambda, ThetaB):                #Plane monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * exp( -(yp-DeltaLambda*np.tan(ThetaB))**2 / (2*Wd**2/12))

def AccMonoPlanWav(DeltaLambda, SigmaYp, ThetaB, Wd):                        #Plane monochromator wave acceptance
    return sqrt(6/pi) * exp( - (DeltaLambda)**2 / (2*(SigmaYp**2+Wd**2/12)*cotan(ThetaB)**2) )

def AccMonoMosaAng(RInt, eta, yp, DeltaLambda, ThetaB):                       #Mosaic monochromator angular acceptance
    return RInt*sqrt(6/pi)/eta * exp(- (yp-DeltaLambda*np.tan(ThetaB))**2 / eta**2 /2 )

def AccMonoMosaWav(DeltaLambda, SigmaYp, eta, ThetaB):                        #Mosaic monochromator wave acceptance
    return sqrt(6/pi) * exp( - (DeltaLambda)**2 / 2 /((SigmaYp**2+eta**2)*cotan(ThetaB)**2))

def AccMonoCurvAng(RMono, RInt, Wd,yp, y, r, ThetaB, Alpha, DeltaLambda ):    #Curved monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * exp( - (yp-y/r/np.sin(ThetaB+Alpha)-DeltaLambda*np.tan(ThetaB))**2/(2*Wd**2/12))

def AccMonoCurvWav(DeltaLambda, ThetaB, SigmaSource, DistanceFromSource, Wd):  #Curved monochromator wave acceptance
    return sqrt(6/pi) * exp( -(DeltaLambda)**2 / (2*cotan(ThetaB)**2*((SigmaSource/DistanceFromSource)**2+Wd**2/12))   )

def AccLensPara(x, Radius, Distance, Sigma):                                          #Compound refractive lense with parabolic holes acceptance
    return exp( -(x**2+Radius*Distance)/2/Sigma**2)

def AccLensCirc(x, Sigma, FWHM, Radius, Distance):                                    #Compound refractive lense with circular holes acceptance
    return exp( -x**2/2/Sigma**2 -x**2*FWHM**2/8/Radius**2/Sigma**2 -x**2*FWHM**4/16/Sigma**2/Radius**4 -Radius*Distance/2/Sigma**2  )

def AccMultAng(Rml, yp, y, Rowland, ThetaML, DeltaLambda, Ws):                #Multilayer angle acceptance
    return Rml*8/3*sqrt(log(2)/pi) * exp( -(-yp-y/Rowland/np.sin(ThetaML)-DeltaLambda*np.tan(ThetaML))**2*8*log(2)/2/Ws**2 )

def AccMultWav(DeltaLambda, SigmaSource, DistanceFromSource, Ws, ThetaML):     #Multilayer wave acceptance
    return sqrt(6/pi) * exp( -(DeltaLambda)**2/2/((SigmaSource/DistanceFromSource)**2+Ws**2/8/log(2))/cotan(ThetaML)**2)


#Testing the functions


    #Testing the matrices

def testMatFlatMono():
    print(">> MatFlatMono in test")
    b=0.8
    ThetaB=1.0
    Mat=np.array([[0.8,0,0],[0,1.25,-0.389352],[0,0,1]])
    print(Mat.shape)
    result = MatFlatMono(b,ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatFlatMono())

def testMatCurvMono():
    print(">> MatCurvMono in test")
    b = 0.8
    Fc = 2
    ThetaB = 1.0
    Mat = np.array([[0.8, 0, 0], [0.5, 1.25, -0.389352], [0, 0, 1]])
    print(Mat.shape)
    result = MatCurvMono(b, Fc, ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatCurvMono())

def testMatMosMono():
    print(">> MatMosMono in test")
    ThetaB = 1.0
    Mat = np.array([[1, 0, 0], [0, -1, 2*np.tan(1)], [0, 0, 1]])
    print(Mat.shape)
    result = MatMosMono(ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatMosMono())

def testMatFlatMirr():
    print(">> MatFlatMirr in test")
    Sigma = 0.02
    IncAng = 30
    Lambda = 1
    Delta = 5
    Mat = np.array([[0.9402, 0, 0], [4.701, 0.9402, 0], [0, 0, 0.9402]])
    print(Mat.shape)
    result = MatFlatMirr(IncAng, Sigma, Lambda, Delta)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatFlatMirr())

def testMatMirr():
    print(">> MatMirr in test")
    Sigma = 0.02
    IncAng = 30
    Lambda = 1
    Delta = 5
    Fm = 10
    S =0.1
    Mat = np.array([[0.9402, 0, 0], [0.56412, 0.9402, 0], [0, 0, 0.9402]])
    print(Mat.shape)
    result = MatMirr(IncAng, Sigma, Lambda, Delta, Fm, S)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatMirr())

def testMatParaLens():
    print(">> MatParaLens in test")
    F=12.3
    Mat = np.array([[1, 0, 0], [0.0813008, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = MatParaLens(F)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatParaLens())

def testMatCircLens():
    print(">> MatCircLens in test")
    F = 12.3
    Coef = 0.4
    Mat = np.array([[1, 0, 0], [0.0325203, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = MatCircLens(Coef, F)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatCircLens())

def testMatMulti():
    print(">> MatMulti in test")
    Fmu = 15.1
    Mat = np.array([[1, 0, 0], [0.0662252, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = MatMulti(Fmu)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatMulti())

    #Testing the source

def testSourceX():
    print(">> SourceX in test")
    x = 0.5
    xp = 2
    Source = 0.605773
    result = SourceX(x, xp, SigmaSourceX=10, SigmaSourceXp=2)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))
print(testSourceX())

def testSourceY():
    print(">> SourceY in test")
    y = 0.5
    yp = 2
    Source = 0.753897
    result = SourceY(y, yp, SigmaSourceY=10, SigmaSourceYp=2, GammaSource=1)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))
print(testSourceY())

def testSourceLambda():
    print(">> SourceLambda in test")
    DeltaLambda = 4
    Source = 0.278037
    result = SourceLambda(DeltaLambda, SigmaSLambda=2.5)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))
print(testSourceLambda())

    #Testing the acceptances

def testAccSlit():
    print(">> AccSlit in test")
    y = 2
    Aperture = 75
    Acceptance = 1.19601
    result = AccSlit(y, Aperture, 0)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
    result2 = AccSlit(y, Aperture, 1)
    Acceptance2=0.867194
    print(Acceptance2, result2)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance2, result2, decimal=5))
print(testAccSlit())

def testAccPinh():
    print(">> AccPinh in test")
    y = 2
    Diameter = 75
    Acceptance = 1.58672
    result = AccPinh(y, Diameter)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccPinh())

def testAccMonoPlanAng():
    print(">> AccMonoPlanAng in test")
    RMono = 0.91
    RInt = 0.2
    Wd = 103.2
    yp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 0.00214474
    result = AccMonoPlanAng(RMono, RInt, Wd, yp, DeltaLambda, ThetaB)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMonoPlanAng())

def testAccMonoPlanWav():
    print(">> AccMonoPlanWav in test")
    Wd = 103.2
    SigmaYp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.21616
    result = AccMonoPlanWav(DeltaLambda, SigmaYp, ThetaB, Wd)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMonoPlanWav())

def testAccMonoMosaAng():
    print(">> AccMonoMosaAng in test")
    RMono = 0.91
    RInt = 0.2
    eta = 103.2
    yp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 0.00264987
    result = AccMonoMosaAng(RInt, eta, yp, DeltaLambda, ThetaB)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMonoMosaAng())

def testAccMonoMosaWav():
    print(">> AccMonoMosaWav in test")
    eta = 103.2
    SigmaYp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.36733
    result = AccMonoMosaWav(DeltaLambda, SigmaYp, eta, ThetaB)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMonoMosaWav())

def testAccMonoCurvAng():
    print(">> AccMonoCurvAng in test")
    RMono = 0.91
    RInt = 0.2
    Wd = 103.2
    yp = 0.001
    y = 8
    r = 9
    Alpha = 0.5
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 0.00218093
    result = AccMonoCurvAng(RMono, RInt, Wd, yp, y, r, ThetaB, Alpha, DeltaLambda)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMonoCurvAng())

def testAccMonoCurvWav():
    print(">> AccMonoCurvWav in test")
    Wd = 103.2
    SigmaSource = 5.5
    DistanceFromSource = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.38197
    result = AccMonoCurvWav(DeltaLambda, ThetaB, SigmaSource, DistanceFromSource, Wd)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMonoCurvWav())

def testAccLensPara():
    print(">> AccLensPara in test")
    x = 2
    Radius = 75
    Distance = 22.2
    Sigma = 33.3
    Acceptance = 0.471162
    result = AccLensPara(x, Radius, Distance, Sigma)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccLensPara())

def testAccLensCirc():
    print(">> AccLensCirc in test")
    x = 2
    Radius = 75
    Distance = 22.2
    FWHM = 43.3
    Sigma = 33.3
    Acceptance = 0.471079
    result = AccLensCirc(x, Sigma, FWHM, Radius, Distance)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccLensCirc())

def testAccMultAng():
    print(">> AccMultAng in test")
    Rml = 200
    yp = 0.03
    y = 22.2
    Rowland = 443.3
    ThetaML = 1
    DeltaLambda = 1
    Ws = 0.7
    Acceptance = 0.0000541429
    result = AccMultAng(Rml, yp, y, Rowland, ThetaML, DeltaLambda, Ws)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMultAng())

def testAccMultWav():
    print(">> AccMultWav in test")
    ThetaML = 1
    DeltaLambda = 1
    Ws = 0.7
    DistanceFromSource = 22
    SigmaSource = 45
    Acceptance = 1.04044
    result = AccMultWav(DeltaLambda, SigmaSource, DistanceFromSource, Ws, ThetaML)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAccMultWav())


#Single slit example
# number of elements in the device : n

z0 = 10000
z1 = 12000
Flight_paths_in_order = [z0,z1]
AperturSlitX0 = 2
AperturSlitY0 = 3
Slit = [AperturSlitX0, AperturSlitY0, np.ones((3, 3)), 0]
MatTab = [MatFlight(z0), Slit[2], MatFlight(z1)]

def MatCreation():
    MatTempX = np.array([x,xp,dl])
    MatTempY = np.array([y, yp, dl])
    for i in range(len(MatTab)-1,0,-1):
        MatTempX = np.dot(MatTab[i],MatTempX)
        MatTempY = np.dot(MatTab[i], MatTempY)
    return [MatTempX,MatTempY]
print(">>>>>>>" , MatCreation()[1], np.shape(MatCreation()[0]),np.shape(MatCreation()[1]))

def NewSource():
    MatTempX = MatCreation()[0]
    MatTempY = MatCreation()[1]
    NewSourceX = SourceX(MatTempX[0], MatTempX[1])
    NewSourceY = SourceY(MatTempY[0], MatTempY[1])
    return [NewSourceX, NewSourceY]
print(NewSource())

del MatTab[0:2]
MatCreation()

def NewSource2():
    NewSourceX = NewSource()[0]*AccSlit(MatCreation()[0][0], Slit[0], Slit[3])
    NewSourceY = NewSource()[1]*AccSlit(MatCreation()[1][0], Slit[1], Slit[3])
    return [NewSourceX, NewSourceY]
print(">>>>>>", NewSource2())


# Integrations
# Beam Size

def IntegralXXP(x, xp):
    NewSourceX = NewSource2()[0]
    if sp.diff(NewSourceX, dl) == 0:
        return 1
    else:
        return sp.integrate(NewSourceX, (dl, -np.inf, np.inf))
print(IntegralXXP(x, xp))

def IntegralYYP(y, yp):
    NewSourceY = NewSource2()[1]
    if sp.diff(NewSourceY, dl) == 0:
        return 1
    else:
        return sp.integrate(NewSourceY, (dl, -np.inf, np.inf))
print(IntegralYYP(y, yp))

def IntegralXYXPYP(x, xp, y, yp):
    SourceI = 6*10^20
    Beam = NewSource2()[0]*NewSource2()[1]*SourceLambda(dl, SigmaSLambda=1e-6)*SourceI
    if sp.diff(Beam, dl)==0:
        return 1
    else:
        return sp.integrate(Beam, (dl, -np.inf, np.inf))
# print(IntegralXYXPYP(x, xp, y, yp))

def IntegralXYXP(x,xp,y, n,    n_digits):
    if sp.diff(IntegralXYXPYP(x, xp, y, yp), yp) == 0:
        return 1
    else:
        v, w = gauss_hermite(n, n_digits)
        result = 0
        for i in range(n):
            result = result + w[i] * IntegralXYXPYP(x, xp, y, yp).subs(yp, v[i])
        return result
# print(IntegralXYXP(x, xp, y, 8, 3))

def IntegralXY(x, y, n, n_digits):
    if sp.diff(IntegralXYXP(x, xp, y, n, n_digits), xp) == 0:
        return 1
    else:
        v, w = gauss_hermite(n, n_digits)
        Int = IntegralXYXP(x, xp, y, n, n_digits)
        result = 0
        for i in range(n):
            result = result + w[i] * Int.subs(xp, v[i])
        return result
# print(IntegralXY(x, y, 5, 2))

def BeamGeoSize(n , n_digits):
    FunctionX = IntegralXY(x, y, n, n_digits).subs(y, 0)                     #todo check the minus in sigma
    FunctionY = IntegralXY(x, y, n, n_digits).subs(x, 0)
    ValueAX = FunctionX.subs(x, 0)
    ValueAY = FunctionY.subs(y, 0)
    ValueExponentX = FunctionX.subs(x, 10^-6)
    ValueExponentY = FunctionY.subs(y, 10^-6)
    SigmaX = 1/2/sqrt(2*np.log(2)) * sp.sqrt(4*log(2)*(10^-4)/sp.log(ValueExponentX/ValueAX))
    SigmaY = 1/2/sqrt(2*np.log(2)) * sp.sqrt(4*log(2)*(10^-4)/sp.log(ValueExponentY/ValueAY))
    return [SigmaX, SigmaY]
# print(BeamGeoSize())

# Angle size of the beam

def IntegralXXPYP(n, n_digits):
    if sp.diff(IntegralXYXPYP(x, xp, y, yp), y) == 0:
        return 1
    else:
        v, w = gauss_hermite(n, n_digits)
        result = 0
        for i in range(n):
            result = result + w[i] * IntegralXYXPYP(x, xp, y, yp).subs(y, v[i])
        return result
# print(IntegralXYXP(x, xp, y, 8, 3))

def IntegralXPYP(n, n_digits):
    if sp.diff(IntegralXXPYP(n, n_digits), x) == 0:
        return 1
    else:
        v, w = gauss_hermite(n, n_digits)
        result = 0
        for i in range(n):
            result = result + w[i] * IntegralXXPYP(n, n_digits).subs(x, v[i])
        return result
print(IntegralXPYP(5,2))

def AngularBeamSize(n, n_digits):
    FunctionXP = IntegralXPYP(n, n_digits).subs(yp,0)
    FunctionYP = IntegralXPYP(n, n_digits).subs(xp,0)
    ValueAXP = FunctionXP.subs(xp,0)
    ValueAYP = FunctionYP.subs(yp,0)
    ValueExponentXP = FunctionXP.subs(xp, 10^-6)
    ValueExponentYP = FunctionYP.subs(yp, 10^-6)
    SigmaXP = 1/2/sqrt(2*np.log(2)) * sp.sqrt(4*log(2)*(10^-4)/sp.log(ValueExponentXP/ValueAXP)/10^-12)
    SigmaYP = 1/2/sqrt(2*np.log(2)) * sp.sqrt(4*log(2)*(10^-4)/sp.log(ValueExponentYP/ValueAYP)/10^-12)
    return [SigmaXP, SigmaYP]
print(AngularBeamSize(4,2))