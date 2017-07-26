import numpy as np
from numpy import exp
from numpy import pi
from numpy import sqrt
from numpy import log
import numpy.testing as nut
import scipy.integrate as si
import mpmath as mp

#DEFINITION OF MATHEMATICAL OBJECTS
#Definition of the transformation matrices

def MatrixFlight(L):                                                   # For a flight path of length
    return np.array([[1,-L,0],[0,1,0],[0,0,1]])

def MatrixMonoPlane(b, ThetaB):                                         # For a perfect flat crystal monochromator
    return np.array([[b,0,0],[0,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatrixMonoBent(b, Fc, ThetaB):                                     # For a perfect curved crystal monochromator (meridionally and sagitally focusing)
    return np.array([[b,0,0],[1/Fc,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatrixMonoMosaic(ThetaB):                                             # For a mosaic monochromator
    return np.array([[1,0,0],[0,-1,2*np.tan(ThetaB)],[0,0,1]])

def MatrixMirrorPlane(IncAng, Sigma, Lambda, Delta):                      # For the plane mirror
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[Delta,1,0],[0,0,1]])

def MatrixMirror(IncAng, Sigma, Lambda, Delta, Fm, S):                      # For the bent and toroidal mirrors
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[(1+S*Fm*Delta)/Fm,1,0],[0,0,1]])

def MatrixCompLensPara(F):                                                 # For the compound refractive lenses with parabolic holes
    return np.array([[1,0,0],[1/F,1,0],[0,0,1]])

def MatrixCompLensCirc(Coef, F):                                           # For the compound refractive lenses with circular holes
    return np.array([[1,0,0],[Coef/F,1,0],[0,0,1]])

def MatrixMultilayer(Fmu):                                                  # For the multilayer
    return np.array([[1,0,0],[1/Fmu,1,0],[0,0,1]])

#Definition of the beam source

def SourceXXP(x, xp, SigmaXSource=1e-1, SigmaXPSource=3e-5):
    return exp(-( (x/SigmaXSource)**2 + (xp/SigmaXPSource)**2) / 2 )
def SourceYYP(y, yp, SigmaYSource=1e-2, SigmaYPSource=1e-4, GammaSource=0):
    return exp( -( (y/SigmaYSource)**2 + ((yp-GammaSource*y)/SigmaYPSource)**2)/2 )
def SourceLambda(dl, SigmaSLambda=1e-3):
    return exp(-(dl)**2/2/SigmaSLambda**2)

#Definition of the acceptance windows of the optical parts

def cotan(x):
    return 1/np.tan(x)

def AcceptanceSlit(y, Aperture, calctype):                                                           #Slit acceptance
    if calctype==0:
        return sqrt(6/pi) / sqrt(6*log(2)/pi) * exp( -(y/Aperture)**2/2*12)
    if calctype==1:
        return 1/sqrt(6*log(2)/pi)*exp(-y**2/(2*Aperture**2/2/pi))
    else:
        return 'the value for calctype is 0 or 1 you pipsqueak !'

def AcceptancePin(y, Diameter):                                                            #Pinhole acceptance
    return sqrt(8/pi) * exp ( -(y/Diameter)**2/2*16 )

def AcceptanceAngleMonoPlane(RMono, RInt, Wd, yp, DeltaLambda, ThetaB):                #Plane monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * exp( -(yp-DeltaLambda*np.tan(ThetaB))**2 / (2*Wd**2/12))

def AcceptanceWaveMonoPlane(DeltaLambda, SigmaYp, ThetaB, Wd):                        #Plane monochromator wave acceptance
    return sqrt(6/pi) * exp( - (DeltaLambda)**2 / (2*(SigmaYp**2+Wd**2/12)*cotan(ThetaB)**2) )

def AcceptanceAngleMonoMosaic(RInt, eta, yp, DeltaLambda, ThetaB):                       #Mosaic monochromator angular acceptance
    return RInt*sqrt(6/pi)/eta * exp(- (yp-DeltaLambda*np.tan(ThetaB))**2 / eta**2 /2 )

def AcceptanceWaveMonoMosaic(DeltaLambda, SigmaYp, eta, ThetaB):                        #Mosaic monochromator wave acceptance
    return sqrt(6/pi) * exp( - (DeltaLambda)**2 / 2 /((SigmaYp**2+eta**2)*cotan(ThetaB)**2))

def AcceptanceAngleMonoBent(RMono, RInt, Wd,yp, y, r, ThetaB, Alpha, DeltaLambda ):    #Curved monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * exp( - (yp-y/r/np.sin(ThetaB+Alpha)-DeltaLambda*np.tan(ThetaB))**2/(2*Wd**2/12))

def AcceptanceWaveMonoBent(DeltaLambda, ThetaB, SigmaSource, DistanceFromSource, Wd):  #Curved monochromator wave acceptance
    return sqrt(6/pi) * exp( -(DeltaLambda)**2 / (2*cotan(ThetaB)**2*((SigmaSource/DistanceFromSource)**2+Wd**2/12))   )

def AcceptanceCompLensPara(x, Radius, Distance, Sigma):                                          #Compound refractive lense with parabolic holes acceptance
    return exp( -(x**2+Radius*Distance)/2/Sigma**2)

def AcceptanceCompLensCirc(x, Sigma, FWHM, Radius, Distance):                                    #Compound refractive lense with circular holes acceptance
    return exp( -x**2/2/Sigma**2 -x**2*FWHM**2/8/Radius**2/Sigma**2 -x**2*FWHM**4/16/Sigma**2/Radius**4 -Radius*Distance/2/Sigma**2  )

def AcceptanceAngleMulti(Rml, yp, y, Rowland, ThetaML, DeltaLambda, Ws):                #Multilayer angle acceptance
    return Rml*8/3*sqrt(log(2)/pi) * exp( -(-yp-y/Rowland/np.sin(ThetaML)-DeltaLambda*np.tan(ThetaML))**2*8*log(2)/2/Ws**2 )

def AcceptanceWaveMulti(DeltaLambda, SigmaSource, DistanceFromSource, Ws, ThetaML):     #Multilayer wave acceptance
    return sqrt(6/pi) * exp( -(DeltaLambda)**2/2/((SigmaSource/DistanceFromSource)**2+Ws**2/8/log(2))/cotan(ThetaML)**2)

#TESTING THE FUNCTIONS
#Testing the matrices

def testMatrixMonoPlane():
    print(">> MatFlatMono in test")
    b=0.8
    ThetaB=1.0
    Mat=np.array([[0.8,0,0],[0,1.25,-0.389352],[0,0,1]])
    print(Mat.shape)
    result = MatrixMonoPlane(b,ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixMonoPlane())

def testMatrixMonoBent():
    print(">> MatCurvMono in test")
    b = 0.8
    Fc = 2
    ThetaB = 1.0
    Mat = np.array([[0.8, 0, 0], [0.5, 1.25, -0.389352], [0, 0, 1]])
    print(Mat.shape)
    result = MatrixMonoBent(b, Fc, ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixMonoBent())

def testMatrixMonoMosaic():
    print(">> MatMosMono in test")
    ThetaB = 1.0
    Mat = np.array([[1, 0, 0], [0, -1, 2*np.tan(1)], [0, 0, 1]])
    print(Mat.shape)
    result = MatrixMonoMosaic(ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixMonoMosaic())

def testMatrixMirrorPlane():
    print(">> MatFlatMirr in test")
    Sigma = 0.02
    IncAng = 30
    Lambda = 1
    Delta = 5
    Mat = np.array([[0.9402, 0, 0], [4.701, 0.9402, 0], [0, 0, 0.9402]])
    print(Mat.shape)
    result = MatrixMirrorPlane(IncAng, Sigma, Lambda, Delta)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixMirrorPlane())

def testMatrixMirror():
    print(">> MatMirr in test")
    Sigma = 0.02
    IncAng = 30
    Lambda = 1
    Delta = 5
    Fm = 10
    S =0.1
    Mat = np.array([[0.9402, 0, 0], [0.56412, 0.9402, 0], [0, 0, 0.9402]])
    print(Mat.shape)
    result = MatrixMirror(IncAng, Sigma, Lambda, Delta, Fm, S)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixMirror())

def testMatrixCompLensPara():
    print(">> MatParaLens in test")
    F=12.3
    Mat = np.array([[1, 0, 0], [0.0813008, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = MatrixCompLensPara(F)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixCompLensPara())

def testMatrixCompLensCirc():
    print(">> MatCircLens in test")
    F = 12.3
    Coef = 0.4
    Mat = np.array([[1, 0, 0], [0.0325203, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = MatrixCompLensCirc(Coef, F)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixCompLensCirc())

def testMatrixMultilayer():
    print(">> MatMulti in test")
    Fmu = 15.1
    Mat = np.array([[1, 0, 0], [0.0662252, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = MatrixMultilayer(Fmu)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatrixMultilayer())

#Testing the source

def testSourceXXP():
    print(">> SourceX in test")
    x = 0.5
    xp = 2
    Source = 0.605773
    result = SourceXXP(x, xp, SigmaXSource=10, SigmaXPSource=2)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))
print(testSourceXXP())

def testSourceYYP():
    print(">> SourceY in test")
    y = 0.5
    yp = 2
    Source = 0.753897
    result = SourceYYP(y, yp, SigmaYSource=10, SigmaYPSource=2, GammaSource=1)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))
print(testSourceYYP())

def testSourceLambda():
    print(">> SourceLambda in test")
    DeltaLambda = 4
    Source = 0.278037
    result = SourceLambda(DeltaLambda, SigmaSLambda=2.5)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))
print(testSourceLambda())

#Testing the acceptances

def testAcceptanceSlit():
    print(">> AccSlit in test")
    y = 2
    Aperture = 75
    Acceptance = 1.19601
    result = AcceptanceSlit(y, Aperture, 0)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
    result2 = AcceptanceSlit(y, Aperture, 1)
    Acceptance2=0.867194
    print(Acceptance2, result2)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance2, result2, decimal=5))
print(testAcceptanceSlit())

def testAcceptancePin():
    print(">> AccPinh in test")
    y = 2
    Diameter = 75
    Acceptance = 1.58672
    result = AcceptancePin(y, Diameter)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptancePin())

def testAcceptanceAngleMonoPlane():
    print(">> AccMonoPlanAng in test")
    RMono = 0.91
    RInt = 0.2
    Wd = 103.2
    yp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 0.00214474
    result = AcceptanceAngleMonoPlane(RMono, RInt, Wd, yp, DeltaLambda, ThetaB)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceAngleMonoPlane())

def testAcceptanceWaveMonoPlane():
    print(">> AccMonoPlanWav in test")
    Wd = 103.2
    SigmaYp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.21616
    result = AcceptanceWaveMonoPlane(DeltaLambda, SigmaYp, ThetaB, Wd)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceWaveMonoPlane())

def testAcceptanceAngleMonoMosaic():
    print(">> AccMonoMosaAng in test")
    RMono = 0.91
    RInt = 0.2
    eta = 103.2
    yp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 0.00264987
    result = AcceptanceAngleMonoMosaic(RInt, eta, yp, DeltaLambda, ThetaB)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceAngleMonoMosaic())

def testAcceptanceWaveMonoMosaic():
    print(">> AccMonoMosaWav in test")
    eta = 103.2
    SigmaYp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.36733
    result = AcceptanceWaveMonoMosaic(DeltaLambda, SigmaYp, eta, ThetaB)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceWaveMonoMosaic())

def testAcceptanceAngleMonoBent():
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
    result = AcceptanceAngleMonoBent(RMono, RInt, Wd, yp, y, r, ThetaB, Alpha, DeltaLambda)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceAngleMonoBent())

def testAcceptanceWaveMonoBent():
    print(">> AccMonoCurvWav in test")
    Wd = 103.2
    SigmaSource = 5.5
    DistanceFromSource = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.38197
    result = AcceptanceWaveMonoBent(DeltaLambda, ThetaB, SigmaSource, DistanceFromSource, Wd)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceWaveMonoBent())

def testAcceptanceCompLensPara():
    print(">> AccLensPara in test")
    x = 2
    Radius = 75
    Distance = 22.2
    Sigma = 33.3
    Acceptance = 0.471162
    result = AcceptanceCompLensPara(x, Radius, Distance, Sigma)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceCompLensPara())

def testAcceptanceCompLensCirc():
    print(">> AccLensCirc in test")
    x = 2
    Radius = 75
    Distance = 22.2
    FWHM = 43.3
    Sigma = 33.3
    Acceptance = 0.471079
    result = AcceptanceCompLensCirc(x, Sigma, FWHM, Radius, Distance)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceCompLensCirc())

def testAcceptanceAngleMulti():
    print(">> AccMultAng in test")
    Rml = 200
    yp = 0.03
    y = 22.2
    Rowland = 443.3
    ThetaML = 1
    DeltaLambda = 1
    Ws = 0.7
    Acceptance = 0.0000541429
    result = AcceptanceAngleMulti(Rml, yp, y, Rowland, ThetaML, DeltaLambda, Ws)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceAngleMulti())

def testAcceptanceWaveMulti():
    print(">> AccMultWav in test")
    ThetaML = 1
    DeltaLambda = 1
    Ws = 0.7
    DistanceFromSource = 22
    SigmaSource = 45
    Acceptance = 1.04044
    result = AcceptanceWaveMulti(DeltaLambda, SigmaSource, DistanceFromSource, Ws, ThetaML)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceWaveMulti())

# todo Simplest example : propagation

z0 = 20000
Flight_paths_in_order = [z0]
MatTab = [MatrixFlight(z0)]
dlBoundaries = [-10**-6,10**-6]
CoefMonoX = 1
CoefMonoY = 1
CoefAtten = 1
SourceI = 1e25

def MatCreation(x, xp, y, yp, dl):
    MatTempX = np.array([x,xp,dl])
    MatTempY = np.array([y, yp, dl])
    MatTempX = np.dot(MatrixFlight(z0),MatTempX)
    MatTempY = np.dot(MatrixFlight(z0), MatTempY)
    return [MatTempX,MatTempY]

def NewSource(x, xp, y, yp, dl):
    MatTempX = MatCreation(x, xp, y, yp, dl)[0]
    MatTempY = MatCreation(x, xp, y, yp, dl)[1]
    NewSourceX = SourceXXP(MatTempX[0], MatTempX[1])
    NewSourceY = SourceYYP(MatTempY[0], MatTempY[1])
    return [NewSourceX, NewSourceY]

def testNewSource(x, xp, y, yp, dl):
    print(">> NewSource in test")
    Math1 = sp.exp(1/2 * (-100*(x - 20000*xp)**2 - 1.11111*10**9*xp**2))
    Math2 = sp.exp(1/2*(-10000 *(y - 20000*yp)**2 - 1.*10**8*yp**2))
    Theory = Math1.subs([(x, 1), (xp, 2)])
    Result = NewSource(x, xp, y, yp , dl)[0].subs([(x, 1), (xp, 2), (y, 3), (yp, 4)])
    nut.assert_almost_equal(Theory, Result, decimal=5)
    Theory = Math2.subs([(y, 3), (yp, 4)])
    Result = NewSource(x, xp, y, yp,dl)[1].subs([(x, 1), (xp, 2), (y, 3), (yp, 4)])
    nut.assert_almost_equal(Theory, Result, decimal=5)
    print('Mistakes ? ')

def IxGeo(xp, dl, x):
    return NewSource(x, xp, 0, 0, dl)[0]
def IyGeo(yp, dl, y):
    return NewSource(0,0, y, yp, dl)[1]
def IXintGeo(xp, dl, x):
    return IxGeo(xp, dl, x)* SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
def IYintGeo(yp, dl, y):
    return IyGeo(yp, dl, y) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI

def BeamGeoSize():

    #integration limit calculation
    IotaXp = 10**-7
    while IxGeo(IotaXp, 0, 0) > 10**-10 * IxGeo(0,0,0) :
        IotaXp = IotaXp * 2
    IotaYp = 10**-7
    while IYintGeo(IotaYp, 0, 0) > 10**-10 * IYintGeo(0,0,0) :
        IotaYp = IotaYp * 2
    IotaYdl = 10**-7
    while IYintGeo(0, IotaYdl, 0) > 10**-10 * IYintGeo(0,0,0) :
        IotaYdl = IotaYdl * 2
    print('The integration limit is :', [IotaXp, IotaYp, IotaYdl])
    # we are going to do 4 integrations : along xp and dl with x = {0,10**-2}, along yp and dl with y = {0,10**-2}
    # we can do that since the integrals can be separated, and it takes less time to calculate. Adding a parameter
    # to a numerical integration multiplies the execution time by approx 300.
    # args(0,10**-2) sets the values of dl and x to 0 and 10**-2
    IxIntegrated_0 = si.quad(IxGeo, -IotaXp, IotaXp, args=(0, 0))[0]    #up is integrated
    IxIntegrated_E2 = si.quad(IxGeo, -IotaXp, IotaXp, args=(0, 10 ** -2))[0]
    IyIntegrated_0 = si.nquad(IYintGeo, [[-IotaYp, IotaYp],[-IotaYdl, IotaYdl]], args=(0,))[0]
    IyIntegrated_E2 = si.nquad(IYintGeo, [[-IotaYp, IotaYp],[-IotaYdl, IotaYdl]], args=(10 ** -2,))[0]
    print(IxIntegrated_0, IxIntegrated_E2, IyIntegrated_0, IyIntegrated_E2)
    ValueAX = IxIntegrated_0 * IyIntegrated_0
    print(ValueAX)
    ValueAY = ValueAX
    ValueExponentX = IyIntegrated_0 * IxIntegrated_E2
    ValueExponentY = IxIntegrated_0 * IyIntegrated_E2
    SigmaX = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) * (10 ** -4) / -log(ValueExponentX / ValueAX))
    SigmaY = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) * (10 ** -4) / -log(ValueExponentY / ValueAY))
    return [SigmaX, SigmaY]
print(BeamGeoSize())

def BeamAngularSize():
    Ixp = lambda x, dl, xp: NewSource(x, xp, 0, 0, dl)[0]
    Iyp = lambda y, dl, yp : NewSource(0,0, y, yp, dl)[1]
    IYpint = lambda y, dl, yp : Iyp(y, dl, yp) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
    #integration limit calculation
    IotaXp = 10**-7
    while Ixp(IotaXp, 0, 0) > 10**-10 * Ixp(0,0,0) :
        IotaXp = IotaXp * 2
    IotaYp = 10**-7
    while IYpint(IotaYp,0,0) > 10**-10 * IYpint(0,0,0):
        IotaYp = IotaYp * 2
    IotaYdl = 10**-7
    while IYpint(0, IotaYdl, 0) > 10**-10 * IYpint(0,0,0):
        IotaYdl = IotaYdl * 2
    print('Integration limits are :', [IotaXp,IotaYp, IotaYdl])
    #Integrations
    IxpIntegrated_0 = si.quad(Ixp, -IotaXp, IotaXp, args=(0, 0))[0]
    IxpIntegrated_E6 = si.quad(Ixp, -IotaXp, IotaXp, args=(0, 10 ** -6))[0]
    IypIntegrated_0 = si.nquad(IYpint, [[-IotaYp, IotaYp], [-IotaYdl, IotaYdl]], args=(0,))[0]
    IypIntegrated_E6 = si.nquad(IYpint, [[-IotaYp, IotaYp], [-IotaYdl, IotaYdl]], args=(10 ** -6,))[0]
    print(IxpIntegrated_0, IxpIntegrated_E6, IypIntegrated_0, IypIntegrated_E6)
    ValueAXp = IxpIntegrated_0 * IypIntegrated_0
    print(ValueAXp)
    ValueAYp = ValueAXp
    ValueExponentXp = IypIntegrated_0 * IxpIntegrated_E6
    ValueExponentYp = IxpIntegrated_0 * IypIntegrated_E6
    SigmaXp = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentXp / ValueAXp))
    SigmaYp = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentYp / ValueAYp))
    return [SigmaXp, SigmaYp]
print(BeamAngularSize())

def Sigma1_MaxFluxL_FluxPhi():
    Ix = lambda x, xp, dl: NewSource(x, xp, 0, 0, dl)[0]
    Iy = lambda y, yp, dl: NewSource(0, 0, y, yp, dl)[1]
    IYint = lambda y, yp, dl: Iy(y, yp, dl) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
    # integration limit calculation
    IotaX = 10 ** -7
    while Ix(IotaX, 0, 0) > 10**-20 * Ix(0, 0, 0):
        IotaX = IotaX * 2
    IotaXp = 10 ** -7
    while Ix(0, IotaXp, 0) > 10**-20 * Ix(0, 0, 0):
        IotaXp = IotaXp * 2
    IotaY = 10 ** -7
    while IYint(IotaY, 0, 0) > 10**-20 * IYint(0, 0, 0):
        IotaY = IotaY * 2
    IotaYp = 10 ** -7
    while IYint(0, IotaYp, 0) > 10**-20 * IYint(0, 0, 0):
        IotaYp = IotaYp * 2
    IotaYdl = 10 ** -7
    while IYint(0, 0, IotaYdl) > 10 ** -20 * IYint(0, 0, 0):
        IotaYdl = IotaYdl * 2
    # IotaX = 10*IotaX
    # IotaY = 10*IotaY
    # IotaYdl = 10*IotaYdl
    # IotaXp = 10*IotaXp
    # IotaYp = 10*IotaYp
    print('Integration limits are :', [IotaX, IotaY, IotaXp, IotaYp, IotaYdl])
    #Integrations
    IxIntegrated = si.nquad(Ix, [[-IotaX, IotaX], [-IotaXp, IotaXp]], args=(0,))[0]
    IyIntegrated_0 = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp, IotaYp]], args=(0,))[0]
    IyIntegrated_E3 = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp,IotaYp]], args=(10 ** -3,))[0]
    print('IxIntegrated is :', IxIntegrated)
    # ValueAL needs dl = 0, so we set it to 0 during the integration
    ValueAL = IxIntegrated * IyIntegrated_0
    print(" ValueAL gives :", ValueAL)
    # ValueExponentL needs dl = 10**-3, so we set it to 10**-3 during the integration
    ValueExponentL = IxIntegrated * IyIntegrated_E3
    print(" ValueExponentL gives :", ValueExponentL)
    Sigma = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentL / ValueAL) / 10 ** 6)
    MaxFluxL = ValueAL
    FluxPhi = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp, IotaYp], [-IotaYdl, IotaYdl]])
    print(FluxPhi)
    FluxPhi = IxIntegrated * FluxPhi[0]
    FluxPhi = CoefAtten * CoefMonoX * CoefMonoY * FluxPhi
    return Sigma, MaxFluxL, FluxPhi
print(Sigma1_MaxFluxL_FluxPhi())