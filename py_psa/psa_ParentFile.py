import numpy as np
from numpy import exp
from numpy import pi
from numpy import sqrt
from numpy import log
from numpy import sin
import numpy.testing as nut
import scipy.integrate as si
import mpmath as mp

#DEFINITION OF MATHEMATICAL OBJECTS
#Definition of the transformation matrices

def MatrixFlight(L):                                                   # For a flight path of length
    return np.array([[1,-L,0],[0,1,0],[0,0,1]])

def MatrixMonoPlane(b, ThetaB):                                         # For a perfect flat crystal monochromator
    return np.array([[b,0,0],[0,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatrixMonoBent(b, ThetaB, Fc):                                     # For a perfect curved crystal monochromator (meridionally and sagitally focusing)
    return np.array([[b,0,0],[1/Fc,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatrixMonoMosaic(ThetaB):                                             # For a mosaic monochromator
    return np.array([[1,0,0],[0,-1,2*np.tan(ThetaB)],[0,0,1]])

def MatrixMirrorPlane(Sigma,IncAng, Lambda, Delta):                      # For the plane mirror
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[Delta,1,0],[0,0,1]])

def MatrixMirror(Fm, Sigma, IncAng, Lambda, Delta, S):                      # For the bent and toroidal mirrors
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[(1+S*Fm*Delta)/Fm,1,0],[0,0,1]])

def MatrixCompLensPara(F):                                                 # For the compound refractive lenses with parabolic holes
    return np.array([[1,0,0],[1/F,1,0],[0,0,1]])

def MatrixCompLensCirc(F, Coef):                                           # For the compound refractive lenses with circular holes
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

def AcceptanceAngleMonoPlane(yp, DeltaLambda, ThetaB, Wd, RMono, RInt):                #Plane monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * exp( -(yp-DeltaLambda*np.tan(ThetaB))**2 / (2*Wd**2/12))

def AcceptanceWaveMonoPlane(DeltaLambda, SigmaYp, ThetaB, Wd):                        #Plane monochromator wave acceptance
    return sqrt(6/pi) * exp( - (DeltaLambda)**2 / (2*(SigmaYp**2+Wd**2/12)*cotan(ThetaB)**2) )

def AcceptanceAngleMonoMosaic(yp, DeltaLambda, ThetaB, eta, RInt):                       #Mosaic monochromator angular acceptance
    return RInt*sqrt(6/pi)/eta * exp(- (yp-DeltaLambda*np.tan(ThetaB))**2 / eta**2 /2 )

def AcceptanceWaveMonoMosaic(DeltaLambda, SigmaYp, ThetaB, eta):                        #Mosaic monochromator wave acceptance
    return sqrt(6/pi) * exp( - (DeltaLambda)**2 / 2 /((SigmaYp**2+eta**2)*cotan(ThetaB)**2))

def AcceptanceAngleMonoBent(y, yp, DeltaLambda, Alpha, ThetaB, Wd, r, RMono, RInt):    #Curved monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * exp( - (yp-y/r/np.sin(ThetaB+Alpha)-DeltaLambda*np.tan(ThetaB))**2/(2*Wd**2/12))

def AcceptanceWaveMonoBent(DeltaLambda, ThetaB, DistanceFromSource, Wd, SigmaSource):  #Curved monochromator wave acceptance
    return sqrt(6/pi) * exp( -(DeltaLambda)**2 / (2*cotan(ThetaB)**2*((SigmaSource/DistanceFromSource)**2+Wd**2/12))   )

def AcceptanceCompLensPara(x, Radius, Distance, Sigma):                                          #Compound refractive lense with parabolic holes acceptance
    return exp( -(x**2+Radius*Distance)/2/Sigma**2)

def AcceptanceCompLensCirc(x, Radius, Distance, Sigma, FWHM):                                    #Compound refractive lense with circular holes acceptance
    return exp( -x**2/2/Sigma**2 -x**2*FWHM**2/8/Radius**2/Sigma**2 -x**2*FWHM**4/16/Sigma**2/Radius**4 -Radius*Distance/2/Sigma**2  )

def AcceptanceAngleMulti(y, yp, DeltaLambda, ThetaML, Rml, Rowland, Ws):                #Multilayer angle acceptance
    return Rml*8/3*sqrt(log(2)/pi) * exp( -(-yp-y/Rowland/np.sin(ThetaML)-DeltaLambda*np.tan(ThetaML))**2*8*log(2)/2/Ws**2 )

def AcceptanceWaveMulti(DeltaLambda, ThetaML, DistanceFromSource, Ws, SigmaSource):     #Multilayer wave acceptance
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
    result = MatrixMonoBent(b, ThetaB, Fc)
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
    result = MatrixMirrorPlane(Sigma, IncAng, Lambda, Delta)
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
    result = MatrixMirror(Fm, Sigma, IncAng, Lambda, Delta, S)
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
    result = MatrixCompLensCirc(F, Coef)
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
    result = AcceptanceAngleMonoPlane(yp, DeltaLambda, ThetaB, Wd, RMono, RInt)
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
    result = AcceptanceAngleMonoMosaic(yp, DeltaLambda,ThetaB, eta, RInt)
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
    result = AcceptanceWaveMonoMosaic(DeltaLambda, SigmaYp, ThetaB, eta)
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
    result = AcceptanceAngleMonoBent(y, yp, DeltaLambda, Alpha, ThetaB, Wd, r, RMono, RInt)
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
    result = AcceptanceWaveMonoBent(DeltaLambda, ThetaB, DistanceFromSource, Wd, SigmaSource)
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
    result = AcceptanceCompLensCirc(x, Radius, Distance, Sigma, FWHM)
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
    result = AcceptanceAngleMulti(y, yp, DeltaLambda, ThetaML, Rml, Rowland, Ws)
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
    result = AcceptanceWaveMulti(DeltaLambda, ThetaML, DistanceFromSource, Ws, SigmaSource)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
print(testAcceptanceWaveMulti())



#Calculation part

z0 = 10000
z1 = 20000
Flight_paths_in_order = [z0, z1]

# Slit data
AperturSlitX0 = 20
AperturSlitY0 = 30
Slit = [AperturSlitX0, AperturSlitY0, np.eye(3), 0]

# Mono data
Alpha0  = 10**-15
Wd0     = 1.000000*10**-6
WidthX0 = 2.000000*10**1
WidthY0 = 2.000000*10**1
RMono0  = 9.600000*10**-1
Rint0   = 1.000000*10**-5
thetaB0 = 1.140288*10**1 * pi /180
bMono0  = sin(thetaB0+Alpha0)/sin(thetaB0-Alpha0)
Mono = [Alpha0, Wd0, WidthX0, WidthY0,RMono0, Rint0, thetaB0, bMono0, MatrixMonoPlane(bMono0, thetaB0), np.eye(3)]

# # Mirror data
Sigma0 = 0
Delta0 = 1.000000*10**-7
IncAng0 = 0
Lambda0 = 1.127140*10**-7
Fm0 = 8.333333*10**3
S = 1
Mirror = [np.eye(3), MatrixMirror(Fm0, Sigma0, IncAng0, Lambda0, Delta0, S)]

ListObject = ['MonoPlaneVertical']
MatTabX = [MatrixFlight(z0), Mono[-1], MatrixFlight(z1)]
MatTabY = [MatrixFlight(z0), Mono[-2], MatrixFlight(z1)]
# MatTabX = [MatrixFlight(z0)]
# MatTabY = [MatrixFlight(z0)]
SourceI = 1e25
SigmaYPSource = 10**-4
SigmaXPSource = 3e-5
bMonoX = 1
bMonoY = 1
CoefMonoX = 1
CoefMonoY = 1
CoefAtten = 1

# class Object:
#     def __init__(self):
#         self.type = "Mirror"
#     Object.Matrix = MatrixMirror



def SourceCreation(x, xp, y, yp, dl):
    # initiating variables, copies of the arrays are created
    MX = MatTabX.copy()
    MY = MatTabY.copy()
    MatTempX = np.array([x,xp,dl])
    MatTempY = np.array([y, yp, dl])
    # the case of single propagation has to be done separately because of a pb of range
    if len(MX)==1:
        MatTempX = np.dot(MX[0], MatTempX)
        MatTempY = np.dot(MY[0], MatTempY)
        NewSourceX = SourceXXP(MatTempX[0], MatTempX[1])
        NewSourceY = SourceYYP(MatTempY[0], MatTempY[1])
    # now we do the general case
    else :
        #first matrix multiplication
        for i in range(len(MX)-1, -1, -1):
            MatTempX = np.dot(MX[i], MatTempX)
            MatTempY = np.dot(MY[i], MatTempY)
        NewSourceX = SourceXXP(MatTempX[0], MatTempX[1])
        NewSourceY = SourceYYP(MatTempY[0], MatTempY[1])
        del MX[0]
        del MY[0]
        k = 0
        #we are going to do our matrix product, then apply the acceptance if needed and calculate the new resulting
        # source. Once this is done, we erase the two first matrices and do this again until our arrays are empty
        while MX!=[] and MY!=[]:
            for i in range(len(MX)-1,-1,-1):
                MatTempX = np.array([x, xp, dl])
                MatTempY = np.array([y, yp, dl])
                MatTempX = np.dot(MX[i],MatTempX)
                MatTempY = np.dot(MY[i], MatTempY)
            if ListObject[k] == 'Slit':
                NewSourceX = NewSourceX * AcceptanceSlit(MatTempX[0], Slit[0], Slit[3])
                NewSourceY = NewSourceY * AcceptanceSlit(MatTempY[0], Slit[1], Slit[3])
            elif ListObject[k] == 'MonoPlaneVertical':
                NewSourceY = NewSourceY * AcceptanceAngleMonoPlane(MatTempY[1], MatTempY[2], thetaB0, Wd0, RMono0, Rint0)
                NewSourceY = NewSourceY * AcceptanceWaveMonoPlane(MatTempY[2], bMonoY*SigmaYPSource, thetaB0, Wd0)
            elif ListObject[k] == 'MonoPlaneHorizontal':
                NewSourceX = NewSourceX * AcceptanceAngleMonoPlane(MatTempX[1], MatTempX[2], thetaB0, Wd0, RMono0, Rint0)
                NewSourceX = NewSourceX * AcceptanceWaveMonoPlane(MatTempX[2], bMonoX*SigmaXPSource, thetaB0, Wd0)
            k = k +1
            del MX[0:2]
            del MY[0:2]
    return [NewSourceX, NewSourceY]

def BeamGeoSize():  # the two next functions are similar to this model
    # creating the functions to integrate, we need to redo this every time because the integration functions integrate
    # using the order in which the parameters are input, like it will integrate Ix in this order : xp, dl and then x
    Ix = lambda xp, dl, x : SourceCreation(x, xp, 0, 0, dl)[0]
    Iy = lambda yp, dl, y : SourceCreation(0,0, y, yp, dl)[1]
    # integration limit calculation. We need to minimize calculation time and maximize precision and avoid empty
    # sampling. So we determine the boundaries of our gaussians with a quick algorithm
    IotaXp = 10 ** -8
    while Ix(IotaXp, 0, 0) > 10 ** -15 * Ix(0, 0, 0):
        IotaXp = IotaXp * 2
    IotaXp = IotaXp * 5               # this line is necessary, we suffer a great lack of precision sometimes without it
    IotaYp = 10 ** -8
    while Iy(IotaYp, 0, 0) > 10 ** -15 * Iy(0, 0, 0):
        IotaYp = IotaYp * 2
    IotaYp = IotaYp * 5
    # Again, to minimize integration, we calculate the integrals of I = Ix*Iy*Idl separately, do do that, we need to
    # that either Ix or Iy does not depend on dl, thus the next test
    if Ix(0,0,0) == Ix(0,100,0): #we check the dependancy on dl, being similar to a gaussian (bijectivity of a gaussian
        # on the interval [0, +inf[ )
        print("Lambda is affected to y")
        #we need to calculate one last integration boundary and create our integration function
        IYint = lambda yp, dl, y : Iy(yp, dl, y) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
        IotaYdl = 10**-8
        while IYint(0, IotaYdl, 0) > 10**-15 * IYint(0,0,0) :
            IotaYdl = IotaYdl * 2
        print('The integration limits are :', [IotaXp, IotaYp, IotaYdl])
        # we are going to do 4 integrations : along xp and dl with x = {0,10**-2}, along yp and dl with y = {0,10**-2}
        # we can do that since the integrals can be separated, and it takes less time to calculate. Adding a parameter
        # to a numerical integration multiplies the execution time by approx 300.
        # args(0,10**-2) sets the values of dl and x to 0 and 10**-2
        XEval = 10**-2
        YEval = 10**-2
        IxIntegrated_0 = si.quad(Ix, -IotaXp, IotaXp, args=(0, 0))[0]    #xp is integrated
        IxIntegrated_E2 = si.quad(Ix, -IotaXp, IotaXp, args=(0, XEval))[0]
        IyIntegrated_0 = si.nquad(IYint, [[-IotaYp, IotaYp],[-IotaYdl, IotaYdl]], args=(0,))[0]
        IyIntegrated_E2 = si.nquad(IYint, [[-IotaYp, IotaYp],[-IotaYdl, IotaYdl]], args=(YEval,))[0]
        #integral denullification section : Python does not have an infinite precision and sometimes the integrals are
        # too small to be calculated and so the result is 0, which is not acceptable. So we calculate them again using
        # a different evaluation
        if IxIntegrated_E2 == 0.0:
            while IxIntegrated_E2 == 0.0:
                XEval = XEval / 2
                IxIntegrated_E2 = si.quad(Ix, -IotaXp, IotaXp, args=(0, XEval))[0]
            XEval = XEval / 4
            IxIntegrated_E2 = si.quad(Ix, -IotaXp, IotaXp, args=(0, XEval))[0]
        if IyIntegrated_E2 == 0.0:
            while IyIntegrated_E2 == 0.0:
                YEval = YEval / 2
                IyIntegrated_E2 = si.nquad(IYint, [[-IotaYp, IotaYp], [-IotaYdl, IotaYdl]], args=(YEval,))[0]
            YEval = YEval / 4
            IyIntegrated_E2 = si.nquad(IYint, [[-IotaYp, IotaYp], [-IotaYdl, IotaYdl]], args=(YEval,))[0]
        print(IxIntegrated_0, IxIntegrated_E2, IyIntegrated_0, IyIntegrated_E2)
        # we now have our integrals, we just need both Sigmas
        ValueAX = IxIntegrated_0 * IyIntegrated_0
        print(ValueAX)
        ValueAY = ValueAX
        ValueExponentX = IyIntegrated_0 * IxIntegrated_E2
        ValueExponentY = IxIntegrated_0 * IyIntegrated_E2
        print('ValueExponentY is :', ValueExponentY)
        SigmaX = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) * XEval**2 / -log(ValueExponentX / ValueAX))
        SigmaY = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) * YEval**2 / -log(ValueExponentY / ValueAY))
        return [SigmaX, SigmaY]
    elif Iy(0,0,0) == Iy(0,100,0):
        # If we arrive here, this means Ix depended of dl, so we do exactly the same thing as in the previous paragraph,
        # but with Iy
        print("Lambda is affected to x")
        IXint = lambda xp, dl, x: Ix(xp, dl, x) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
        IotaXdl = 10 ** -8
        while IXint(0, IotaXdl, 0) > 10 ** -15 * IXint(0, 0, 0):
            IotaXdl = IotaXdl * 2
        IotaXdl = IotaXdl * 5
        print('The integration limits are :', [IotaXp, IotaYp, IotaXdl])
        # integrations
        XEval = 10**-2
        YEval = 10**-2
        IyIntegrated_0 = si.quad(Iy, -IotaYp, IotaYp, args=(0, 0))[0]
        IyIntegrated_E2 = si.quad(Iy, -IotaYp, IotaYp, args=(0, YEval))[0]
        IxIntegrated_0 = si.nquad(IXint, [[-IotaXp, IotaXp], [-IotaXdl, IotaXdl]], args=(0,))[0]
        IxIntegrated_E2 = si.nquad(IXint, [[-IotaXp, IotaXp], [-IotaXdl, IotaXdl]], args=(XEval,))[0]
        # integral denullification section
        if IxIntegrated_E2 == 0.0:
            while IxIntegrated_E2 == 0.0:
                XEval = XEval / 2
                IxIntegrated_E2 = si.nquad(IXint, [[-IotaXp, IotaXp], [-IotaXdl, IotaXdl]], args=(XEval,))[0]
            XEval = XEval / 4
            IxIntegrated_E2 = si.nquad(IXint, [[-IotaXp, IotaXp], [-IotaXdl, IotaXdl]], args=(XEval,))[0]
        if IyIntegrated_E2 == 0.0:
            while IyIntegrated_E2 == 0.0:
                YEval = YEval / 2
                IyIntegrated_E2 = si.quad(Iy, -IotaYp, IotaYp, args=(0, YEval))[0]
            YEval = YEval / 4
            IyIntegrated_E2 = si.quad(Iy, -IotaYp, IotaYp, args=(0, YEval))[0]
        print(IyIntegrated_0, IyIntegrated_E2, IxIntegrated_0, IxIntegrated_E2)
        #simple calculation part
        ValueAX = IxIntegrated_0 * IyIntegrated_0
        print(ValueAX)
        ValueAY = ValueAX
        ValueExponentX = IyIntegrated_0 * IxIntegrated_E2
        ValueExponentY = IxIntegrated_0 * IyIntegrated_E2
        SigmaX = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) * XEval**2 / -log(ValueExponentX / ValueAX))
        SigmaY = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) * YEval**2 / -log(ValueExponentY / ValueAY))
        return [SigmaX, SigmaY]
    else:
        print("Computation time too long, DeltaLambda variation on both axis -> not possible")
        return 0
print(BeamGeoSize())

def BeamAngularSize():
    #creating the functions to integrate
    Ixp = lambda x, dl, xp: SourceCreation(x, xp, 0, 0, dl)[0]
    Iyp = lambda y, dl, yp : SourceCreation(0,0, y, yp, dl)[1]
    # integration limit calculation
    IotaX = 10 ** -8
    while Ixp(IotaX, 0, 0) > 10 ** -15 * Ixp(0, 0, 0):
        IotaX = IotaX * 2
    IotaX = IotaX * 10
    IotaY = 10 ** -8
    while Iyp(IotaY, 0, 0) > 10 ** -15 * Iyp(0, 0, 0):
        IotaY = IotaY * 2
    IotaY = IotaY * 10
    #dl dependancy check on Ixp
    if Ixp(0, 0, 0) == Ixp(0, 100, 0):
        print("Lambda is affected to y")
        IYpint = lambda y, dl, yp : Iyp(y, dl, yp) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
        IotaYdl = 10**-8
        while IYpint(0, IotaYdl, 0) > 10**-15 * IYpint(0,0,0):
            IotaYdl = IotaYdl * 2
        IotaYdl = IotaYdl * 10
        print('Integration limits are :', [IotaX,IotaY, IotaYdl])
        #Integrations
        XpEval = 10**-6
        YpEval = 10**-6
        IxpIntegrated_0 = si.quad(Ixp, -IotaX, IotaX, args=(0, 0))[0]
        IxpIntegrated_E6 = si.quad(Ixp, -IotaX, IotaX, args=(0, XpEval))[0]
        IypIntegrated_0 = si.nquad(IYpint, [[-IotaY, IotaY], [-IotaYdl, IotaYdl]], args=(0,))[0]
        IypIntegrated_E6 = si.nquad(IYpint, [[-IotaY, IotaY], [-IotaYdl, IotaYdl]], args=(YpEval,))[0]
        # integral denullification section
        if IxpIntegrated_E6 == 0.0:
            while IxpIntegrated_E6 == 0.0:
                XpEval = XpEval / 2
                IxpIntegrated_E6 = si.quad(Ixp, -IotaX, IotaX, args=(0, XpEval))[0]
            XpEval = XpEval / 4
            IxpIntegrated_E6 = si.quad(Ixp, -IotaX, IotaX, args=(0, XpEval))[0]
        if IypIntegrated_E6 == 0:
            while IypIntegrated_E6 == 0:
                YpEval = YpEval / 2
                IypIntegrated_E6 = si.nquad(IYpint, [[-IotaY, IotaY], [-IotaYdl, IotaYdl]], args=(YpEval,))[0]
            YpEval = YpEval / 4
            IypIntegrated_E6 = si.nquad(IYpint, [[-IotaY, IotaY], [-IotaYdl, IotaYdl]], args=(YpEval,))[0]
        print(IxpIntegrated_0, IxpIntegrated_E6, IypIntegrated_0, IypIntegrated_E6)
        #Simple calculation part
        ValueAXp = IxpIntegrated_0 * IypIntegrated_0
        print(ValueAXp)
        ValueAYp = ValueAXp
        ValueExponentXp = IypIntegrated_0 * IxpIntegrated_E6
        ValueExponentYp = IxpIntegrated_0 * IypIntegrated_E6
        SigmaXp = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentXp / ValueAXp) * XpEval**2 )
        SigmaYp = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentYp / ValueAYp) * YpEval**2 )
        SigmaXpMilliRad = SigmaXp * 10**6
        SigmaYpMilliRad = SigmaYp * 10**6
        return [SigmaXpMilliRad, SigmaYpMilliRad]
    #Dl dependancy check on Iyp
    elif Iyp(0,0,0) == Iyp(0,100,0):
        print("Lambda is affected to x")
        IXpint = lambda x, dl, xp: Ixp(x, dl, xp) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
        IotaXdl = 10 ** -8
        while IXpint(0, IotaXdl, 0) > 10 ** -15 * IXpint(0, 0, 0):
            IotaXdl = IotaXdl * 2
        IotaXdl = IotaXdl * 10
        print('Integration limits are :', [IotaX, IotaY, IotaXdl])
        # Integrations
        XpEval = 10**-6
        YpEval = 10**-6
        IypIntegrated_0 = si.quad(Iyp, -IotaY, IotaY, args=(0, 0))[0]
        IypIntegrated_E6 = si.quad(Iyp, -IotaY, IotaY, args=(0, YpEval))[0]
        IxpIntegrated_0 = si.nquad(IXpint, [[-IotaX, IotaX], [-IotaXdl, IotaXdl]], args=(0,))[0]
        IxpIntegrated_E6 = si.nquad(IXpint, [[-IotaX, IotaX], [-IotaXdl, IotaXdl]], args=(XpEval,))[0]
        # integral denullification section
        if IxpIntegrated_E6 == 0.0:
            while IxpIntegrated_E6 == 0.0:
                XpEval = XpEval / 2
                IxpIntegrated_E6 = si.nquad(IXpint, [[-IotaX, IotaX], [-IotaXdl, IotaXdl]], args=(XpEval,))[0]
            XpEval = XpEval / 4
            IxpIntegrated_E6 = si.quad(Ixp, -IotaX, IotaX, args=(0, XpEval))[0]
        if IypIntegrated_E6 == 0:
            while IypIntegrated_E6 == 0:
                YpEval = YpEval / 2
                IypIntegrated_E6 = si.quad(Iyp, -IotaY, IotaY, args=(0, YpEval))[0]
            YpEval = YpEval / 4
            IypIntegrated_E6 = si.quad(Iyp, -IotaY, IotaY, args=(0, YpEval))[0]
        print(IxpIntegrated_0, IxpIntegrated_E6, IypIntegrated_0, IypIntegrated_E6)
        # Simple calculation part
        ValueAXp = IxpIntegrated_0 * IypIntegrated_0
        print(ValueAXp)
        ValueAYp = ValueAXp
        ValueExponentXp = IypIntegrated_0 * IxpIntegrated_E6
        ValueExponentYp = IxpIntegrated_0 * IypIntegrated_E6
        SigmaXp = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentXp / ValueAXp) * XpEval**2 )
        SigmaYp = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentYp / ValueAYp) * YpEval**2 )
        SigmaXpMilliRad = SigmaXp * 10**6
        SigmaYpMilliRad = SigmaYp * 10**6
        return [SigmaXpMilliRad, SigmaYpMilliRad]
    else:
        print("Computation time too long, DeltaLambda variation on both axis -> not possible")
        return 0
print(BeamAngularSize())

def Sigma1_MaxFluxL_FluxPhi():
    #creating the functions to integrate
    Ix = lambda x, xp, dl: SourceCreation(x, xp, 0, 0, dl)[0]
    Iy = lambda y, yp, dl: SourceCreation(0, 0, y, yp, dl)[1]
    #integration limit calculation
    IotaX = 10 ** -8
    while Ix(IotaX, 0, 0) > 10 ** -15 * Ix(0, 0, 0):
        IotaX = IotaX * 2
    IotaX = IotaX * 2
    IotaXp = 10 ** -8
    while Ix(0, IotaXp, 0) > 10 ** -15 * Ix(0, 0, 0):
        IotaXp = IotaXp * 2
    IotaXp = IotaXp * 2
    IotaY = 10 ** -8
    while Iy(IotaY, 0, 0) > 10 ** -15 * Iy(0, 0, 0):
        IotaY = IotaY * 2
    IotaY = IotaY * 2
    IotaYp = 10 ** -8
    while Iy(0, IotaYp, 0) > 10 ** -15 * Iy(0, 0, 0):
        IotaYp = IotaYp * 2
    IotaYp = IotaYp * 2
    if Ix(0, 0, 0) == Ix(0, 0, 100):       #first dependancy check
        print("Lambda is affected to y")
        # last integration boundary needed and function creation
        IYint = lambda y, yp, dl: Iy(y, yp, dl) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
        IotaYdl = 10 ** -8
        while IYint(0, 0, IotaYdl) > 10 ** -15 * IYint(0, 0, 0):
            IotaYdl = IotaYdl * 2
        IotaYdl = IotaYdl * 10
        print('Integration limits are :', [IotaX, IotaY, IotaXp, IotaYp, IotaYdl])
        #Integrations
        DlEval = 10**-3
        IxIntegrated = si.nquad(Ix, [[-IotaX, IotaX], [-IotaXp, IotaXp]], args=(0,))[0]
        IyIntegrated_0 = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp, IotaYp]], args=(0,))[0]
        IyIntegrated_E3 = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp,IotaYp]], args=(DlEval,))[0]
        # integral denullification section
        if IyIntegrated_E3 == 0.0:
            while IyIntegrated_E3 == 0.0:
                DlEval = DlEval / 2
                IyIntegrated_E3 = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp, IotaYp]], args=(DlEval,))[0]
            DlEval = DlEval / 4
            IyIntegrated_E3 = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp, IotaYp]], args=(DlEval,))[0]
        print('Lambda will be evaluated at :', DlEval)
        print('IxIntegrated is :', IxIntegrated, 'Other values :', IyIntegrated_0,IyIntegrated_E3 )
        #simple calculation part
        ValueAL = IxIntegrated * IyIntegrated_0
        print(" ValueAL gives :", ValueAL)
        ValueExponentL = IxIntegrated * IyIntegrated_E3
        print(" ValueExponentL gives :", ValueExponentL)
        Sigma = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentL / ValueAL) * DlEval**2)
        MaxFluxL = ValueAL
        # calculation of the flux
        FluxPhi = si.nquad(IYint, [[-IotaY, IotaY], [-IotaYp, IotaYp], [-IotaYdl, IotaYdl]])
        print(FluxPhi)
        FluxPhi = IxIntegrated * FluxPhi[0]
        FluxPhi = CoefAtten * CoefMonoX * CoefMonoY * FluxPhi
        return Sigma, MaxFluxL, FluxPhi
    elif Iy(0, 0, 0) == Iy(0, 0, 100):
        IXint = lambda x, xp, dl: Ix(x, xp, dl) * SourceLambda(dl, SigmaSLambda=1e-3) * SourceI
        IotaXdl = 10 ** -7
        while IXint(0, 0, IotaXdl) > 10 ** -15 * IXint(0, 0, 0):
            IotaXdl = IotaXdl * 2
        IotaXdl = IotaXdl * 10
        print('Integration limits are :', [IotaX, IotaY, IotaXp, IotaYp, IotaXdl])
        # Integrations
        DlEval = 10 ** -3
        IyIntegrated = si.nquad(Iy, [[-IotaY, IotaY], [-IotaYp, IotaYp]], args=(0,))[0]
        IxIntegrated_0 = si.nquad(IXint, [[-IotaX, IotaX], [-IotaXp, IotaXp]], args=(0,))[0]
        IxIntegrated_E3 = si.nquad(IXint, [[-IotaX, IotaX], [-IotaXp, IotaXp]], args=(DlEval,))[0]
        # integral denullification section
        if IxIntegrated_E3 == 0.0:
            while IxIntegrated_E3 == 0.0:
                DlEval = DlEval / 2
                IxIntegrated_E3 = si.nquad(IXint, [[-IotaX, IotaX], [-IotaXp, IotaXp]], args=(DlEval,))[0]
            DlEval = DlEval / 4
            IxIntegrated_E3 = si.nquad(IXint, [[-IotaX, IotaX], [-IotaXp, IotaXp]], args=(DlEval,))[0]
        print('Lambda will be evaluated at :', DlEval)
        print('IxIntegrated is :', IxIntegrated_0)
        # simple calculation part
        ValueAL = IxIntegrated_0 * IyIntegrated
        print(" ValueAL gives :", ValueAL)
        ValueExponentL = IyIntegrated * IxIntegrated_E3
        print(" ValueExponentL gives :", ValueExponentL)
        Sigma = 1 / 2 / sqrt(2 * np.log(2)) * sqrt(4 * log(2) / -log(ValueExponentL / ValueAL) * DlEval**2 )
        MaxFluxL = ValueAL
        # calculation of the flux
        FluxPhi = si.nquad(IXint, [[-IotaX, IotaX], [-IotaXp, IotaXp], [-IotaXdl, IotaXdl]])
        print(FluxPhi)
        FluxPhi = IyIntegrated * FluxPhi[0]
        FluxPhi = CoefAtten * CoefMonoX * CoefMonoY * FluxPhi
        return Sigma, MaxFluxL, FluxPhi
    else:
        print("Computation time too long, DeltaLambda variation on both axis -> not possible")
        return 0
print(Sigma1_MaxFluxL_FluxPhi())




# i dont really know what these functions are for
#
# # Integrations
# # Beam Size
#
# # def IntegralXXP():
# #     NewSourceX = NewSource2()[0]
# #     if sp.diff(NewSourceX, dl) == 0:
# #         return 1
# #     else:
# #         return sp.integrate(NewSourceX, (dl, -np.inf, np.inf))
# # print(IntegralXXP())
# #
# # def IntegralYYP():
# #     NewSourceY = NewSource2()[1]
# #     if sp.diff(NewSourceY, dl) == 0:
# #         return 1
# #     else:
# #         return sp.integrate(NewSourceY, (dl, -np.inf, np.inf))
# # print(IntegralYYP())
