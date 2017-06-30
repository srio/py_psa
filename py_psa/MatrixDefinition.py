import numpy as np
from numpy import exp
from numpy import pi
from numpy import sqrt
from numpy import log
import numpy.testing as nut

#Definition of the transformation matrices

def MatFlight(L):                                                   # For a flight path of length
    return np.array([[1,-L,0],[0,1,0],[0,0,1]])

def MatFirstMono(c):
    return np.array([[c,0,0],[0,1,0],[0,0,1]])                    #For the first monochromator

def MatFlatMono(b,c, ThetaB):                                         # For a perfect flat crystal monochromator
    return np.array([[b*c,0,0],[0,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatCurvMono(b, Fc, ThetaB, c):                                     # For a perfect curved crystal monochromator (meridionally and sagitally focusing)
    return np.array([[b*c,0,0],[1/Fc,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def MatMosMono(ThetaB, c):                                             # For a mosaic monochromator
    return np.array([[c,0,0],[0,-1,2*np.tan(ThetaB)],[0,0,1]])

def MatFlatMirr(IncAng, Sigma, Lambda, Delta):                      # For the plane mirror
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[Delta,1,0],[0,0,1]])

def MatMirr(IncAng, Sigma, Lambda, Delta, Fm):                      # For the bent and toroidal mirrors
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[(1+S*Fm*Delta)/Fm,1,0],[0,0,1]])

def MatParaLens(F):                                                 # For the compound refractive lenses with parabolic holes
    return np.array([[1,0,0],[1/F,1,0],[0,0,1]])

def MatCircLens(Coef, F):                                           # For the compound refractive lenses with circular holes
    return np.array([[1,0,0],[Coef/F,1,0],[0,0,1]])

def MatMulti(Fmu):                                                  # For the multilayer
    return np.array([[1,0,0],[exp(1)/Fmu,1,0],[0,0,1]])


#Definition of the beam source

def SourceX(x, SigmaSourceX, xp, SigmaSourceXp):
    return exp(-( (x/SigmaSourceX)**2 + (xp/SigmaSourceXp)**2) / 2 )
def SourceY(y, SigmaSourceY, yp, GammaSource, SigmaSourceYp):
    return exp( -( (y/SigmaSourceY)**2 + ((yp-GammaSource*y)/SigmaSourceYp)**2)/2 )
def SourceLambda(DLambda, Lambda, SigmaSourceLambda):
    return exp(-(DLambda/Lambda)**2/2/SigmaSourceLambda**2)

#Definition of the acceptance windows of the optical parts

def cotan(x):
    return 1/np.tan(x)

    #Slit acceptance
def AccSlit(y, Aperture):
    return sqrt(6/pi) / sqrt(6*log(2)/pi) * exp( -(y/Aperture)**2/2*12) #todo check sqrt(6*log(2)/pi)

    #Pinhole acceptance
def AccPinh(y, Diameter):
    return sqrt(8/pi) * exp ( -(y/Diameter)**2/2*16 )

    #Plane monochromator angle acceptance
def AccMonoPlanAng(RMono, RInt, Wd, yp, DeltaLambda, Lambda, ThetaB):
    return RMono*RInt*sqrt(6/pi)/Wd * exp( -(yp-DeltaLambda/Lambda*np.tan(ThetaB))**2 / (2*Wd**2/12))

    #Plane monochromator wave acceptance
def AccMonoPlanWav(DeltaLambda, Lambda, SigmaYp, ThetaB, Wd):
    return sqrt(6/pi) * exp( - (DeltaLambda/Lambda)**2 / (2*(SigmaYp**2+Wd**2/12)*cotan(ThetaB)**2) )

    #Mosaic monochromator angular acceptance
def AccMonoMosaAng(RInt, eta, yp, DeltaLambda, Lambda, ThetaB):
    return RInt*sqrt(6/pi)/eta * exp(- (yp-DeltaLambda/Lambda*np.tan(ThetaB))**2 / eta**2 /2 ) # todo : check /2 and RInt

    #Mosaic monochromator wave acceptance
def AccMonoMosaWav(DeltaLambda, Lambda, SigmaYp, eta, ThetaB):
    return sqrt(6/pi) * exp( - (DeltaLambda/Lambda)**2 / 2 /((SigmaYp**2+eta**2)*cotan(ThetaB))) # todo : check Sigma Yp

    #Curved monochromator angle acceptance
def AccMonoCurvAng(RMono, RInt, Wd,yp, y, r, ThetaB, Alpha, DeltaLambda, Lambda ):
    return RMono*RInt*sqrt(6/pi)/Wd * exp( - (yp-y/r/np.sin(ThetaB+Alpha)-DeltaLambda/Lambda*np.tan(ThetaB))**2/(2*Wd**2/12)) # todo check yp-y

    #Curved monochromator wave acceptance
def AccMonoCurvWav(DeltaLambda,Lambda, ThetaB, SigmaSource, DistanceFromSource, Wd):
    return sqrt(6/pi) * exp( -(DeltaLambda/Lambda)**2 / (2*cotan(ThetaB)**2*((SigmaSource/DistanceFromSource)**2+Wd**2/12))   ) #todo check too

    #Compound refractive lense with parabolic holes acceptance
def AccLensPara(x, Radius, Distance, Sigma):
    return exp( -(x**2+Radius*Distance)/2/Sigma**2)

    #Compound refractive lense with circular holes acceptance
def AccLensCirc(x, Sigma, FWHM, Radius, Distance):
    return exp( -x**2/2/Sigma**2 -x**2*FWHM**2/8/Radius**2/Sigma**2 -x**2*FWHM**4/16/Sigma**2/Radius**4 -Radius*Distance/2/Sigma**2  ) #todo why FWHM

    #Multilayer angle acceptance
def AccMultAng(Rml, yp, y, Rowland, ThetaML, DeltaLambda, Lambda, Ws):
    return Rml*8/3*sqrt(log(2)/pi) * exp( -(-yp-y/Rowland/np.sin(ThetaML)-DeltaLambda/Lambda*np.tan(ThetaML))**2/2/Ws**2 )

    #Multilayer wave acceptance
def AccMultWav(DeltaLambda,Lambda, SigmaSource, DistanceFromSource, Ws, ThetaML):
    return sqrt(6/pi) * exp( -(DeltaLambda/Lambda)**2/2/((SigmaSource/DistanceFromSource)**2+Ws**2/8/log(2))/cotan(ThetaML)**2)

#Testing the functions

def testMatFlatMono():
    print(">> MatFlatMono in test")
    b=0.8
    c=0.25
    ThetaB=1.0
    Mat=np.array([[0.2,0,0],[0,1.25,-0.389352],[0,0,1]])
    print(Mat.shape)
    result = MatFlatMono(b,c,ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatFlatMono())

def testMatCurvMono():
    print(">> MatCurvMono in test")
    b = 0.8
    Fc = 2
    c = 0.25
    ThetaB = 1.0
    Mat = np.array([[0.2, 0, 0], [0.5, 1.25, -0.389352], [0, 0, 1]])
    print(Mat.shape)
    result = MatCurvMono(b, Fc, ThetaB, c)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))
print(testMatCurvMono())

def testMatMosMono():
    print(">> MatMosMono in test")
    c = 0.25
    ThetaB = 1.0
    Mat = np.array([[0.25, 0, 0], [0, -1, 2*np.tan(1)], [0, 0, 1]])
    print(Mat.shape)
    result = MatMosMono(ThetaB, c)
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