import numpy as np
from math import exp
from math import pi
from math import sqrt
from math import log

# todo: add help
#Definition of the transformation matrices
    # For a flight path of length
    def MatFlight(L):
        return np.array([[1,-L,0],[0,1,0],[0,0,1]])

    # For the first monochromator (?)
    # MatFirstMono=np.array([[c,0,0],[0,1,0],[0,0,1]])

    # For a perfect flat crystal monochromator
    def MatFlatMono(b, ThetaB):
        return np.array([[b,0,0],[0,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

    # For a perfect curved crystal monochromator (meridionally and sagitally focusing)
    def MatCurvMono(b, Fc, ThetaB):
        return np.array([[b,0,0],[1/Fc,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

    # For a mosaic monochromator
    def MatMosMono(ThetaB):
        return np.array([[1,0,0],[0,-1,2*np.tan(ThetaB)],[0,0,1]])

    # For the plane mirror
    def MatFlatMirr(IncAng, Sigma, Lambda, Delta):
        return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[Delta,1,0],[0,0,1]])

    # For the bent and toroidal mirrors
    def MatMirr(IncAng, Sigma, Lambda, Delta, Fm):
        return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[(1+S*Fm*Delta)/Fm,1,0],[0,0,1]])

    # For the compound refractive lenses with parabolic holes
    def MatParaLens(F):
        return np.array([[1,0,0],[1/F,1,0],[0,0,1]])

    # For the compound refractive lenses with circular holes
    def MatCircLens(Coef, F):
        return np.array([[1,0,0],[Coef/F,1,0],[0,0,1]])

    # For the multilayer
    def MatMulti(Fmu):
        return np.array([[1,0,0],[exp(1)/Fmu,1,0],[0,0,1]])

#Definition of the beam source

    def SourceX(x, SigmaSourceX, xp, SigmaSourceXp):
        return exp(-( (x/SigmaSourceX)**2 + (xp/SigmaSourceXp)**2) / 2 )
    def SourceY(y, SigmaSourceY, yp, GammaSource, SigmaSourceYp):
        return exp( -( (y/SigmaSourceY)**2 + ((yp-GammaSource*y)/SigmaSourceYp)**2)/2 )
    def SourceLambda(DLambda, Lambda, SigmaSourceLambda):
        return exp(-(DLambda/Lambda)**2/2/SigmaSourceLambda**2)

#Definition of the acceptance windows of the optical parts (problem line 44)

    #Slit acceptance
    def AccSlit(y, Aperture):
        return sqrt(6/pi) / sqrt(6*log(2)/pi) * exp( -(y/Aperture)**2/2*12)

    #Pinhole acceptance
    def AccPinh(y, Diameter):
        return sqrt(8/pi) * exp ( -(y/Diameter)**2/2*16 )

    #Plane monochromator angle acceptance
    def AccMonoPlanAng(RMono, RInt, Wd, DeltaLambda, Lambda, ThetaB):
        return RMono*RInt*sqrt(6/pi)/Wd * exp( -(yp-DeltaLambda/Lambda*np.tan(ThetaB))**2 / (2*Wd**2/12))

    #Plane monochromator wave acceptance
    def AccMonoPlanWav(DeltaLambda, Lambda, SigmaYp, ThetaB):
        return sqrt(6/pi) * exp( - (DeltaLambda/Lambda)**2 / (2*(SigmaYp**2+Wd**2/12)*np.cot(ThetaB)**2) )

    #Mosaic monochromator angular acceptance
    def AccMonoMosaAng(RInt, eta, yp, DeltaLambda, Lambda, ThetaB):
        return RInt*sqrt(6/pi)/eta * exp(- (yp-DeltaLambda/Lambda*np.tan(ThetaB))**2 / eta**2 /2 ) # /2 ? RInt ?

    #Mosaic monochromator wave acceptance
    def AccMonoMosaWav(DeltaLambda, Lambda, SigmaYp, eta, ThetaB):
        return sqrt(6/pi) * exp( - (DeltaLambda/Lambda)**2 / 2 /((SigmaYp**2+eta**2)*np.cot(ThetaB))) #Sigma Yp ???

    #Curved monochromator angle acceptance
    def AccMonoCurvAng(RMono,RInt):
        return RMono*RInt
