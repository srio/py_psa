import numpy as np
from numpy import exp
from numpy import pi
from numpy import sqrt
from numpy import log
import numpy.testing as nut
import scipy.integrate as si
import mpmath as mp
import sympy as sp
import matplotlib.pylab as plt

# This file is only there to be able to see the source functions IXXP and IYYP once the propagation is over

x, xp, y, yp, dl = sp.symbols('x xp y yp dl')

#DEFINITION OF MATHEMATICAL OBJECTS
#Definition of the transformation matrices

def matrixFlightSymbolic(L):                                                       # For a flight path of length
    return np.array([[1,-L,0],[0,1,0],[0,0,1]])

def matrixMonoPlaneSymbolic(b, ThetaB):                                         # For a perfect flat crystal monochromator
    return np.array([[b,0,0],[0,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def matrixMonoBentSymbolic(b, Fc, ThetaB):                                     # For a perfect curved crystal monochromator (meridionally and sagitally focusing)
    return np.array([[b,0,0],[1/Fc,1/b,(1-1/b)*np.tan(ThetaB)],[0,0,1]])

def matrixMonoMosaicSymbolic(ThetaB):                                             # For a mosaic monochromator
    return np.array([[1,0,0],[0,-1,2*np.tan(ThetaB)],[0,0,1]])

def matrixMirrorPlaneSymbolic(IncAng, Sigma, Lambda, Delta):                      # For the plane mirror
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[Delta,1,0],[0,0,1]])

def matrixMirrorSymbolic(IncAng, Sigma, Lambda, Delta, Fm, S):                      # For the bent and toroidal mirrors
    return exp(-(4*pi*np.sin(IncAng)*Sigma/Lambda)**2)*np.array([[1,0,0],[(1+S*Fm*Delta)/Fm,1,0],[0,0,1]])

def matrixCompLensParaSymbolic(F):                                                 # For the compound refractive lenses with parabolic holes
    return np.array([[1,0,0],[1/F,1,0],[0,0,1]])

def matrixCompLensCircSymbolic(Coef, F):                                           # For the compound refractive lenses with circular holes
    return np.array([[1,0,0],[Coef/F,1,0],[0,0,1]])

def matrixMultilayerSymbolic(Fmu):                                                  # For the multilayer
    return np.array([[1,0,0],[1/Fmu,1,0],[0,0,1]])

#Definition of the beam source

def sourceXXPSymbolic(x, xp, SigmaXSource, SigmaXPSource):
    return sp.exp(-( (x/SigmaXSource)**2 + (xp/SigmaXPSource)**2) / 2 )
def sourceYYPSymbolic(y, yp, SigmaYSource, SigmaYPSource, GammaSource):
    return sp.exp( -( (y/SigmaYSource)**2 + ((yp-GammaSource*y)/SigmaYPSource)**2)/2 )
def sourceLambdaSymbolic(dl, SigmaSLambda):
    return sp.exp(-dl**2/2/SigmaSLambda**2)

#Definition of the acceptance windows of the optical parts

def cotanSymbolic(x):
    return 1/np.tan(x)

def acceptanceSlitSymbolic(y, Aperture, calctype):                                                           #Slit acceptance
    if calctype==0:
        return sqrt(6/pi) / sqrt(6*log(2)/pi) * sp.exp( -(y/Aperture)**2/2*12)
    if calctype==1:
        return 1/sqrt(6*log(2)/pi)*sp.exp(-y**2/(2*Aperture**2/2/pi))
    else:
        return 'the value for calctype is 0 or 1 you pipsqueak !'

def acceptancePinSymbolic(y, Diameter):                                                            #Pinhole acceptance
    return sqrt(8/pi) * sp.exp ( -(y/Diameter)**2/2*16 )

def acceptanceAngleMonoPlaneSymbolic(yp, DeltaLambda, ThetaB, Wd, RMono, RInt):                #Plane monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * sp.exp( -(yp-DeltaLambda*np.tan(ThetaB))**2 / (2*Wd**2/12))

def acceptanceWaveMonoPlaneSymbolic(DeltaLambda, SigmaYp, ThetaB, Wd):                        #Plane monochromator wave acceptance
    return sqrt(6/pi) * sp.exp( - (DeltaLambda)**2 / (2*(SigmaYp**2+Wd**2/12)*cotanSymbolic(ThetaB)**2) )

def acceptanceAngleMonoMosaicSymbolic(yp, DeltaLambda, ThetaB, eta, RInt):                       #Mosaic monochromator angular acceptance
    return RInt*sqrt(6/pi)/eta * sp.exp(- (yp-DeltaLambda*np.tan(ThetaB))**2 / eta**2 /2 )

def acceptanceWaveMonoMosaicSymbolic(DeltaLambda, SigmaYp, ThetaB, eta):                        #Mosaic monochromator wave acceptance
    return sqrt(6/pi) * sp.exp( - (DeltaLambda)**2 / 2 /((SigmaYp**2+eta**2)*cotanSymbolic(ThetaB)**2))

def acceptanceAngleMonoBentSymbolic(y, yp, DeltaLambda, Alpha, ThetaB, Wd, r, RMono, RInt):    #Curved monochromator angle acceptance
    return RMono*RInt*sqrt(6/pi)/Wd * sp.exp( - (yp-y/r/np.sin(ThetaB+Alpha)-DeltaLambda*np.tan(ThetaB))**2/(2*Wd**2/12))

def acceptanceWaveMonoBentSymbolic(DeltaLambda, ThetaB, DistanceFromSource, Wd, SigmaSource):  #Curved monochromator wave acceptance
    return sqrt(6/pi) * sp.exp( -(DeltaLambda)**2 / (2*cotanSymbolic(ThetaB)**2*((SigmaSource/DistanceFromSource)**2+Wd**2/12))   )

def acceptanceCompLensParaSymbolic(x, Radius, Distance, Sigma):                                          #Compound refractive lense with parabolic holes acceptance
    return sp.exp( -(x**2+Radius*Distance)/2/Sigma**2)

def acceptanceCompLensCircSymbolic(x, Radius, Distance, Sigma, FWHM):                                    #Compound refractive lense with circular holes acceptance
    return sp.exp( -x**2/2/Sigma**2 -x**2*FWHM**2/8/Radius**2/Sigma**2 -x**2*FWHM**4/16/Sigma**2/Radius**4 -Radius*Distance/2/Sigma**2  )

def acceptanceAngleMultiSymbolic(y, yp, DeltaLambda, ThetaML, Rml, Rowland, Ws):                #Multilayer angle acceptance
    return Rml*8/3*sqrt(log(2)/pi) * sp.exp( -(-yp-y/Rowland/np.sin(ThetaML)-DeltaLambda*np.tan(ThetaML))**2*8*log(2)/2/Ws**2 )

def acceptanceWaveMultiSymbolic(DeltaLambda, ThetaML, DistanceFromSource, Ws, SigmaSource):     #Multilayer wave acceptance
    return sqrt(6/pi) * sp.exp( -(DeltaLambda)**2/2/((SigmaSource/DistanceFromSource)**2+Ws**2/8/log(2))/cotanSymbolic(ThetaML)**2)


# Useful functions for the calculation part

def buildMatTabSymbolic(ListObject, ListDistance):
    n = len(ListObject)
    if ListDistance == []:
        return 'Error, ListDistance is empty'
    else :
        MatTabX = [ListDistance[0]]
        MatTabY = [ListDistance[0]]
        if n == 0:
            return [MatTabX, MatTabY]
        else :
            for k in range(n):
                    MatTabX.append(ListObject[k][-1])
                    MatTabX.append(matrixFlightSymbolic(ListDistance[k+1]))
                    MatTabY.append(ListObject[k][-2])
                    MatTabY.append(matrixFlightSymbolic(ListDistance[k+1]))
            return [MatTabX, MatTabY]

def propagateMatrixListSymbolic(x, xp, y, yp, dl, SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, bMonoX, bMonoY):
    # initiating variables, copies of the arrays are created
    MX = MatTabX.copy()
    MY = MatTabY.copy()
    MatTempX = np.array([x, xp, dl])
    MatTempY = np.array([y, yp, dl])
    # the case of single propagation has to be done separately because of a pb of range
    if len(MX)==1:
        MatTempX = np.dot(MX[0], MatTempX)
        MatTempY = np.dot(MY[0], MatTempY)
        NewSourceX = sourceXXPSymbolic(MatTempX[0], MatTempX[1], SigmaXSource, SigmaXPSource)
        NewSourceY = sourceYYPSymbolic(MatTempY[0], MatTempY[1], SigmaYSource, SigmaYPSource, GammaSource)
    # now we do the general case
    else :
        #first matrix multiplication
        for i in range(len(MX)-1, -1, -1):
            MatTempX = np.dot(MX[i], MatTempX)
            MatTempY = np.dot(MY[i], MatTempY)
        NewSourceX = sourceXXPSymbolic(MatTempX[0], MatTempX[1], SigmaXSource, SigmaXPSource)
        NewSourceY = sourceYYPSymbolic(MatTempY[0], MatTempY[1], SigmaYSource, SigmaYPSource, GammaSource)
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
            if ListObject[k][0] == 'Slit':
                NewSourceX = NewSourceX * acceptanceSlitSymbolic(MatTempX[0], ListObject[k][1], ListObject[k][3])
                NewSourceY = NewSourceY * acceptanceSlitSymbolic(MatTempY[0], ListObject[k][2], ListObject[k][3])
            elif ListObject[k][0] == 'Pinhole' :
                NewSourceX = NewSourceX * acceptancePinSymbolic(MatTempX[0], ListObject[k][1])
                NewSourceY = NewSourceY * acceptancePinSymbolic(MatTempX[0], ListObject[k][1])
            elif ListObject[k][0] == 'MonoPlaneVertical':
                NewSourceY = NewSourceY * acceptanceAngleMonoPlaneSymbolic(MatTempY[1], MatTempY[2], ListObject[k][1], ListObject[k][2], ListObject[k][3], ListObject[k][4])
                NewSourceY = NewSourceY * acceptanceWaveMonoPlaneSymbolic(MatTempY[2], bMonoY*SigmaYPSource, ListObject[k][1], ListObject[k][2])
            elif ListObject[k][0] == 'MonoPlaneHorizontal':
                NewSourceX = NewSourceX * acceptanceAngleMonoPlaneSymbolic(MatTempX[1], MatTempX[2], ListObject[k][1], ListObject[k][2], ListObject[k][3], ListObject[k][4])
                NewSourceX = NewSourceX * acceptanceWaveMonoPlaneSymbolic(MatTempX[2], bMonoX*SigmaXPSource, ListObject[k][1], ListObject[k][2])
            elif ListObject[k][0] == 'MultiHorizontal':
                NewSourceX = NewSourceX * acceptanceAngleMultiSymbolic(MatTempX[0], MatTempX[1], MatTempX[2], ListObject[k][1], ListObject[k][2], ListObject[k][3], ListObject[k][4])
                NewSourceX = NewSourceX * acceptanceWaveMultiSymbolic(MatTempX[2], ListObject[k][1], ListObject[k][5], ListObject[k][4], SigmaXSource)
            elif ListObject[k][0] == 'MultiVertical' :
                NewSourceY = NewSourceY * acceptanceAngleMultiSymbolic(MatTempY[0], MatTempY[1], MatTempY[2], ListObject[k][1], ListObject[k][2], ListObject[k][3], ListObject[k][4])
                NewSourceY = NewSourceY * acceptanceWaveMultiSymbolic(MatTempY[2], ListObject[k][1], ListObject[k][5], ListObject[k][4], SigmaYSource)
            elif ListObject[k][0] == 'MonoBentHorizontal':
                NewSourceX = NewSourceX * acceptanceAngleMonoBentSymbolic(MatTempX[0], MatTempX[1], MatTempX[2], ListObject[k][1], ListObject[k][2], ListObject[k][3], ListObject[k][4], ListObject[k][5], ListObject[k][6])
                NewSourceX = NewSourceX * acceptanceWaveMonoBentSymbolic(MatTempX[2], ListObject[k][2], ListObject[k][7], ListObject[k][3], SigmaXSource)
            elif ListObject[k][0] == 'MonoBentVertical' :
                NewSourceY = NewSourceY * acceptanceAngleMonoBentSymbolic(MatTempY[0], MatTempY[1], MatTempY[2],ListObject[k][1], ListObject[k][2], ListObject[k][3], ListObject[k][4], ListObject[k][5], ListObject[k][6])
                NewSourceY = NewSourceY * acceptanceWaveMonoBentSymbolic(MatTempY[2], ListObject[k][2], ListObject[k][7], ListObject[k][3], SigmaXSource)
            elif ListObject[k][0] == 'MonoMosaicHorizontal':
                NewSourceX = NewSourceX * acceptanceAngleMonoMosaicSymbolic(MatTempX[1], MatTempX[2], ListObject[k][1], ListObject[k][2], ListObject[k][3])
                NewSourceX = NewSourceX * acceptanceWaveMonoMosaicSymbolic(MatTempX[2], SigmaXPSource, ListObject[k][1], ListObject[k][2])
            elif ListObject[k][0] == 'MonoMosaicVertical' :
                NewSourceY = NewSourceY * acceptanceAngleMonoMosaicSymbolic(MatTempY[1], MatTempY[2], ListObject[k][1], ListObject[k][2], ListObject[k][3])
                NewSourceY = NewSourceY * acceptanceWaveMonoMosaicSymbolic(MatTempY[2], SigmaYPSource, ListObject[k][1], ListObject[k][2])
            elif ListObject[k][0] == 'LensParabolicHorizontal':
                NewSourceX = NewSourceX * acceptanceCompLensParaSymbolic(MatTempX[0], ListObject[k][1], ListObject[k][2],
                                                                 ListObject[k][3])
            elif ListObject[k][0] == 'LensParabolicVertical':
                NewSourceY = NewSourceY * acceptanceCompLensParaSymbolic(MatTempY[0], ListObject[k][1], ListObject[k][2],
                                                                 ListObject[k][3])
            elif ListObject[k][0] == 'LensParabolic2D':
                NewSourceX = NewSourceX * acceptanceCompLensParaSymbolic(MatTempX[0], ListObject[k][1], ListObject[k][2],
                                                                 ListObject[k][3])
                NewSourceY = NewSourceY * acceptanceCompLensParaSymbolic(MatTempY[0], ListObject[k][1], ListObject[k][2],
                                                                 ListObject[k][3])
            k = k +1
            del MX[0:2]
            del MY[0:2]
    return [NewSourceX, NewSourceY]

def sourceFinaleSymbolic(SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, SigmaSLambda, GammaSource, MatTabX, MatTabY, ListObject, SourceI, bMonoX, bMonoY):
    # IXXP = lambda x, xp, dl : propagateMatrixList(x, xp, 0, 0, dl, SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, bMonoX, bMonoY)[0]
    # IYYP = lambda y, yp, dl : propagateMatrixList(0, 0, y, yp, dl, SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, bMonoX, bMonoY)[1]
    # ISigma = lambda dl : sourceLambda(dl, SigmaSLambda) * SourceI
    IXXP = propagateMatrixListSymbolic(x, xp, 0, 0, dl, SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, bMonoX, bMonoY)[0]
    IYYP = propagateMatrixListSymbolic(0, 0, y, yp, dl, SigmaXSource, SigmaXPSource, SigmaYSource, SigmaYPSource, GammaSource, MatTabX, MatTabY, ListObject, bMonoX, bMonoY)[1]
    ISigma = sourceLambdaSymbolic(dl, SigmaSLambda) * SourceI
    return IXXP,IYYP, ISigma


