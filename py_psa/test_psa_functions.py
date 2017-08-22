from psa_functions_numeric import *

#TESTING THE FUNCTIONS
#Testing the matrices

def testMatrixMonoPlane():
    print(">> MatFlatMono in test")
    b=0.8
    ThetaB=1.0
    Mat=np.array([[0.8,0,0],[0,1.25,-0.389352],[0,0,1]])
    print(Mat.shape)
    result = matrixMonoPlane(b,ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

def testMatrixMonoBent():
    print(">> MatCurvMono in test")
    b = 0.8
    Fc = 2
    ThetaB = 1.0
    Mat = np.array([[0.8, 0, 0], [0.5, 1.25, -0.389352], [0, 0, 1]])
    print(Mat.shape)
    result = matrixMonoBent(b, Fc, ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

def testMatrixMonoMosaic():
    print(">> MatMosMono in test")
    ThetaB = 1.0
    Mat = np.array([[1, 0, 0], [0, -1, 2*np.tan(1)], [0, 0, 1]])
    print(Mat.shape)
    result = matrixMonoMosaic(ThetaB)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

def testMatrixMirrorPlane():
    print(">> MatFlatMirr in test")
    Sigma = 0.02
    IncAng = 30
    Lambda = 1
    Delta = 5
    Mat = np.array([[0.9402, 0, 0], [4.701, 0.9402, 0], [0, 0, 0.9402]])
    print(Mat.shape)
    result = matrixMirrorPlane(IncAng, Sigma, Lambda, Delta)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

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
    result = matrixMirror(IncAng, Sigma, Lambda, Delta, Fm, S)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

def testMatrixCompLensPara():
    print(">> MatParaLens in test")
    F=12.3
    Mat = np.array([[1, 0, 0], [0.0813008, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = matrixCompLensPara(F)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

def testMatrixCompLensCirc():
    print(">> MatCircLens in test")
    F = 12.3
    Coef = 0.4
    Mat = np.array([[1, 0, 0], [0.0325203, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = matrixCompLensCirc(Coef, F)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

def testMatrixMultilayer():
    print(">> MatMulti in test")
    Fmu = 15.1
    Mat = np.array([[1, 0, 0], [0.0662252, 1, 0], [0, 0, 1]])
    print(Mat.shape)
    result = matrixMultilayer(Fmu)
    print("Mistakes ?", nut.assert_array_almost_equal(Mat, result, decimal=5))

#Testing the source

def testSourceXXP():
    print(">> SourceX in test")
    x = 0.5
    xp = 2
    Source = 0.605773
    result = sourceXXP(x, xp, SigmaXSource=10, SigmaXPSource=2)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))

def testSourceYYP():
    print(">> SourceY in test")
    y = 0.5
    yp = 2
    Source = 0.753897
    result = sourceYYP(y, yp, SigmaYSource=10, SigmaYPSource=2, GammaSource=1)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))

def testSourceLambda():
    print(">> SourceLambda in test")
    DeltaLambda = 4
    Source = 0.278037
    result = sourceLambda(DeltaLambda, SigmaSLambda=2.5)
    print(Source, result)
    print("Mistakes ?", nut.assert_almost_equal(Source, result, decimal=5))

#Testing the acceptances

def testAcceptanceSlit():
    print(">> AccSlit in test")
    y = 2
    Aperture = 75
    Acceptance = 1.19601
    result = acceptanceSlit(y, Aperture, 0)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))
    result2 = acceptanceSlit(y, Aperture, 1)
    Acceptance2=0.867194
    print(Acceptance2, result2)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance2, result2, decimal=5))

def testAcceptancePin():
    print(">> AccPinh in test")
    y = 2
    Diameter = 75
    Acceptance = 1.58672
    result = acceptancePin(y, Diameter)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceAngleMonoPlane():
    print(">> AccMonoPlanAng in test")
    RMono = 0.91
    RInt = 0.2
    Wd = 103.2
    yp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 0.00214474
    result = acceptanceAngleMonoPlane(yp, DeltaLambda, ThetaB, Wd, RMono, RInt)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceWaveMonoPlane():
    print(">> AccMonoPlanWav in test")
    Wd = 103.2
    SigmaYp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.21616
    result = acceptanceWaveMonoPlane(DeltaLambda, SigmaYp, ThetaB, Wd)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceAngleMonoMosaic():
    print(">> AccMonoMosaAng in test")
    RMono = 0.91
    RInt = 0.2
    eta = 103.2
    yp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 0.00264987
    result = acceptanceAngleMonoMosaic(yp, DeltaLambda,ThetaB, eta, RInt)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceWaveMonoMosaic():
    print(">> AccMonoMosaWav in test")
    eta = 103.2
    SigmaYp = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.36733
    result = acceptanceWaveMonoMosaic(DeltaLambda, SigmaYp, ThetaB, eta)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

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
    result = acceptanceAngleMonoBent(y, yp, DeltaLambda, Alpha, ThetaB, Wd, r, RMono, RInt)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceWaveMonoBent():
    print(">> AccMonoCurvWav in test")
    Wd = 103.2
    SigmaSource = 5.5
    DistanceFromSource = 0.001
    DeltaLambda = 0.2
    ThetaB = 33
    Acceptance = 1.38197
    result = acceptanceWaveMonoBent(DeltaLambda, ThetaB, DistanceFromSource, Wd, SigmaSource)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceCompLensPara():
    print(">> AccLensPara in test")
    x = 2
    Radius = 75
    Distance = 22.2
    Sigma = 33.3
    Acceptance = 0.471162
    result = acceptanceCompLensPara(x, Radius, Distance, Sigma)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceCompLensCirc():
    print(">> AccLensCirc in test")
    x = 2
    Radius = 75
    Distance = 22.2
    FWHM = 43.3
    Sigma = 33.3
    Acceptance = 0.471079
    result = acceptanceCompLensCirc(x, Radius, Distance, Sigma, FWHM)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

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
    result = acceptanceAngleMulti(y, yp, DeltaLambda, ThetaML, Rml, Rowland, Ws)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

def testAcceptanceWaveMulti():
    print(">> AccMultWav in test")
    ThetaML = 1
    DeltaLambda = 1
    Ws = 0.7
    DistanceFromSource = 22
    SigmaSource = 45
    Acceptance = 1.04044
    result = acceptanceWaveMulti(DeltaLambda, ThetaML, DistanceFromSource, Ws, SigmaSource)
    print(Acceptance, result)
    print("Mistakes ?", nut.assert_almost_equal(Acceptance, result, decimal=5))

print(testMatrixMonoPlane())
print(testMatrixMonoBent())
print(testMatrixMonoMosaic())
print(testMatrixMirrorPlane())
print(testMatrixMirror())
print(testMatrixCompLensPara())
print(testMatrixCompLensCirc())
print(testMatrixMultilayer())

print(testSourceXXP())
print(testSourceYYP())
print(testSourceLambda())

print(testAcceptanceSlit())
print(testAcceptancePin())
print(testAcceptanceAngleMonoPlane())
print(testAcceptanceWaveMonoPlane())
print(testAcceptanceAngleMonoMosaic())
print(testAcceptanceWaveMonoMosaic())
print(testAcceptanceAngleMonoBent())
print(testAcceptanceWaveMonoBent())
print(testAcceptanceCompLensPara())
print(testAcceptanceCompLensCirc())
print(testAcceptanceAngleMulti())
print(testAcceptanceWaveMulti())