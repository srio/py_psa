# This file is to convert results from Mathematica, given as FWHM, to Sigmas
import numpy as np
x = 2.48375
xp = 70.644600
y = 0.00848213
yp = 6.0451100
dl = 2.322650e-3
Tab = np.array([ x ,    xp  ,    y ,    yp ,     dl  ])
print(Tab)
print(1/2/np.sqrt(2*np.log(2)) * Tab)