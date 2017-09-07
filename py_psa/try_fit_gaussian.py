
import numpy as np
import scipy


# https://gist.github.com/andrewgiessel/6122739
def gaussian(height, center_x, center_y, width_x, width_y, rotation):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)

    rotation = np.deg2rad(rotation)
    center_x = center_x * np.cos(rotation) - center_y * np.sin(rotation)
    center_y = center_x * np.sin(rotation) + center_y * np.cos(rotation)


    def rotgauss(x,y):
        xp = x * np.cos(rotation) - y * np.sin(rotation)
        yp = x * np.sin(rotation) + y * np.cos(rotation)
        g = height*np.exp(
            -(((center_x-xp)/width_x)**2+
              ((center_y-yp)/width_y)**2)/2.)
        return g
    return rotgauss

def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X*data).sum()/total
    y = (Y*data).sum()/total
    col = data[:, int(y)]
    width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
    height = data.max()
    return height, x, y, width_x, width_y, 0.0


def fitgaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = scipy.optimize.leastsq(errorfunction, params)
    return p

def IXXP(x,xp,z=0):  # xp in microradiants?????
    return np.exp(-(-5.0*x + 3.63797880709171e-11*xp)**2/2 - (5.0*x - 66666.6666666667*xp)**2/2)

def get_sigmas(function=IXXP, limit1=100, limit2=100e-6, npoints=500):


    x = np.linspace(-limit1, limit1, npoints)
    xp = np.linspace(-limit2, limit2, npoints)
    X = np.outer(x, np.ones_like(xp))
    XXP = np.zeros_like(X)

    for ix, xval in enumerate(x):
        for ixp, xpval in enumerate(xp):
            XXP[ix, ixp] = IXXP(xval, xpval)

    # plot_image(XXP,x,xp,aspect='auto',title="data")

    p = fitgaussian(XXP)

    return [p[3],p[4]]



if __name__ == "__main__":

    from psa_functions_numeric import *
    from srxraylib.plot.gol import plot_image

    flag_initial_data = 0



    #
    # initial data
    #

    if flag_initial_data == 0:

        IotaX = 3.4359738368
        IotaXp = 0.0004194304
        NumPoints = 200

        plotAB(IXXP, IotaX, IotaXp, 500, title="XXP", xtitle="X", ytitle="XP")


        x = np.linspace(-IotaX, IotaX, NumPoints)
        xp = np.linspace(-IotaXp, IotaXp, NumPoints)
        X = np.outer(x, np.ones_like(xp))
        XXP = np.zeros_like(X)

        for ix, xval in enumerate(x):
            for ixp, xpval in enumerate(xp):
                XXP[ix, ixp] = IXXP(xval, xpval)

        plot_image(XXP,x,xp,aspect='auto',title="data")

        p = fitgaussian(XXP)
        print("moments h x y wx wy 0 in pixels: ",moments(XXP))
        print("fit h x y wx wy 0 in pixels: ",p)

        print("fitted center: ",x[0]+p[2]*(x[1]-x[0]),xp[0]+p[3]*(xp[1]-xp[0]))
        print("fitted width: ",p[4]*(x[1]-x[0]),p[5]*(xp[1]-xp[0]))


        pp = get_sigmas(IXXP,IotaX,IotaXp,NumPoints)
        print("******Integral from fit: ",1e-6 * 2*np.pi*pp[0]*pp[1]*(x[1]-x[0])*(xp[1]-xp[0]))

    else:
        NumPoints = 400

        IotaX = 100
        IotaXp = 800e-6

        # plotAB(IXXP, IotaX, IotaXp, 500, title="XXP", xtitle="X", ytitle="XP")

        x = np.linspace(-IotaX, IotaX, NumPoints)
        xp = np.linspace(-IotaXp, IotaXp, NumPoints)
        X = np.outer(x, np.ones_like(xp))
        XXP = np.zeros_like(X)

        pixelX = 2*IotaX/NumPoints
        pixelY = 2*IotaXp/NumPoints

        p0 = [1.0, 0.5*NumPoints, 0.5*NumPoints, IotaX/10/pixelX, IotaXp/25/pixelY, 0]
        pfunc = gaussian(p0[0],p0[1],p0[2],p0[3],p0[4],p0[5])
        for ix, xval in enumerate(x):
            for ixp, xpval in enumerate(xp):
                XXP[ix, ixp] = pfunc(ix, ixp)

        plot_image(XXP,x,xp,aspect='auto',title="data")

        p = fitgaussian(XXP)
        print("moments h x y wx wy 0 in pixels: ",moments(XXP))
        print("fit h x y wx wy 0 in pixels: ",p)

        print("fitted center: ",x[0]+p[2]*(x[1]-x[0]),xp[0]+p[3]*(xp[1]-xp[0]))
        print("fitted width: ",p[4]*(x[1]-x[0]),p[5]*(xp[1]-xp[0]))

        print("Integral source: ",2*np.pi*p0[3]*p0[4]*(x[1]-x[0])*(xp[1]-xp[0]))

    #
    # fit data
    #

    # FIT = np.zeros_like(X)
    # fittedgaussian = gaussian(p[0],p[1],p[2],p[3],p[4],p[5])
    # for ix, xval in enumerate(x):
    #     for ixp, xpval in enumerate(xp):
    #         FIT[ix, ixp] = fittedgaussian(ix, ixp)
    #
    # plot_image(FIT,aspect='auto',title="Fitted")
    #
    # print("X step: ",x[1]-x[0],"XP step: ",xp[1]-xp[0])
    # print("Integral fit: ", 2*np.pi*p[3]*p[4]*(x[1]-x[0])*(xp[1]-xp[0]))

