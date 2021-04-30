
from scipy import optimize
import numpy as np
import mpmath as mp
mp.dps = 200


class ComputeCurvature:
    def __init__(self):
        """ Initialize some variables """
        self.xc = 0  # X-coordinate of circle center
        self.yc = 0  # Y-coordinate of circle center
        self.r =  0   # Radius of the circle
        self.xx = np.array([])  # Data points
        self.yy = np.array([])  # Data points

    def calc_r(self, xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((self.xx-xc)**2 + (self.yy-yc)**2, dtype=np.longdouble)

    def f(self, c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        ri = self.calc_r(*c)
        return ri - ri.mean()

    def df(self, c):
        """ Jacobian of f_2b
        The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
        xc, yc = c
        df_dc = np.empty((len(c), self.xx.size))

        ri = self.calc_r(xc, yc)
        df_dc[0] = (xc - self.xx)/ri                   # dR/dxc
        df_dc[1] = (yc - self.yy)/ri                   # dR/dyc
        df_dc = df_dc - df_dc.mean(axis=1)[:, np.newaxis]
        return df_dc

    def fit(self, xx, yy):
        self.xx = xx
        self.yy = yy
        # print(type(xx[1]),type(yy[1]))
        center_estimate = np.array([np.mean(xx), np.mean(yy)])
        center = optimize.leastsq(self.f, center_estimate, Dfun=self.df, col_deriv=True, gtol=1e-8)[0]

        self.xc, self.yc = center
        ri = self.calc_r(*center)
        self.r = ri.mean()
        # self.r = mp.agm(ri)

        return self.r  # Return the radius of curvature

    def curvature(self, xx, yy):
        return 1/self.fit(xx,yy)


if __name__=="__main__":
    import matplotlib.pyplot as plt

    x = np.r_[15, 16, 20, 22, 30, 35]
    y = np.r_[10, 16, 28, 31, 37, 39]
    comp_curv = ComputeCurvature()
    curvature = comp_curv.fit(x, y)
    print(curvature)

    plt.plot(x,y)
    plt.plot(comp_curv.xc, comp_curv.yc, 'ro')
    plt.show()

