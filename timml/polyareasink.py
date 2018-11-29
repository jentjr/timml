"""PolyAreaSink class

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from .element import Element
from .linesink import LineSinkHoBase
from .linedoublet import LineDoubletHoBase
from .besselaesnew import *

besselaesnew.initialize()


class PolyAreaSink(Element):
    """Uniform infiltration over a polygon area
    
    Parameters
    ----------
    model: Model object
        Model to which the element is added
    xy : array or list
        list or array of (x,y) pairs of coordinates of polygon [L]
    N:  numeric
        Infiltration rate [L/T]
    layer: integer, list, or array
        layer to which the element is added
    name: str, or None
        name of the element
    label: str, or None (default: None)
        label of the element
    
    """
    
    tiny = 1e-6

    def __init__(self, model, xy, N, layer=0, name="PolyAreaSink", label=None):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layer, name=name, label=label
        )
        self.xylist = xy
        self.xy = np.atleast_2d(xy).astype("d")
        self.x = self.xy[:, 0]
        self.y = self.xy[:, 1]
        self.N = float(N)
        self.model.add_element(self)

    def __repr__(self):
        return "PolyAreaSink" + str(self.xylist) + " with infiltration " + str(self.N)

    def initialize(self):
        self.parameters = np.array([[self.N]])
        self.Nsides = len(self.xy)
        self.xp = np.zeros(self.Nsides + 1, "d")
        self.yp = np.zeros(self.Nsides + 1, "d")
        self.z1 = np.zeros(self.Nsides, "D")
        self.z2 = np.zeros(self.Nsides, "D")
        for i in range(self.Nsides):
            self.xp[i] = self.xy[i][0]
            self.yp[i] = self.xy[i][1]
            self.z1[i] = complex(self.xp[i], self.yp[i])
        self.xp[-1] = self.xy[0][0]
        self.yp[-1] = self.xy[0][1]
        self.z2[:-1] = self.z1[1:]
        self.z2[-1] = self.z1[0]
        self.lengths = abs(self.z2 - self.z1)
        
        # find center of area-sink
        self.xc = np.average(self.xp[:-1])
        self.yc = np.average(self.yp[:-1])
        self.xmin = min(self.xp)
        self.ymin = min(self.yp)
        self.xmax = max(self.xp) + 2.0 * self.tiny
        self.ymax = max(self.yp) + 2.0 * self.tiny
        self.xpmin = np.zeros(self.Nsides, "d")
        self.ypmin = np.zeros(self.Nsides, "d")
        self.xpmax = np.zeros(self.Nsides, "d")
        self.ypmax = np.zeros(self.Nsides, "d")
        for i in range(self.Nsides):
            self.xpmin[i] = min(self.xp[i], self.xp[i + 1])
            self.xpmax[i] = max(self.xp[i], self.xp[i + 1])
            self.ypmin[i] = min(self.yp[i], self.yp[i + 1])
            self.ypmax[i] = max(self.yp[i], self.yp[i + 1])
        
        self.area = polygon_area(self.xy)    
        
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.aq.add_element(self)
        
    def zetaj(self, x, y):
        z = x + 1j*y
        zetaj = (2*z - self.z1 - self.z2)/(self.z2 - self.z1)
        return zetaj

    def zetajplus(self, x, y):
        z = x + 1j*y
        zetaplus = 2.0*(z - self.z1)/(self.z2 - self.z1)
        return zetaplus
    
    def zetajminus(self, x, y):
        z = x + 1j*y
        zetaminus = 2.0*(z - self.z2)/(self.z2 - self.z1)
        return zetaminus
    
    def zetamplus(self, x, y):
        z = x + 1j*y
        zetamplus = 2*(z - self.z1[1:])/(self.z2[1:] - self.z1[1:])
        return zetamplus
    
    def zetamminus(self, x, y):
        z = x + 1j*y
        zetamminus = 2*(z - self.z2[1:])/(self.z2[1:] - self.z1[1:])
        return zetamminus
    
    def etaj(self, x, y):
        zjp = self.zetajplus(x, y)
        zjm = self.zetajminus(x, y)
        zmm = self.zetamminus(x, y)
        zmp = self.zetamplus(x, y)
        etaj = zjp*np.log(zjm/zjp) + 2*np.sum(np.log(zmm/zmp) + 2)
        return etaj

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        
        z = x + 1j*y
        yj = self.zetaj(x, y).conj() - self.zetaj(x, y)
        Ljsq = (self.z2 - self.z1) * (self.z2 - self.z1).conj()
        ej = self.etaj(x, y)
        
        isinside, x, y = self.isinside(x, y)
        
        if aq == self.model.aq:
            rv[0] = -self.N/(32*np.pi)*np.sum(yj*Ljsq*(ej + ej.conj()) + self.N*self.area/(4*np.pi)*(np.log(z - self.z1) + np.log(z.conjugate() - self.z1.conj())))
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.find_aquifer_data(x, y)
        disx = np.zeros(aq.naq)
        disy = np.zeros(aq.naq)
        isinside, x, y = self.isinside(x, y)
        for ls in self.laplacels:
            pass
        for ld in self.besselld:
            pass
        for ld in self.lapalceld:
            pass
        if isinside:
            disx[0] = disx[0] + (x - self.xc)
        return disx, disy

    def isinside(self, x, y):
        """Check to see if point x,y is inside element and returns new points if at corner point"""

        isinside = 0
        if x >= self.xmin and self.xmax and y >= self.ymin and y <= self.ymax:
            z = complex(x, y)
            zdis = abs(z-self.z1) / self.lengths
            if min(zdis) < self.tiny:
                for i in range(len(zdis)):
                    if zdis[i] < self.tiny:
                        z = z + self.tiny * self.lengths[i]
                        x = z.real
                        y = z.imag
                        break
            bigZ = (2.0 * z - (np.array(self.z1) + np.array(self.z2))) / (
                np.array(self.z2) - np.array(self.z1)
            )
            bigZmin1 = bigZ - 1.0
            bigZplus1 = bigZ + 1.0
            angles = np.log(bigZmin1 / bigZplus1).imag
            angle = np.sum(angles)
            if angle > np.pi:
                isinside = True

        return isinside, x, y

    def plot(self):
        fig = plt.gcf()
        ax = plt.gca()
        poly = Polygon(self.xy, closed=True, fill=False)
        ax.add_patch(poly)
        plt.plot()


def polygon_area(xy):
    n = len(xy)
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += xy[i][0] * xy[j][1]
        area -= xy[j][0] * xy[i][1]
    area = abs(area) / 2.0
    return area
