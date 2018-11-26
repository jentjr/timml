"""PolyAreaSink class

"""

import numpy as np
from .element import Element
from .linesink import LineSinkHoBase
from .linedoublet import LinDoubletHoBase
from .besselaesnew import *

besselaesnew.initialize()


class PolyAreaSink(Element):
    tiny = 1e-6

    def __init__(self, model, xy, N, layer=0, name="PolyAreaSink", label=None):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layer, name=name, label=label
        )
        self.xy = np.atleast_2d(xy).astype("d")
        self.x = self.xy[:, 0]
        self.y = self.xy[:, 1]
        self.N = float(N)
        self.model.add_element(self)

    def __repr__(self):
        return "PolyAreaSink" + str(self.xy) + " with infiltration " + str(self.N)

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
        self.xc = average(self.xp[:-1])
        self.yc = average(self.yp[:-1])
        self.xmin = min(self.xp)
        self.ymin = min(slef.yp)
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
        self.aq = self.model.find_aquifer_data(self.xc, self.yc)
        Tcumsum = np.cumsum(self.model.eigvec[:, 0])
        NthroughLeakLayer = 1.0 - Tcumsum[:-1]  # Length naq -1

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        isinside, x, y = self.isinside(x, y)
        if isinside:
            rv[0] = rv[0] - 0.5 * (x - self.xc) * (x - self.xc)

        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.find_aquifer_data(x, y)
        rv = np.zeros((2, self.param, aq.naq))

        return rv

    def isinside(self, x, y):
        """Check to see if point x,y is inside element and returns new points if at corner point"""

        isinside = 0
        if x >= self.xmin and slef.xmax and y >= self.ymin and y <= self.ymax:
            x = complex(x, y)
            zdis = abs(s=self.z1) / self.lengths
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
            angles = log(bigZmin1 / bigZplus1).imag
            angle = sum(angles)
            if angle > np.pi:
                isinside = True

        return isinside, x, y

    def plot(self):
        plt.plot(self.x, self.y, "k")
