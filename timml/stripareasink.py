import numpy as np
import inspect  # Used for storing the input
from .element import Element

__all__ = ["StripAreaSink"]


class StripAreaSink(Element):
    """Uniform infiltration over a strip of uniform width
    
    Parameters
    ----------
    model: Model object
        Model to which the element is added
    xleft: numeric
        x-coordinate for left side of strip [L]
    xright: numeric
        x-coordinate for right side of strip [L]
    N:  numeric
        Infiltration rate [L/T]
    layer: integer, list, or array
        layer to which the element is added
    name: string, or None
        name of the element
    label: string, or None (default: None)
        label of the element
    
    """

    def __init__(
        self,
        model,
        xleft=-1,
        xright=1,
        N=0.001,
        layer=0,
        name="StripAreaSink",
        label=None,
    ):
        Element.__init__(
            self, model, nparam=1, nunknowns=0, layers=layer, name=name, label=label
        )
        self.xleft = xleft
        self.xright = xright
        self.N = N
        self.model.add_element(self)

    def __repr__(self):
        return self.name + " at " + str((self.xc, self.yc))

    def initialize(self):
        self.xc = 0.5 * (self.xleft + self.xright)
        self.yc = 0
        self.L = self.xright - self.xleft
        self.aq = self.model.aq.find_aquifer_data(self.xc, self.yc)
        self.aq.add_element(self)
        self.parameters = np.array([[self.N]])
        if self.aq.ilap:
            self.lab = self.aq.lab[1:]
            self.A = -self.aq.coef[self.layers, 1:] * self.lab ** 2 / 2
            self.B = self.A * (np.exp(-self.L / self.lab) - 1)
            self.plabsq = self.aq.coef[self.layers, 1:] * self.lab ** 2
        else:
            print("StripAreaSink cannot be added to semi-confined system")

    def potinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((self.nparam, aq.naq))
        if aq == self.aq:
            if x < self.xleft:
                rv[0, 0] = -(self.xleft - self.xc) * (x - self.xc) + self.L ** 2 / 8
                rv[0, 1:] = self.B * np.exp((x - self.xleft) / self.lab)
            elif x > self.xright:
                rv[0, 0] = -(self.xright - self.xc) * (x - self.xc) + self.L ** 2 / 8
                rv[0, 1:] = self.B * np.exp(-(x - self.xright) / self.lab)
            else:
                rv[0, 0] = -0.5 * (x ** 2 - 2 * self.xc * x + self.xc ** 2)
                rv[0, 1:] = (
                    self.A
                    * (
                        np.exp(-(x - self.xleft) / self.lab)
                        + np.exp((x - self.xright) / self.lab)
                    )
                    + self.plabsq
                )
        return rv

    def disvecinf(self, x, y, aq=None):
        if aq is None:
            aq = self.model.aq.find_aquifer_data(x, y)
        rv = np.zeros((2, self.nparam, aq.naq))
        if aq == self.aq:
            if x < self.xleft:
                rv[0, 0, 0] = self.xleft - self.xc
                rv[0, 0, 1:] = -self.B / self.lab * np.exp((x - self.xleft) / self.lab)
            elif x > self.xright:
                rv[0, 0, 0] = self.xright - self.xc
                rv[0, 0, 1:] = self.B / self.lab * np.exp(-(x - self.xright) / self.lab)
            else:
                rv[0, 0, 0] = x - self.xc
                rv[0, 0, 1:] = (
                    self.A
                    / self.lab
                    * (
                        np.exp(-(x - self.xleft) / self.lab)
                        - np.exp((x - self.xright) / self.lab)
                    )
                )
        return rv
