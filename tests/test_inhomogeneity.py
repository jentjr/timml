from timml import *

def test_counterclockwise_xy():
    ccxy = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
    ml =  ModelMaq(z=[1, 0, -1], c=[100], topboundary="semi", hstar=1.0)
    inhom = PolygonInhomMaq(ml, ccxy)
    w = Well(ml)
    ml.solve()
    
    ml = ModelMaq(kaq=[1], z=[10, 0])
    uf = Uflow(ml, slope=0.01, angle=0)
    ld = PolygonInhomMaq(ml, xy=ccxy, c=1, order=3)
    ml.solve()

def test_clockwise_xy():
    cwxy = [(-1, -1), (-1, 1), (1, 1), (1, -1)]
    ml =  ModelMaq(z=[1, 0, -1], c=[100], topboundary="semi", hstar=1.0)
    inhom = PolygonInhomMaq(ml, xy=cwxy)
    w = Well(ml)
    ml.solve()
    
    ml = ModelMaq(kaq=[1], z=[10, 0])
    uf = Uflow(ml, slope=0.01, angle=0)
    ld = PolygonInhomMaq(ml, xy=cwxy, c=1, order=3)
    ml.solve()
