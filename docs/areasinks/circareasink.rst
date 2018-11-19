Circular Area-Sink
------------------
The comprehensive potential for a circular pond of infiltration :math:`N`, centered at :math:`(xc, yc)` 
and of radius :math:`R` is:
     
     .. math:: 
        
        r = \sqrt{(x - xc)^2 + (y - yc)^2}
        
        \Theta=-\frac{1}{4}N(r^2 - R^2),\quad\quad r \leq R
        
        \Theta=-\frac{NR^2}{2}\ln(\frac{r}{R}),\quad\quad r > R
     
.. autoclass:: timml.circareasink.CircAreaSink
