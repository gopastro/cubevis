"""Liang-Barsky algorithm to find the intersection of
straight line with rectangle.
http://www.cs.helsinki.fi/group/goa/viewing/leikkaus/intro.html
"""
import math

def yofline(m, a, b, x):
    """given slope m, and point (a,b) that line
    goes through and a new point x, determines new y"""
    return m*(x-a)+b

def dotproduct(n, Q):
    return n[0]*Q[0] + n[1]*Q[1]

def getT(x1, x2, A):
    if (x2-x1) == 0.0:
        return 10.0
    else:
        return (A-x1)/(x2-x1)
    


def line_rectangle_intersection(L, R, B, T, a, b, angle):
    """Given a rectangle with pixel edges defined by L(left),
    R(right), B(bottom), and T(top), and a line given by a point
    (a, b) with angle, angle wrt to X-axis, this routine will
    return the two intersection points in pixels to the rectangle
    @param L: left edge of rectangle in pixels
    @type L: float
    @param R: right edge of rectangle in pixels
    @type R: float
    @param B: bottom edge of rectangle in pixels
    @type B: float
    @param T: top edge of rectangle in pixels
    @type T: float
    @param a: x point of line
    @type a: float
    @param a: y point of line
    @type a: float
    @param angle: angle with respect to x-xaxis in degrees
    @type angle: float
    @return two tuples of (x,y) of intersecting points
    """

    L = float(L) #left edge
    R = float(R) #right edge
    B = float(B) #bottom edge
    T = float(T) #top edge

    #point in line
    a = float(a)
    b = float(b)
    m = math.tan(math.radians(float(angle)))

    #first point in line
    x1 = a - (R-L)
    y1 = yofline(m, a, b, x1)

    #second point in line
    x2 = a + (R-L)
    y2 = yofline(m, a, b, x2)

    tmin = 0.0
    tmax = 1.0

    QminusP = ((x2-x1), (y2-y1))
    #print "QminusP", QminusP
    #n->normal vector (outward of edge)
    #with edge defined by equation ax+by+c = 0
    #then n = (a,b)
    #t-values and determine entering or exiting..

    #left edge defined by x-L=0, hence n = (-1,0)

    tL = getT(x1, x2, L)
    if tL>= tmin and tL<=tmax:
        #valid point
        if dotproduct(QminusP, (-1,0)) < 0.0:
            #entering
            tmin = tL
        else:
            #exiting
            tmax = tL
    #print "L", tmin, tmax
    #right edge defined by x-R=0, hence n = (1,0)
    tR = getT(x1, x2, R)
    if tR>= tmin and tR<=tmax:
        #valid point
        if dotproduct(QminusP, (1,0)) < 0.0:
            #entering
            tmin = tR
        else:
            tmax = tR
    #print "R", tmin, tmax            
    #bottom edge defined by y-B=0, hence n = (0, -1)
    tB = getT(y1, y2, B)
    if tB>= tmin and tB<=tmax:
        #valid point
        if dotproduct(QminusP, (0,-1)) < 0.0:
            #entering
            tmin = tB
        else:
            tmax = tB
    #print "B", tmin, tmax
    #top edge defined by y-T=0, hence n = (0, 1)
    tT = getT(y1, y2, T)
    if tT>= tmin and tT<=tmax:
        #valid point
        #print dotproduct(QminusP, (0,1))
        if dotproduct(QminusP, (0,1)) < 0.0:
            #entering
            tmin = tT
        else:
            tmax = tT    
    #print "T", tmin, tmax
    #print tL, tR, tB, tT
    #print tmin, tmax
    minP = ((x1 + (x2-x1)*tmin), (y1+(y2-y1)*tmin))
    maxP = ((x1 + (x2-x1)*tmax), (y1+(y2-y1)*tmax))
    return minP, maxP
    
if __name__ == "__main__":
    i1, i2 = line_rectangle_intersection(0.0, 180., 0.0, 150., 75., 80., -89.9999999)
    print(i1, i2)
    
