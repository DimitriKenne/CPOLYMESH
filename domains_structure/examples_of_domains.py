# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 02:40:17 2023

Updated on June 30, 2024

@author: Dimitri Jordan Kenne

Authors:
-------
    - D. J. Kenne (Code Implementation) <dimitri.kenne@doctoral.uj.edu.pl>
    - Alvise Sommariva <alvise@math.unipd.it>
    - Marco Vianello  <marcov@math.unipd.it>

Thank you to:
------------ 
    Prof. Dr. hab. Leokadia Biales-Ciez for her help in developing this work.

License:
-------
This project is licensed under the MIT License - 
see the LICENSE file for details.

"""


# ------------------------------------------------------------------------------
# ------------------------------import packages---------------------------------
# pylint: disable=invalid-name
from scipy.interpolate import CubicSpline, BPoly
import cmath
import numpy as np
from domains_structure.polynomial_curves import PolyCurve, UnionPolyCurves
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------


def define_domain(example):
    '''
    Choose a prebuilt domain amongs the following:
        0. Unit circle
        1. Segment [-1,1]
        2. Polygon M
        3. Sun
        4. Ellipse 
        5. Union of circles
        6. Lune
        7. Cardioid
        8. 4 lenses
        9. Curve polygon
        10. Limacon
        11. Lissajous
        12. Egg
        13. Rhodonea
        14. Habenitch clover
        15. Bifolium
        16. Torpedo
        17. Double egg
        18. L. Sautereau Butterfly 1
        19. L. sautereau Butterfly  2
        20. Borromean curve

    Parameters
    ----------
    example : int
        A number between 0 and 10.

    Returns
    -------
    domain : PolyCurve
        See domain_structure.polynomial_curves.py for description.

    '''
    if example == 0:
        domain = circle(0, 1)
    elif example == 1:
        domain = segment(-1, 1)
    elif example == 2:
        sub_example = 9
        vertices, name = gallery_polygons(sub_example)
        domain = polygon(vertices, name)
    elif example == 3:
        domain = sun(0, 1)
    elif example == 4:
        domain = ellipse(0, 2, 1)
    elif example == 5:
        sub_example = 1
        centers, radii, name = gallery_union_circles(sub_example)
        domain = union_circles(centers, radii, name)
    elif example == 6:
        centers = [-1, 1]
        radii = [1.5, 1.5]
        domain = lune(centers, radii)
    elif example == 7:
        domain = cardioid(1)
    elif example == 8:
        domain = union_lenses(name='Butterfly')
    elif example == 9:
        sub_example = 3  # see gallery_curvpolygon below
        Sx, Sy, vertices, degrees, _ = gallery_curvpolygon(sub_example)
        domain = curve_polygon(Sx, Sy, vertices, degrees)
    elif example == 10:
        a = 7
        b = 5
        domain = limacon(a, b)
    elif example == 11:
        # Note: A lissajous curve is defined by:
        #     z(t) =x(t)+1j*y(t) where
        #     x(t) = av[0]*np.sin(wV[0]*t+phiV[0])
        #     y(t) = av[1]*np.sin(wV[1]*t+phiV[1])
        # with wV vector of positive integers
        aV = [1, 1]
        wV = [1, 2]
        phiV = [0, 0]
        domain = lissajous(aV, wV, phiV)
    elif example == 12:
        # NOTE: An egg is defined by:
        #     z(t) = a*(cos(t))**n*exp(1j*t), t in [0, 2*pi]
        # with n positive interger
        a = 1
        n = 3
        domain = egg(a, n)
    elif example == 13:
        # NOTE:
        # rhodonea: z(t)=a*cos(k*t)*exp(i*t), t in [0,2*pi]
        # 1. choose "k > 0" as integer.
        # 2. choose "a > 0"
        # * https://en.wikipedia.org/wiki/Rose_(mathematics)
        a = 1
        k = 4
        domain = rhodonea(a, k)
    elif example == 14:
        # NOTE:
        # n: parameter defining the boundary, with
        #          z(t)=(1+cos(n*t)+(sin(n*t))**2)*exp(i*t),  t in [0,2*pi]
        # with "n > 0" integer; e.g. choose "n=3"
        # It is a clover like domain.
        n = 3
        domain = habenicht_clover(n)
    elif example == 15:
        # NOTE:
        # the boundary is z(t)=cos(t)*(sin(t))^2)*exp(i*t),  t in [0,2*pi]
        domain = bifolium()
    elif example == 16:
        # NOTE:
        # the boundary is z(t)=a*cos(t)*cos(2*t).*exp(i*t),  t in [0,2*pi]
        a = 1
        domain = torpedo(a)
    elif example == 17:
        # NOTE:
        # the boundary is z(t)=(a/2)*(cos(2t)+1)*exp(i*t),  t in [0,2*pi]
        a = 1 
        domain = double_egg(a)
    elif example == 18:
        # Reference:
        #https://mathcurve.com/courbes2d/ornementales/ornementales.shtml
        domain = butterflysautereau(1)
    elif example == 19:
        domain = butterflysautereau(2)
    elif example == 20:
        domain = borromean(5)
    return domain
# ------------------------------------------------------------------------------


def segment(a=-1, b=1):
    """
    Create the segment [a,b]. 

    Parameters
    ----------
    a : complex, optional
        The first vertex of the segment. The default is -1.
    b : complex, optional
        The second vertex of the segment. The default is 1.

    Returns
    -------
    P: PolyCurve
      It represent the segment [a,b]. 
      See domain_structure.polynomial_curves.py for description.

    """
    name = f'segment [{a},{b}]'
    interval = (a, b)
    degree = 1
    curve_type = 'alg'

    # define # define a parametrization
    def param_segment(t):
        return t
    S = PolyCurve(param_segment, interval, degree, curve_type, name)
    return S


def circle(z, r):
    """
    Create the circle of radius r and center z

    Parameters
    ----------
    z : complex
        The center of the circle.
    r : float
        The radius.

    Returns
    -------
    P: PolyCurve
        The circle of radius r and center z. 
        See domain_structure.polynomial_curves.py for description.

    """
    name = f'circle of center {z} and radius {r}'
    interval = (0, 2*np.pi)
    degree = 1
    curve_type = "trig"

    # define a parametrization
    def param_circle(t, z=z, r=r):
        C = z+r*complex(np.cos(t), np.sin(t))
        return C

    return PolyCurve(param_circle, interval, degree, curve_type, name)


def ellipse(z, a, b):
    """
    Create a function that is the parametrisation of the complex ellipse of 
    centerz and semi axes a, b and defined on [0,2*np.pi]

    Parameters
    ----------
    z : complex
        The center of the circle.
    a, b : positive float
        semi axes.

    Returns
    -------
    P: PolyParamCurve
        The ellipse of center z and semi axes a, b. 
        See domain_structure.polynomial_curves.py for description.

    """
    name = f'Ellipse of center {z} and semi-axis {a}, {b}'
    interval = (0, 2*np.pi)
    degree = 1
    curve_type = "trig"

    # define a parametrization
    def param_ellipse(t, z=z, a=a, b=b):
        C = z+complex(a*np.cos(t), b*np.sin(t))
        return C

    return PolyCurve(param_ellipse, interval, degree, curve_type, name)


def cardioid(c=1):
    """
    Creat the cardioid of parameter c 

    Parameters
    ----------
    c : float, optional

    Returns
    -------
    P: PolyCurve
        It is the cardioid of parameter c. 
        See domain_structure.polynomial_curves.py for description.

    """
    name = 'Cardioid'
    interval = (0, 2*np.pi)
    degree = 2
    curve_type = "trig"

    # define a parametrization
    def param_cardioid(t, c=c):
        x_t = c*(1 - np.cos(t))*np.cos(t)
        y_t = c*(1 - np.cos(t))*np.sin(t)
        return complex(x_t, y_t)

    P = PolyCurve(param_cardioid, interval, degree, curve_type, name=name)
    return P


def polygon(vertices, name=None):
    """
    Create a polygon with vertices 

    Parameters
    ----------
    vertices : list or array
        Contains the vertices of the polygon. The first and the last element
        must be equal unless the domain is just a union of segments 

    Returns
    -------
    curves: UnionPolyCurves
       A polygon. See domain_structure.polynomial_curves.py for description.

    """
    if name == None:
        name = 'Polygon'
    curves = []
    N = len(vertices)
    for i in range(1, N):
        curves.append(segment(vertices[i-1], vertices[i]))

    return UnionPolyCurves(curves, name=name)


def curve_polygon(Sx, Sy, intervals_splines, degrees, name=None):
    '''
    Create a curve polygon

    Parameters
    ----------
    Sx : list of splines for the x-direction
        DESCRIPTION.
    Sy : list of splines for the x-direction
        DESCRIPTION.
    intervals_splines : list of tuple
        Each tuple intervals_splines[i] gives the breakpoints used to define 
        Sx+1j*Sy.
    degrees : list of int
        degrees[i] is the degree of the splines Sx[i] and Sy[i].
    name : TYPE, optional
        DESCRIPTION. The default is 'Curve polygon'.

    Returns
    -------
    PolyCurve
        See domain_structure.polynomial_curves.py for description..

    '''
    if name == None:
        name = 'Curve polygon'
    curves = []
    n = len(intervals_splines)
    for i in range(n):
        interval = intervals_splines[i]

        # Define the polynomial from the splines Sx and Sy
        def polynomial(t):
            return Sx[i](t) + 1j*Sy[i](t)

        degree = degrees[i]
        curve = PolyCurve(polynomial, interval, degree, 'alg')
        curves.append(curve)
    return UnionPolyCurves(curves, name)


def union_circles(Z, R, name='Union of circles'):
    """
    Create a union of circles

    Parameters
    ----------
    Z : list
        Its elements are the center of the circles.
    R : list
        Its elements are the radius of the circles.

    Returns
    -------
    P : UnionPolyParamCurve
         Contains the parameterisations of the different circles. 

    """

    curves = []
    N = len(R)
    for k in range(N):
        curves.append(circle(Z[k], R[k]))

    return UnionPolyCurves(curves, name)


def sun(center, r, beams=8, beams_lenght=0.5):
    '''
    Create a sun by drawing a circle centered at the center with radius r, 
    and attach beams as equispaced segments of length beams_length, positioned 
    orthogonally around the circle.

    Parameters
    ----------
    center : complex
        The sun's center.
    r : float
        The circle radious.
    beams : int, optional
        The number of beams. The default is 8.
    beams_lenght : float, optional
        The lenght of the beams. The default is 0.5.

    Returns
    -------
    UnionPolyCurves
        See domain_structure.polynomial_curves.py for description.

    '''
    name = 'sun'
    curves = [circle(center, r)]
    theta = np.linspace(0, 2*np.pi, beams+1)
    theta = theta[:-1]

    for k in range(beams):
        vertex_1 = center + r*np.exp(1j*theta[k])
        vertex_2 = center + (r+beams_lenght)*np.exp(1j*theta[k])
        curves.append(segment(vertex_1, vertex_2))
    return UnionPolyCurves(curves, name=name)


def lune(X, R):
    """
    Generate the parametrisations of the arcs of circles forming a moon

    Parameters
    ----------
    X : list
        Contains the circle centers .
    R : list
        Contains the circle radius.

    Returns
    -------
    P: UnionPolyCurve
        See domain_structure.polynomial_curves.py for description.

    """

    name = 'Lune'
    x, y = circle_intersect(X[0], R[0], X[1], R[1])
    t1 = min(arg_cplx_number(x-X[0]), arg_cplx_number(y-X[0]))
    t2 = max(arg_cplx_number(x-X[0]), arg_cplx_number(y-X[0]))
    s1 = min(arg_cplx_number(x-X[1]), arg_cplx_number(y-X[1]))
    s2 = max(arg_cplx_number(x-X[1]), arg_cplx_number(y-X[1]))

    if abs(t1-t2) < np.pi:
        c = t1
        t1 = t2
        t2 = c+2*np.pi
    if abs(s1-s2) > np.pi:
        c = s1
        s1 = s2
        s2 = c+2*np.pi
    f = circle(X[0], R[0])
    g = circle(X[1], R[1])
    f.interval = (t1, t2)
    g.interval = (s1, s2)

    return UnionPolyCurves([f, g], name)

# -----------------------------------------------------------------------------
# --------------------------A lens---------------------------------------------


def lens(X, R):
    """
    Create a lens

    Parameters
    ----------
    X : list
        Contains the circle centers .
    R : list
        Contains the circle radius.

    Returns
    -------
    P: UnionPolyParamCurve
        See domain_structure.polynomial_curves.py for description.

    """

    name = 'Lens'
    x, y = circle_intersect(X[0], R[0], X[1], R[1])
    t1 = min(arg_cplx_number(x-X[0]), arg_cplx_number(y-X[0]))
    t2 = max(arg_cplx_number(x-X[0]), arg_cplx_number(y-X[0]))
    s1 = min(arg_cplx_number(x-X[1]), arg_cplx_number(y-X[1]))
    s2 = max(arg_cplx_number(x-X[1]), arg_cplx_number(y-X[1]))

    if abs(t1-t2) > np.pi:
        c = t1
        t1 = t2
        t2 = c+2*np.pi
    if abs(s1-s2) > np.pi:
        c = s1
        s1 = s2
        s2 = c+2*np.pi
    f = circle(X[0], R[0])
    g = circle(X[1], R[1])
    f.interval = (t1, t2)
    g.interval = (s1, s2)

    return UnionPolyCurves([f, g], name)

# -----------------------------------------------------------------------------
# ------------------------Four lens--------------------------------------------


def union_lenses(X=[-1, -1j, 1, 1j, -1], R=[1, 1, 1, 1, 1], name=None):
    """
    Create many lenses

    Parameters
    ----------
    X : list
        Contains the circle centers . The default [-1, -1j, 1, 1j, -1].
    R : list
        Contains the circle radius. The default is [1,1, 1,1, 1]

    Returns
    -------
    P: UnionPolyParamCurve
        See domain_structure.polynomial_curves.py for description.

    """
    if name == None:
        name = f'{len(X)} lenses'
    curves = []
    for i in range(len(X)-1):
        f = lens(X[i:i+2], R[i:i+2]).curves
        curves = curves+f
    return UnionPolyCurves(curves, name)


def limacon(a, b, name='limacon'):
    '''
    Create a Limacon

    Parameters
    ----------
    a : float
    b : float   
    name : TYPE, optional
        The default is 'limacon'.

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The Limacon of parameters a and b

    Note:  1. for avoiding auto-intersections choose "b <= a";
         2. the domain is convex if "a >= 2*b";
         3. the domain is  concave if "b <= a < 2*b".
    '''
    def polynomial(t): return b/2 + a*np.exp(1j*t) + (b/2)*np.exp(2j*t)
    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, 2, 'trig')


def lissajous(aV, wV, phiV, name=None):
    '''
    Create a Lissajous curve

    Parameters
    ----------
    aV : list
    wV : list of positive intergers.
    phiV : list
    name : TYPE, optional
         The default is None.

    Returns
    -------
    PolyCurve (See domain_structure.polynomial_curves.py for description.)
        The Lissajous curve.

    '''
    if name is None:
        name = 'Lissajous'

    def polynomial(t):
        x = aV[1] * np.sin(wV[0]*t+phiV[0])
        y = aV[1] * np.sin(wV[1]*t+phiV[1])
        return complex(x, y)

    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, 1, 'trig', name)


def egg(a, n, name='Egg'):
    '''
    Create an Egg domain

    Parameters
    ----------
    a : float
    n : positive int    
    name : TYPE, optional
        The default is 'Egg'.

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The Egg of parameters a and n
    '''
    def polynomial(t): return a*(np.cos(t))**n*np.exp(1j*t)
    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, n+1, 'trig', name)


def double_egg(a, name='Double egg'):
    '''
    Create a double egg domain

    Parameters
    ----------
    a : float  
    name : TYPE, optional
        The default is 'Egg'.

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The double egg of parameters a and n
    '''
    def polynomial(t): return (a/2)*(np.cos(2*t)+1)*np.exp(1j*t)
    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, 3, 'trig', name)


def rhodonea(a, n, name='Rhodonea'):
    '''
    Create a Rhodonea

    Parameters
    ----------
    a : float
    n : positive int    
    name : TYPE, optional
        The default is 'Rhodonea'.

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The rhodonea of parameters a and n
    '''
    def polynomial(t): return a*np.cos(n*t)*np.exp(1j*t)
    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, n+1, 'trig', name)


def habenicht_clover(n, name='Habenicht clover'):
    '''
    Create a Habenicht_clover

    Parameters
    ----------
    n : positive int    
    name : TYPE, optional
        The default is 'Habenicht_clover'.

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The Habenicht_clover of parameters n
    '''
    def polynomial(t): return (1+np.cos(n*t) + np.sin(n*t)**2)*np.exp(1j*t)
    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, 2*n+1, 'trig', name)


def bifolium(name='Bifolium'):
    '''
    Create a Bifolium

    Parameters
    ----------  
    name : TYPE, optional
        The default is 'Bifolium'.

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The Bifolium of parameters n
    '''
    def polynomial(t): return np.cos(t)*np.sin(t)**2*np.exp(1j*t)
    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, 4, 'trig', name)


def torpedo(a, name='Torpedo'):
    '''
    Create a Torpedo

    Parameters
    ----------  
    a: float
    name : TYPE, optional
        The default is 'Torpedo'.

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The Torpedo of parameters n
    '''
    def polynomial(t): return a*np.cos(t)*np.cos(2*t)*np.exp(1j*t)
    interval = (0, 2*np.pi)

    return PolyCurve(polynomial, interval, 4, 'trig', name)


def butterflysautereau(a, name='L. sautereau Butterfly'):
    '''
    Create a butterflysautereau 
    Reference: 
    https://mathcurve.com/courbes2d/ornementales/ornementales.shtml

    parameters
    a : 1 or 2

    Returns
    -------
    Polycurve (See domain_structure.polynomial_curves.py for description.)
        The L. Sautereau butterfly of parameter a
    '''

    interval = (0, 2*np.pi)
    if a==1:
        def polynomial1(t): return (np.sin(5*t)+ 3*np.cos(t))*np.exp(1j*t)
        def polynomial2(t): return (np.sin(5*t)- 3*np.cos(t))*np.exp(1j*t)
        left_wing = PolyCurve(polynomial1, interval, 6, 'trig')
        right_wing = PolyCurve(polynomial2, interval, 6, 'trig')
        domain = UnionPolyCurves([left_wing, right_wing], name)
    elif a==2:
        def polynomial(t): return (-3*np.cos(2*t)+ np.sin(7*t)-1)*np.exp(1j*t) 
        domain = PolyCurve(polynomial, interval, 8, 'trig', name)
    return domain


def borromean(ndisks):
    '''
    Create a Borromean ring as an instance of PolyCurve
    
    Reference  https://en.wikipedia.org/wiki/Borromean_rings
    ---------
    
    Parameters
    ----------
    n : int
        number of disks making the borromean curve.
    return
    ------
    domain : PolyCurve (See domain_structure.polynomial_curves.py for description.)
    '''
    
    th=2*np.pi/ndisks
    th0=(np.pi-th)/2
    curves = []
    name = 'Borromean curve'
    
    for k in range(ndisks):
        thL=th0+k*th
        interval = (0, 2*np.pi)
        polynomial = lambda t: np.sqrt(3)*np.exp(1j*t)+complex(np.cos(thL), np.sin(thL))
        curves.append(PolyCurve(polynomial, interval, 1, 'trig'))
    print( f'\n {len(curves)} \n \n' )
    domain = UnionPolyCurves(curves, name)
    return domain
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


def arg_cplx_number(x):
    """
    Find the argument of a complex number in [0,2*np.pi)

    Parameters
    ----------
    x : complex

    Returns
    -------
    t : float

    """

    t = cmath.phase(x)
    if t < 0:
        t = 2*np.pi + t
    return t
# ----------------


def circle_intersect(c1, r1,  c2, r2):
    """
    Find the intersection between two complax planar circles

    Parameters
    ----------
    c1 : complex
        The first circle center.
    r1 : float
        The first circle radius.
    c2 : complex
        The second circle center.
    r2 : float
        The second circle radius.

    Returns
    -------
    z1 : complex
        The first intersection point.
    z2 : complex
        The second intersection point.

    """

    u = c2 - c1
    d = abs(u)

    if d > r1+r2 or d == 0:
        print('There is no intersection')
        return None

    v = u/d
    x = (d**2 - r2**2 + r1**2)/(2*d)
    a = ((-d+r2-r1)*(-d-r2+r1)*(-d+r2+r1)*(d+r2+r1))**0.5/d
    y = a/2

    z1 = c1 + v*(complex(x, y))
    z2 = c1 + v*(complex(x, -y))
    return z1, z2

# -----------------------------------------------------------------------------


def gallery_union_circles(example=0):
    # COMPLEX DOMAIN AS UNION OF DISKS.
    # IMPORTANT: centers elements are complex numbers.

    name = None
    if example == 0:
        disks = np.array([
            [4.172670690843695e-01, 6.443181301936917e-01, 3.630885392869130e-01],
            [4.965443032574213e-02, 3.786093826602684e-01, 5.468057187389680e-01],
            [9.027161099152811e-01, 8.115804582824772e-01, 5.211358308040015e-01],
            [9.447871897216460e-01, 5.328255887994549e-01, 2.315943867085238e-01],
            [5.085086553811270e-01, 4.886089738035791e-01, 3.786806496411588e-01],
            [5.107715641721097e-01, 5.785250610234389e-01, 7.126944716789141e-01],
            [8.176277083222621e-01, 2.372835797715215e-01, 5.004716241548430e-01],
            [7.948314168834530e-01, 4.588488281799311e-01, 4.710883745419393e-01]
        ])

        centers = [complex(x[0], x[1]) for x in disks[:, 0:2]]
        radii = disks[:, 2]

    elif example == 1:
        disks = np.array([
            [5.085086553811270e-01, 4.886089738035791e-01, 9.786806496411588e-01],
            [5.107715641721097e-01, 5.785250610234389e-01, 7.126944716789141e-01],
            [8.176277083222621e-01, 2.372835797715215e-01, 5.004716241548430e-01],
            [7.948314168834530e-01, 4.588488281799311e-01, 4.710883745419393e-01],
            [4.172670690843695e-01, 6.443181301936917e-01, 9.630885392869130e-01],
            [4.965443032574213e-02, 3.786093826602684e-01, 5.468057187389680e-01],
            [9.027161099152811e-01, 8.115804582824772e-01, 5.211358308040015e-01],
            [9.447871897216460e-01, 5.328255887994549e-01, 2.315943867085238e-01],
            [4.908640924680799e-01, 3.507271035768833e-01, 4.888977439201669e-01],
            [4.892526384000189e-01, 9.390015619998868e-01, 6.240600881736895e-01],
            [3.377194098213772e-01, 8.759428114929838e-01, 6.791355408657477e-01],
            [9.000538464176620e-01, 5.501563428984222e-01, 3.955152156685930e-01],
            [3.692467811202150e-01, 6.224750860012275e-01, 3.674366485444766e-01],
            [1.112027552937874e-01, 5.870447045314168e-01, 9.879820031616328e-01],
            [7.802520683211379e-01, 2.077422927330285e-01, 3.773886623955214e-02],
            [3.897388369612534e-01, 3.012463302794907e-01, 8.851680082024753e-01],
            [2.416912859138327e-01, 4.709233485175907e-01, 9.132868276392390e-01],
            [4.039121455881147e-01, 2.304881602115585e-01, 7.961838735852120e-01],
            [9.645452516838859e-02, 8.443087926953891e-01, 9.871227865557430e-02],
            [1.319732926063351e-01, 1.947642895670493e-01, 2.618711838707161e-01],
            [9.420505907754851e-01, 2.259217809723988e-01, 3.353568399627965e-01],
            [9.561345402298023e-01, 1.707080471478586e-01, 6.797279513773380e-01],
            [5.752085950784656e-01, 2.276642978165535e-01, 1.365531373553697e-01],
            [5.977954294715582e-02, 4.356986841038991e-01, 7.212274985817402e-01],
            [2.347799133724063e-01, 3.111022866504128e-01, 1.067618616072414e-01],
            [3.531585712220711e-01, 9.233796421032439e-01, 6.537573486685596e-01],
            [8.211940401979591e-01, 4.302073913295840e-01, 4.941739366392701e-01],
            [1.540343765155505e-02, 1.848163201241361e-01, 7.790517232312751e-01],
            [4.302380165780784e-02, 9.048809686798929e-01, 7.150370784006941e-01],
            [1.689900294627044e-01, 9.797483783560852e-01, 9.037205605563163e-01],
            [6.491154749564521e-01, 4.388699731261032e-01, 8.909225043307892e-01],
            [7.317223856586703e-01, 1.111192234405988e-01, 3.341630527374962e-01],
            [6.477459631363067e-01, 2.580646959120669e-01, 6.987458323347945e-01],
            [4.509237064309449e-01, 4.087198461125521e-01, 1.978098266859292e-01],
            [5.470088922863450e-01, 5.948960740086143e-01, 3.054094630463666e-02],
            [2.963208056077732e-01, 2.622117477808454e-01, 7.440742603674624e-01],
            [7.446928070741562e-01, 6.028430893820830e-01, 5.000224355902009e-01],
            [1.889550150325445e-01, 7.112157804336829e-01, 4.799221411460605e-01],
            [6.867754333653150e-01, 2.217467340172401e-01, 9.047222380673627e-01],
            [1.835111557372697e-01, 1.174176508558059e-01, 6.098666484225584e-01],
            [3.684845964903365e-01, 2.966758732183269e-01, 6.176663895884547e-01],
            [6.256185607296904e-01, 3.187783019258823e-01, 8.594423056462123e-01],
            [7.802274351513768e-01, 4.241667597138072e-01, 8.054894245296856e-01],
            [8.112576886578526e-02, 5.078582846611182e-01, 5.767215156146851e-01],
            [9.293859709687300e-01, 8.551579709004398e-02, 1.829224694149140e-01],
            [7.757126786084023e-01, 2.624822346983327e-01, 2.399320105687174e-01],
            [4.867916324031724e-01, 8.010146227697388e-01, 8.865119330761013e-01],
            [4.358585885809191e-01, 2.922027756214629e-02, 2.867415246410610e-02],
            [4.467837494298063e-01, 9.288541394780446e-01, 4.899013885122239e-01],
            [3.063494720165574e-01, 7.303308628554529e-01, 1.679271456822568e-01],
            [5.085086553811270e-01, 4.886089738035791e-01, 9.786806496411588e-01],
            [5.107715641721097e-01, 5.785250610234389e-01, 7.126944716789141e-01],
            [8.176277083222621e-01, 2.372835797715215e-01, 5.004716241548430e-01],
            [7.948314168834530e-01, 4.588488281799311e-01, 4.710883745419393e-01]
        ])

        Ndisks = 15
        centers = [complex(x[0], x[1]) for x in disks[:Ndisks, 0:2]]
        radii = disks[:Ndisks, 2]

    elif example == 2:
        M = 20
        angles = np.linspace(0, 2 * np.pi, M)

        r1 = 2
        cents_1 = r1 * (np.cos(angles) + 1j * np.sin(angles))
        rs1 = (r1 / 4) * np.ones(len(cents_1))

        r2 = 4
        cents_2 = r2 * (np.cos(angles) + 1j * np.sin(angles))
        rs2 = (r2 / 4) * np.ones(len(cents_2))

        centers = np.concatenate((cents_1, cents_2))
        radii = np.concatenate((rs1, rs2))

    elif example == 3:
        M = 45
        t = np.linspace(0, 5, M)
        y = 2 * t
        x = 2.5 * np.cos(2 * t)
        cents_1 = x + y*1j
        x = 2.5 * np.sin(2 * t)
        cents_2 = x + y*1j
        centers = np.concatenate((cents_1, cents_2))
        radii = 0.3 * np.ones(len(centers))

    elif example == 4:
        name = 'Unit circle'
        centers = np.array([[0, 0]])
        radii = np.array([1])

    else:
        M = 100
        cents = np.random.rand(M, 2)
        radii = np.random.rand(M)
        centers = centers[:, 0] + 1j * cents[:, 1]

    return centers, radii, name


def gallery_polygons(example):
    # IMPORTANT: first and last vertices are equal.

    vertices = []
    name = None
    if example == 0:
        name = 'Unit square [0,1]^2'
        k = 1
        vertices = k * np.array([[0, 0], [1, 0], [1, 1], [0, 1]])

    elif example == 1:
        name = 'convex polygon'
        vertices = np.array([[0.1, 0], [0.7, 0.2], [1, 0.5], [0.75, 0.85],
                             [0.5, 1], [0, 0.25]])

    elif example == 2:
        name = 'Non convex polygon'
        vertices = (1/4) * np.array([[1, 0], [3, 2], [3, 0], [4, 2],
                                     [3, 3], [3, 0.85*4], [2, 4], [0, 3],
                                     [1, 2]])

    elif example == 3:
        name = 'Polygon'
        vertices = (1/4) * np.array([[1, 0], [3, 2], [3, 0], [4, 2], [3, 3],
                                     [3, 4], [2, 4], [0, 3], [1, 2]])

    elif example == 4:
        name = 'Polygon'
        vertices = (1/4) * np.array([[0, 0], [1, 2], [2, 0], [3, 2], [4, 0],
                                     [4, 4], [3, 2], [2, 4], [1, 2], [0, 4]])

    elif example == 5:
        name = 'Unit square [-1,1]^2'
        vertices = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])

    elif example == 6:
        name = 'Unit triangle: [0 0; 1 0; 1 1]'
        vertices = np.array([[0, 0], [1, 0], [1, 1]])

    elif example == 7:
        name = 'Pseudo-circle'
        N = 20
        a, b = 0, 2 * np.pi
        h = (b - a) / N
        theta = np.arange(a, b, h)
        x = np.cos(theta)
        y = np.sin(theta)
        x = 0.5 * x + 0.5
        y = 0.5 * y + 0.5
        vertices = np.column_stack((x, y))

    elif example == 8:
        name = 'M'
        vertices = np.array([
            [-1, 1], [-1, 0], [-1, -1], [-0.6, -1], [-0.6, 0], [-0.6, 0.6],
            [0, 0], [0.6, 0.6], [0.6, 0], [0.6, -1], [1, -1], [1, 0], [1, 1],
            [0.6, 1], [0, 0.4], [-0.6, 1], [-1, 1]
        ])

    else:
        name = 'Non convex polygon'
        vertices = (1/4)*np.array([[1, 0], [3, 2], [3, 0], [4, 2], [3, 3],
                                   [3, 0.85*4], [2, 4], [0, 3], [1, 2]])

    # Adding the first row to the end of the array
    vertices = np.vstack((vertices, vertices[0]))

    # Converting vertices to complex numbers
    vertices = vertices[:, 0] + 1j * vertices[:, 1]

    return vertices, name


def gallery_curvpolygon(subexample):
    P = []
    vertices = None
    domain_name = ''

    # Case 1
    if subexample == 0:
        domain_name = 'Unit square [0,1]^2'
        k = 1
        vertices = k * np.array([[0, 0], [1, 0], [1, 1], [0, 1]])
        # The first and the last vertices coincide
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, vertices.shape[0]]])

    # Case 2
    elif subexample == 1:
        domain_name = 'Convex polygon'
        vertices = np.array([[0.1, 0], [0.7, 0.2], [1, 0.5], [0.75, 0.85],
                             [0.5, 1], [0, 0.25]])
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, vertices.shape[0]]])

    # Case 3
    elif subexample == 2:
        domain_name = 'Non convex polygon'
        vertices = (1/4) * np.array([[1, 0], [3, 2], [3, 0], [4, 2], [3, 3],
                                     [3, 0.85*4], [2, 4], [0, 3], [1, 2]])
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, 7], [2, vertices.shape[0]]])

    # Case 4
    elif subexample == 3:
        domain_name = 'Polygon'
        vertices = (1/4) * np.array([[1, 0], [3, 2], [3, 0], [4, 2], [3, 3],
                                     [3, 4], [2, 4], [0, 3], [1, 2]])
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[2, 6], [4, vertices.shape[0]]])

    # Case 5
    elif subexample == 4:
        domain_name = 'Polygon'
        vertices = (1/4) * np.array([[0, 0], [1, 2-0.1], [2, 0], [3, 2-0.1],
                                     [4, 0], [4, 4], [3, 2], [2, 4],
                                     [1, 2+0.1], [0, 4]])
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[2, 6], [3, 10], [2, vertices.shape[0]]])

    # Case 6
    elif subexample == 5:
        domain_name = 'Unit square [-1,1]^2'
        vertices = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, vertices.shape[0]]])

    # Case 7
    elif subexample == 6:
        domain_name = 'Unit triangle: [0 0; 1 0; 1 1]'
        vertices = np.array([[0, 0], [1, 0], [1, 1]])
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, vertices.shape[0]]])

    # Case 8
    elif subexample == 7:
        domain_name = 'Pseudo-circle'
        N = 40
        a, b = 0, 2 * np.pi
        h = (b - a) / N
        theta = np.arange(a, b, h)
        x = np.cos(theta)
        y = np.sin(theta)
        x = 0.5 * x + 0.5
        y = 0.5 * y + 0.5
        vertices = np.column_stack((x, y))
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, vertices.shape[0]]])

    # Case 9
    elif subexample == 8:
        domain_name = 'Wings'
        X0 = 0.3 * np.array([0, 1, 2, 3, 4])
        X = np.concatenate((X0, [3.5, 2, -1]))
        Y = np.array([0, 2, 4.5, 2, 0, 3.5, 5, 3.5])
        vertices = np.column_stack((X, Y))
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, vertices.shape[0]]])

    else:
        domain_name = 'Non convex polygon'
        vertices = (1/4) * np.array([[1, 0], [3, 2], [3, 0], [4, 2], [3, 3],
                                     [3, 0.85*4], [2, 4], [0, 3], [1, 2]])
        vertices = np.vstack((vertices, vertices[0]))

        spline_type = 'not-a-knot'
        spline_parms = np.array([[4, vertices.shape[0]]])

    X = vertices[:, 0]
    Y = vertices[:, 1]

    # Determine the spline boundary via the vectors of splines "Sx" and "Sy".
    Sx, Sy, intervals_splines, spline_orders = compute_spline_boundary(X, Y)

    vertices = X + 1j*Y

    degrees = spline_orders - 1
    return Sx, Sy, intervals_splines, degrees, domain_name


# ------------------------------------------------------------------------------


def compute_spline_boundary(X, Y, spline_parms=None, spline_type=None):
    """
    This function computes the parametrical description of the boundary in the 
    "x" and "y" variable, via many splines, with possible different orders.

    Args:
    X (array_like): Sampling points in the x-direction.
    Y (array_like): Sampling points in the y-direction.
    spline_parms (ndarray, optional): 
        Array specifying spline order and block. 
        [spline_parms(i,0)=2] : piecewise linear splines.
            [spline_parms(i,0)=4] : cubic splines
                         (depending on "spltypestring" in input).
            [spline_parms(i,0)=k] : in the case k is not 2 or 4 it
                          chooses the k-th order splines (for additional
                          help digit "help spapi" in matlab shell).

            [spline_parms(:,2)]: vector of final components of a block.

            example:

            "spline_parms= np.array([[2, 31], [4, 47], [8, 67]])" means that 
            from the 1st vertex to the 31th vertex we have an order 2 spline 
            (piecewise linear), from the 32th vertex to the 47th we use a 4th 
            order spline (i.e. a cubic and periodic spline by default), from
             the 48th to the 67th (and final!) we use an 8th order spline.

        Default is None.

    spline_type (str, optional): Type of spline for cubic splines. 
    Default is None.

    Returns:
    Sx (list): List of splines describing the boundary in the x-direction.
    Sy (list): List of splines describing the boundary in the y-direction.
    """

    # Set default spline parameters
    if spline_parms is None:
        spline_parms = np.array([[4, len(X)]])
    if spline_parms.shape[0] == 0:
        spline_parms = np.array([[4, len(X)]])

    if spline_type is None:
        spline_type = 'periodic'

    Sx = []
    Sy = []
    intervals_splines = []  # will contain the breakpoints for each spline inside Sx

    for i in range(spline_parms.shape[0]):
        if i == 0:
            imin = 0
        else:
            imin = spline_parms[i - 1, 1]
        imax = spline_parms[i, 1]

        tL = np.arange(imin, imax)
        xL = X[imin:imax]
        yL = Y[imin:imax]

        SxL, SyL = compute_parametric_spline(tL, xL, yL, spline_parms[i, 0],
                                             spline_type)

        Sx.append(SxL)
        Sy.append(SyL)
        intervals_splines.append((tL[0], tL[-1]))

    return Sx, Sy, intervals_splines, spline_parms[:, 0]


def compute_parametric_spline(s, x, y, spline_order, spline_type='periodic'):
    """
    Compute parametric spline relevant parameters so that a point at the 
    boundary of the domain has coordinates (x(s), y(s)).

    Args:
    s (array_like): Parameter data.
    x (array_like): Determines spline x(s) interpolating (s, x).
    y (array_like): Determines spline y(s) interpolating (s, y).
    spline_order (int): Spline order (degree+1).
    spline_type (str, optional): Type of spline. Default is 'periodic'.

    Returns:
    ppx (CubicSpline or BPoly): Spline x(t) data.
    ppy (CubicSpline or BPoly): Spline y(t) data.
    """

    if spline_order == 4:
        # Cubic splines: using CubicSpline
        ppx = CubicSpline(s, x, bc_type=spline_type)
        ppy = CubicSpline(s, y, bc_type=spline_type)
    else:
        # Non-cubic splines: using BPoly
        ppx = BPoly.from_derivatives(s, x.reshape(1, -1), k=spline_order-1)
        ppy = BPoly.from_derivatives(s, y.reshape(1, -1), k=spline_order-1)

    return ppx, ppy
