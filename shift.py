import itertools
import pylab as plt
import matplotlib.cm as cm
import numpy as np
from   scipy.ndimage.filters import gaussian_filter
from   numpy.linalg import solve
import math
import matplotlib.pyplot


def vimage(cat1, cat2, dmax, psize, fwhm):
    """Given two position catalogues of stars, each a numpy array of the form
    [[x0,y0],[x1,y1], ...], assumed to be linearly translated from each this
    works out the (x,y) translation vector that gets from cat1 to cat2.
    Method: displacement vectors between stars in each field are stored in a
    square grid centred on (0,0). The idea is the real displacement should
    show up as a cluster of pixels with vector hits. The size of the grid used
    is determined by the maximum shift allowed, dmax, and the pixel size,
    psize, (expressed in whatever units are used for the catalogues, usually
    CCD pixel numbers presumably.)
    Arguments::
      cat1 : first position catalogue
      cat2 : second position catalogue
      dmax : maximum distance shift in either x or y
      psize : size of pixels used for sub-image used to compute shift
      fwhm : FWHM of blurring used to locate peak in image
    """

    NHALF  = int(dmax/psize)
    NSIDE  = 2*NHALF+1
    mshift = (NHALF+0.5)*psize
    img = np.zeros((NSIDE,NSIDE))
    x2s, y2s = cat2[:,0], cat2[:,1]
    for x1, y1 in cat1:
        ok = (x2s > x1-mshift) & (x2s < x1+mshift) & \
             (y2s > y1-mshift) & (y2s < y1+mshift)
        for x2, y2 in cat2[ok]:
            ix = NHALF+int(round((x2-x1)/psize))
            iy = NHALF+int(round((y2-y1)/psize))
            img[iy,ix] += 1

    # smooth image
    img = gaussian_filter(img,fwhm/psize/2.3548,mode='constant')

    # identify maximum pixel
    ind = np.arange(NSIDE)
    ix, iy = np.meshgrid(ind, ind)
    peak = img == img.max()
    #if len(ix[peak]) > 1:
    #    raise Exception("Found more than one maximum pixel")

    # now have first approximation to the shift
    ixp = ix[peak][0]
    iyp = iy[peak][0]
    xp = psize*(ixp-NHALF)
    yp = psize*(iyp-NHALF)
    if ixp == 0 or ixp == NSIDE-1 or iyp == 0 or iyp == NSIDE-1:
        # max pixel at edge of array. Just return pixel position
        # as "refined" position
        xr = xp
        yr = yp

    else:
        # Make a quadratic approx to refine the peak position.
        # Estimate first and second partial derivatives from
        # 3x3 pixels centred on peak
        fx  = (img[iyp,ixp+1] - img[iyp,ixp-1])/2.
        fy  = (img[iyp+1,ixp] - img[iyp-1,ixp])/2.
        fxx = img[iyp,ixp-1] + img[iyp,ixp+1] - 2*img[iyp,ixp]
        fyy = img[iyp-1,ixp] + img[iyp+1,ixp] - 2*img[iyp,ixp]
        fxy = (img[iyp+1,ixp+1] + img[iyp-1,ixp-1] -
               img[iyp+1,ixp-1] - img[iyp-1,ixp+1])/4.
        b   = np.array((fx,fy)).T
        A   = np.array(((fxx,fxy),(fxy,fyy)))
        x   = solve(A,b)
        xr  = xp - psize*x[0]
        yr  = yp - psize*x[1]
    return (img, xp,yp,xr,yr)

def match(cat1, cat2, xs, ys, mmax):
    """
    Carries out pair matching between two catalogues
    with cat2 shifted relative to cat1 by (xs,ys).
    Arguments::
      cat1 : first catalogues of (x,y) pairs
      cat2 : second catalogues of (x,y) pairs
      xs   : Shift in X of cat2 relative to cat1
      ys   : Shift in Y of cat2 relative to cat1
      mmax : maximum shift in X and Y to look for matches
    Returns::
      nmatch : number of unique matches found
      inds   : indices of matching stars in cat1 for each cat2 star; -1
               if not matched. Hence len(inds[inds > -1]) == nmatch.
    """

    nmatch = 0
    x1s, y1s  = cat1[:,0], cat1[:,1]
    ind1 = np.arange(len(cat1))
    ind2 = np.empty(len(cat2),dtype=int)
    ind2.fill(-1)
    for i2, p2 in enumerate(cat2):
        x2, y2 = p2
        ok = (x1s > x2-xs-mmax) & (x1s < x2-xs+mmax) & \
             (y1s > y2-ys-mmax) & (y1s < y2-ys+mmax)
        if len(x1s[ok]) == 1:
            nmatch += 1
            ind2[i2] = ind1[ok][0]

    return (nmatch, ind2)

if __name__ == '__main__':

    # generate artificial catalogue

    psize  = 0.5
    fwhm   = 4.
    dmax   = 20.
    mmax   = 3.

    NSTARS = 1000
    x0 = np.random.uniform(high=1000., size=NSTARS)
    y0 = np.random.uniform(high=1000., size=NSTARS)
    

    # catalogue 1: all of the stars
    x1 = x0
    y1 = y0
    cat1 = np.array(zip(x1, y1))
    print(cat1)
    print(cat1.shape)

    d1 = []
    d2 = []
    for i in range(100):

        #print '\n---------------------------------'
        # catalogue 2: a different 50-100% of the stars
        n2 = int(len(x0)*(0.5+0.5*np.random.uniform()))
        x2 = x0[:n2].copy()
        y2 = y0[:n2].copy()

        # Perturb catalogue 2 with a displacement and some
        # small jitter [should also add some rotation]
        xs, ys = np.random.uniform(-5.,5.), np.random.uniform(-5.,5.)
        x2 += xs + np.random.normal(scale=0.7, size=len(x2))
        y2 += ys + np.random.normal(scale=0.7, size=len(y2))
        print(x2)
        cat2 = np.array(zip(x2, y2))

        xp,yp,xr,yr = vimage(cat1, cat2, dmax, psize, fwhm)
        print('x=',xs, xp, xr, 'y=', ys, yp, yr)

        d1.append(np.sqrt((xp-xs)**2+(yp-ys)**2))
        d2.append(np.sqrt((xr-xs)**2+(yr-ys)**2))

        nmatch, inds = match(cat1, cat2, xp, yp, mmax)
        #print nmatch, inds

    d1 = np.array(d1)
    d2 = np.array(d2)

    print(d1.mean(), d2.mean())