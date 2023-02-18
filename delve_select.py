#!/usr/bin/env python
"""
Select objects from DELVE catalogs at Fermilab.
"""
__author__ = "Alex Drlica-Wagner"
import os
import pandas as pd
import numpy as np
import pylab as plt
import healpy as hp
import fitsio

NSIDE = 32
BANDS = ['g','r','i','z']
DIRNAME = '/data/des91.b/data/kadrlica/projects/delve/cat'
BASENAME = 'cat_hpx_%05d.fits'
HPXNAME = '%s/hpx_%s_%05d.fits'
RELEASES = {
    'dr1' :os.path.join(DIRNAME,'dr1'),
    'dr2': os.path.join(DIRNAME,'dr2'),
}
COLUMNS = None
#COLUMNS = ['QUICK_OBJECT_ID','RA','DEC']

def ang2disc(nside, lon, lat, radius, inclusive=False, fact=4, nest=False):
    """
    Wrap `query_disc` to use lon, lat, and radius in degrees.
    
    Parameters
    ----------
    nside : pixel resolution
    lon   : longitude (deg)
    lat   : latitude (deg)
    radius: disc radius (deg)
    """
    vec = hp.ang2vec(lon,lat,lonlat=True)
    return hp.query_disc(nside,vec,np.radians(radius),inclusive,fact,nest)

def angsep(lon1,lat1,lon2,lat2):
    """
    Angular separation (deg) between two sky coordinates.
    Borrowed from astropy (www.astropy.org)

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1],
    which is slighly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    lon1,lat1 = np.radians([lon1,lat1])
    lon2,lat2 = np.radians([lon2,lat2])
    
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1,num2), denominator))

def find_filenames(ra,dec,radius=1.0,release='dr1',style='cat'):
    """ Find catalog filenames. """
    pixels = ang2disc(NSIDE,ra,dec,radius,inclusive=True,nest=False)
    dirname = RELEASES[release]
    if style == 'cat':
        filenames = [os.path.join(dirname,style,BASENAME%p) for p in pixels]
    elif style == 'hpx':
        filenames = []
        for band in BANDS:
            for p in pixels:
                filenames.append(os.path.join(dirname,style,HPXNAME%(band,band,p)))

    filenames = [f for f in filenames if os.path.exists(f)]
    return filenames

def load_catalogs(ra,dec,radius=1.0,release='dr1',style='cat',verbose=False):
    """ Load catalogs."""
    filenames = find_filenames(ra,dec,radius,release,style=style)
    if args.verbose:
        print("Catalog filenames:")
        for f in filenames:
            print("  "+f)

    catalogs = []
    for f in filenames:
        if verbose: print("Loading %s..."%os.path.basename(f))
        cat = fitsio.read(f,columns=COLUMNS)
        sel = np.ones(len(cat),dtype=bool)
        if radius is not None:
            sep = angsep(ra,dec,cat['RA'],cat['DEC'])
            sel = sep < radius
        catalogs.append(cat[sel])

    return np.concatenate(catalogs)

def write_catalog(outfile,catalog):
    if args.verbose: print("Writing %s..."%outfile)

    if outfile.endswith(('.fits','.fz')):
        fitsio.write(outfile,catalog,clobber=True)
    elif outfile.endswith(('.csv','.txt','.dat')):
        df = pd.DataFrame(catalog.byteswap().newbyteorder())
        sep = ',' if outfile.endswith('.csv') else ' '
        df.to_csv(args.outfile,index=False,sep=sep)
    else:
        raise Exception("Unrecognized output file format: %s"%outfile)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ra',type=float,
                        help = 'Right ascension (deg)')
    parser.add_argument('dec',type=float,
                        help = 'Declination (deg)')
    parser.add_argument('-r','--radius',default=180, type=float,
                        help = 'selection radius (arcsec)')
    parser.add_argument('--release',choices=RELEASES.keys(),default='dr1',
                        help='data releases')
    parser.add_argument('--style',choices=['cat','hpx'],default='cat',
                        help='catalog style')
    parser.add_argument('-v','--verbose',action='store_true',
                        help='verbosity')
    parser.add_argument('-o','--outfile',
                        help='output filename (.fits,.csv,...)')
    args = parser.parse_args()

    if args.verbose:
        print("Selecting catalogs from release: %s"%args.release)
        print("RA, DEC, RADIUS = (%.5f deg, %.5f deg, %.1f arcsec)"%(args.ra,args.dec,args.radius))

    ra,dec = args.ra,args.dec # deg
    radius = args.radius/3600. # arcsec to deg

    catalog = load_catalogs(ra,dec,radius,release=args.release,
                            style=args.style,verbose=args.verbose)
    print("Catalog with %s objects."%len(catalog))

    if args.outfile:
        write_catalog(args.outfile,catalog)
    else:
        df = pd.DataFrame(catalog.byteswap().newbyteorder())
        if args.style == 'cat':
            columns=['RA','DEC','MAG_AUTO_G','MAG_AUTO_R','MAG_AUTO_I','MAG_AUTO_Z']
        elif args.style == 'hpx':
            columns = ['RA','DEC','MAG_AUTO','MJD_OBS']
        print(df[columns])
