# Select DELVE Catalog Data at Fermilab

Select catalog data from a specific region of the DELVE catalog at Fermilab.

## Installation

Just git clone.

## Usage

```
usage: delve_select.py [-h] [-r RADIUS] [--release {dr1,dr2}]
                       [--style {cat,hpx}] [-v] [-o OUTFILE]
                       ra dec

Select objects from DELVE catalogs at Fermilab.

positional arguments:
  ra                    Right ascension (deg)
  dec                   Declination (deg)

optional arguments:
  -h, --help            show this help message and exit
  -r RADIUS, --radius RADIUS
                        selection radius (arcsec)
  --release {dr1,dr2}   data releases
  --style {cat,hpx}     catalog style
  -v, --verbose         verbosity
  -o OUTFILE, --outfile OUTFILE
                        output filename (.fits,.csv,...)
```

So for example:
```
./delve_select.py 170 -20 -r 360 --release dr2 --style cat -v -o out.fits
Selecting catalogs from release: dr2
RA, DEC, RADIUS = (170.00000 deg, -20.00000 deg, 360.0 arcsec)
Catalog filenames:
  /data/des91.b/data/kadrlica/projects/delve/cat/dr2/cat/cat_hpx_08188.fits
Loading cat_hpx_08188.fits...
Catalog with 1457 objects.
Writing out.fits...
```
