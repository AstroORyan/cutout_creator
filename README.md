# cutout_creator
## Author: David O'Ryan
## Date: 31/01/2022

I plan for this algorithm to be a very basic galaxy cutout/thumbnail creator. The input will be a .csv of Right Ascension, Declination and (optionally) the galaxy name. The
algorithm will then use the RA and DEC to call the CAS Server of SDSS ahnd download the required FITs files. With these FITs files, it will then pull out the g, r and i colour
bands. These will then be used in the Lupton_rgb function of astropy to create a small cutout. It will then delete the full fits files from the memory/hard drive, and only leave
the cutout image(s).

So, to make this I will need to:

    1. Learn how to call the SDSS servers remotely from a function and pull down the relevant fits files.
    2. Learn how to tell if my cutouts have lots of 'blank' pixels (i.e. my galaxy is right on the edge of the image).
    3. Be able to stitch the images together in a relaible fashion.

Will update this README appropriately with development of the code.