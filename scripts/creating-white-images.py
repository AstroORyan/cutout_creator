'''
Author: David O'Ryan (10.02.2022)
This script will be responsible for creating the white images that will be used as my observational dataset. It 
will read in my positions, and go through all the fits files in the observations folder and pick out the ones
which contain the target.

If there are multiple fits files, they will be mosaicked together, and then a cutout applied. Might need some user
input which allows me to check whether the cutout is good enough. 
'''
# Imports
import typer
import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import sys
from tqdm import tqdm
import math

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy import wcs
from astropy.nddata import Cutout2D, block_reduce
import astropy.units as u
from astropy.coordinates import SkyCoord

from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
from reproject import reproject_interp

from regions import PolygonSkyRegion


# Functions

# Main Function
def main(
    manifest_path: str = typer.Option(...,'-M',help='Please specify the path to the target manifest containing their RAs and DECs. Must be a .csv containing column headers of RA and DEC.'),
    fits_folder: str = typer.Option(...,'-F',help='Please specify the path to the folder containing all needed fits folders.'),
    save_path: str = typer.Option(...,'-S',help='Please specify a folder path to save the example white images and data in.')
):
    print('Loading in Data and Initialising...')
    # Download the RA and DEC list with galaxy Name:
    manifest_df = pd.read_csv(manifest_path)[['Names','Prim_RA','Prim_DEC','Sec_RA','Sec_DEC']]

    # Getting the fits file paths:
    fits_list = glob.glob(fits_folder + '*.fits.bz2')
    #fits_list = fits_list[:10]

    # Initialising dictionaries:
    coords_dict = {}
    sorted_dict = {}

    print('Initialisation Complete.')

    print('Loading in all file WCS...')
    # Finding which fits files go with each galaxy.
    for i in tqdm(fits_list):
        w = wcs.WCS(fits.open(i)[0].header)
        filename = os.path.basename(i)
        coords_dict[str(filename)] = w
    print('WCS extracted.')

    print('Assigning files to galaxy...')
    for i in tqdm(range(len(manifest_df))):
        gal_name = manifest_df.Names.iloc[i]
        prim_ra = manifest_df.Prim_RA.iloc[i]
        prim_dec = manifest_df.Prim_DEC.iloc[i]
        sec_ra = manifest_df.Sec_RA.iloc[i]
        sec_dec = manifest_df.Sec_DEC.iloc[i]

        sorted_dict[str(gal_name)] = []

        target_coords_prim = SkyCoord(ra=prim_ra * u.degree, dec=prim_dec * u.degree, frame = 'icrs')
        target_coords_sec = SkyCoord(ra=sec_ra * u.degree, dec=sec_dec * u.degree, frame = 'icrs')

        for j in coords_dict.keys():
            w = coords_dict[j]
            vertices = SkyCoord(w.calc_footprint(),unit='deg',frame='icrs')
            region_sky = PolygonSkyRegion(vertices=vertices)
            if region_sky.contains(target_coords_prim,w)[0] or region_sky.contains(target_coords_sec,w)[0]:
                sorted_dict[str(gal_name)].append(j)
        
        if len(sorted_dict[str(gal_name)]) == 0:
            sorted_dict[str(gal_name)].append('Pass')
    print('Completed.')

    # Creating the White Images:
    print('Building white images (this may take a while if a lot of mosaicking must occur)...')
    reduction_dict = {}
    for i in tqdm(sorted_dict.keys()):
        fits_files = sorted_dict[i]
        white_image_list = []
        if fits_files[0] == 'Pass':
            continue
        elif len(fits_files) > 5:
            colour_dict = {'u':[],'g':[],'r':[],'i':[],'z':[]}
            for j in fits_files:
                colour_dict[j.split('-')[1]].append(j)
            for j in colour_dict.keys():
                galaxy_hdus = [fits.open(fits_folder + file)[0] for file in colour_dict[j]]
                wcs_out,shape_out = find_optimal_celestial_wcs(galaxy_hdus)
                array,footprint = reproject_and_coadd(galaxy_hdus,wcs_out,shape_out=shape_out,reproject_function=reproject_interp)
                position = SkyCoord(manifest_df.query('Names == @i').Prim_RA * u.degree, manifest_df.query('Names == @i').Prim_DEC * u.degree, frame = 'icrs')
                position_sec = SkyCoord(manifest_df.query('Names == @i').Sec_RA * u.degree, manifest_df.query('Names == @i').Sec_DEC * u.degree, frame = 'icrs')

                position_pix = position.to_pixel(wcs_out)
                position_sec_pix = position_sec.to_pixel(wcs_out)

                dist = 2.5 * np.sqrt((np.float64(position_pix[0]) - np.float64(position_sec_pix[0]))**2 + (np.float64(position_pix[1]) - np.float64(position_sec_pix[1]))**2)

                cutoff = int(math.ceil(dist / 100)) * 100

                tmp = Cutout2D(array,position,(cutoff,cutoff),wcs=wcs_out,mode='partial').data * 1e-26 * 1e-3
                white_image_list.append(tmp)
        else:
            for j in range(5):
                data = fits.open(fits_folder + fits_files[j])[0].data
                w = coords_dict[fits_files[j]]
                position = SkyCoord(manifest_df.query('Names == @i').Prim_RA * u.degree, manifest_df.query('Names == @i').Prim_DEC * u.degree, frame = 'icrs')
                position_sec = SkyCoord(manifest_df.query('Names == @i').Sec_RA * u.degree, manifest_df.query('Names == @i').Sec_DEC * u.degree, frame = 'icrs')

                position_pix = position.to_pixel(w)
                position_sec_pix = position_sec.to_pixel(w)

                dist = 2.5 * np.sqrt((np.float64(position_pix[0]) - np.float64(position_sec_pix[0]))**2 + (np.float64(position_pix[1]) - np.float64(position_sec_pix[1]))**2)

                cutoff = int(math.ceil(dist / 100)) * 100

                tmp = Cutout2D(data,position,(cutoff,cutoff),wcs=w,mode='partial').data * 1e-26 * 1e-3
                white_image_list.append(tmp)

        if len(white_image_list) > 5:
            print('Something went very wrong in white image loop!!')

        white_image = np.zeros(white_image_list[0].shape)
        
        for p in range(len(white_image_list)):
            white_image += white_image_list[p]

        reduction = white_image.shape[0]/50

        reduction_dict[i] = reduction
        
        white_image_reduced = block_reduce(white_image, reduction)

        np.save(f'{save_path}{i}.npy',white_image_reduced)
        plt.figure(figsize=(12,8))
        plt.imshow(-2.5*np.log10(white_image) - 48.6)
        plt.colorbar()
        plt.title(f'{i} with block_reduce = {reduction}')
        plt.savefig(save_path+i+'.jpg')
        plt.close()
    
    reduction_df = (
        pd.DataFrame([reduction_dict])
        .T
        .reset_index()
        .rename(columns={'index':'Names',0:'block_reduce'})
    )

    reduction_df.to_csv(save_path+'block_reductions.csv')

    print('Images created.')
    print('Algorithm Complete.')

    
# Initialisation
if __name__ == '__main__':
    typer.run(main)
