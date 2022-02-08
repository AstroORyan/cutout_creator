'''
author: oryan

This script will use the new SciServer package that I've just installed to download the cutouts from SkySerer from a list of RA and DECs. Let's go!!
'''

# Imports
import typer
from pathlib import Path
import pandas as pd
from SciServer import SkyServer
from astropy import units as u
from astropy.coordinates import SkyCoord
from PIL import Image
import uuid

# Functions

# Main Function
def main(
    manifest: Path = typer.Option(..., '-M',help='Please specify the path to the manifest of objects to get cutouts of. Requires columns of RA and DEC.'),
    save_folder: Path = typer.Option(...,'-S',help='Please specify the save folder for the cutouts you are downloading.'),
    dataRelease: str = typer.Option('DR14', '-D',help='Please specify the SDSS release you want to get the cutouts from.'),
    opt: str = typer.Option('','-O',help='Please specify any options you want for the cutout (see https://rdrr.io/github/sciserver/SciScript-R/man/SkyServer.getJpegImgCutout.html) for all options. Default: blank.'),
    query: str = typer.Option('','-Q',help = 'Please specify any extra SQL commands you want. Default: blank.')
    ):

    df = pd.read_csv(manifest,index_col=0)

    output_dict = {'OBJID':[],'filepath':[]}

    for i in range(len(df)):
        ra = df.RA.iloc[i]
        dec = df.DEC.iloc[i]
        
        img_arr = SkyServer.getJpegImgCutout(ra=ra, dec=dec, width=512,height=512,scale=0.4,dataRelease=dataRelease,opt=opt,query=query)

        img = Image.fromarray(img_arr)

        img.thumbnail((224,224))

        filename = f'{df.OBJID.iloc[i]}.png'
        save_path = Path.joinpath(save_folder,filename)

        img.save(save_path,quality=95)

        output_dict['OBJID'].append(df.OBJID.iloc[i])
        output_dict['filepath'].append(save_path)

        del ra, dec, img
    
    output_df = pd.DataFrame(output_dict)

    export_df = df.merge(output_df,on='OBJID',how='left')

    manifest_save = Path.joinpath(save_folder,'combined-sample-manifest.csv')
    export_df.to_csv(manifest_save)



# Initialisation
if __name__ == '__main__':
    typer.run(main)