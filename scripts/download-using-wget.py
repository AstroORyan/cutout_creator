'''
Author: David O'Ryan
Date: 13/02/2022
Point of this script is to use wget to download observations and fits files, but without have to re-run it tons of times whenever it crashes. Should intake the .txt file downloaded from SAS.
'''
import typer
import wget
import sys
import pandas as pd
import os
import time
import glob
from tqdm import tqdm

def main(
    files_path: str = typer.Option(...,'-P',help='Please specify the path to the downloaded file paths from SAS. Should be a .txt file.'),
    save_path: str = typer.Option(...,'-S',help='Please specify the folder you wish to save everything in.')
):
    manifest = pd.read_csv(files_path,delimiter = ' ')
    if len(manifest) < 1.5:
        manifest_df = manifest.T.reset_index()
        manifest_df = manifest_df.rename(columns={'index':'filepaths'})
    else:
        manifest_df = manifest.reset_index()
        manifest_df = manifest_df.rename(columns={'index':'filepaths'})
    
    manifest_df_filenames = (
        manifest_df
        .assign(filenames = manifest_df.filepaths.apply(lambda x: os.path.basename(x)))
    )

    existing_filepaths = glob.glob(save_path + '*.fits.bz2')
    existing_files = []
    for i in existing_filepaths:
        existing_files.append(os.path.basename(i))

    for i in tqdm(range(len(manifest_df_filenames))):
        row = manifest_df_filenames.iloc[i]
        SDSS_path = row.filepaths
        local_save_filename = row.filenames

        if os.path.basename(SDSS_path) in existing_files:
            continue
        
        downloaded = False
        attempts = 0
        while not downloaded and attempts < 10:
            try:
                wget.download(SDSS_path,out=save_path+local_save_filename)
                downloaded = True
            except:
                print(f'Error in {attempts} download attempt. {10 - attempts} attempts left. Pausing download for 3 seconds...')
                time.sleep(3)
                attempts += 1

        if not downloaded:
            print('WARNING: Attempted to connect to SDSS 10 times and still failed. Aborting for safety. Check url causing issue.')
            print(f'Errored url = {SDSS_path}')
            sys.exit()



if __name__ == '__main__':
    typer.run(main)