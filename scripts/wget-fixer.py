'''
author: David O'Ryan

This script will be used to take in the download.txt document and to update it with the correct spacing to be used in wget. Will take in a .txt document provided by the user and output a file named "TXT_File_Name_transposed.txt".

'''
# Imports
import typer
from pathlib import Path
# Functions

# Main Script
def main(
    df_path: Path = typer.Option(...,'-I',help='Please provide the full path to the wget file from SAS.')
):
    pass

# Initialisation
if __name__ == '__name__':
    typer.run(main)