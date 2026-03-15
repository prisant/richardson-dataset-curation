#This is a piece of sample code for iterating over the directory structure in
#the top 2018 dataset
#You will have to specify correct paths for your machine.

import os
import sys
import subprocess
import urllib.request
from pathlib import Path

def download_pdb(
        pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: 
    :the full path to the downloaded PDB file or 
    :None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None


if __name__ == "__main__":

  #path to top level dir of pdb storage
  top_pdb_dir = 'top2018_pdbs_full_filtered_hom70'

  #chain_file_path = "top2018_chains_hom70_mcfilter_60pct_complete.txt"
  chain_file_path = "top2018_chains_hom70_fullfiltered_60pct_complete.txt"
  with open( chain_file_path ) as chain_file:
      pdbchain_linelist = chain_file.readlines()

  #each line is one chain to look up
  #19HC_A
  #1A2Z_C
  #1A4I_B

  for idx, pdbchain_line in enumerate(pdbchain_linelist):

    #Uncomment for testing
    #if idx > 0:
        #break
 
    full = pdbchain_line.strip()
    pdb, subdir, chain = full[0:4], full[0:2], full[5]
    pdb_dir = os.path.join( top_pdb_dir, subdir, pdb )  
    pdb_pruned_filepath = os.path.join( pdb_dir, full+"_pruned_all.pdb" )
    #e.g. 'top2018_pdbs_full_filtered_hom70/1a/1a2z/1a2z_C_pruned_all.pdb'
    pdb_origin_filepath =  os.path.join( pdb_dir, pdb+".pdb" )
    #e.g. 'top2018_pdbs_full_filtered_hom70/1a/1a2z/1a2z.pdb'
    pdb_dssp_filepath = os.path.join( pdb_dir, pdb+".dssp" )
    
    if not Path(pdb_origin_filepath).is_file():
        try:
            funrtn = download_pdb(pdb,pdb_dir)
        except:
            #...
            print( "!", pdb, full[-1:], pdb_origin_filepath )
        else:
            #...
            print( "+", pdb, full[-1:], pdb_origin_filepath )
    else:
        #...
        print( "-", pdb, full[-1:], pdb_origin_filepath )

    if (
        not Path(pdb_dssp_filepath).is_file() and
        Path(pdb_origin_filepath).is_file()
    ):
        subprocess.run(
            [ "dssp", "-i", pdb_origin_filepath, "-o", pdb_dssp_filepath ]
        )
        print( "*", pdb, full[-1:], pdb_dssp_filepath )

