#!/usr/bin/env python

import glob
import os

import pandas as pd
from tqdm import tqdm

from ppmp.protein import Protein

JSON_PATH = './data/json/'
PDB_PATH = './data/pdb/'
CSV_PATH = './data/csv/'

protein_name = list()


# Read all protein names form the filesystem
for path in tqdm(glob.glob(JSON_PATH + '*.json'), desc='Pre-loading protein IDs'):
    file_name = os.path.basename(path)
    protein_name.append(os.path.splitext(file_name)[0])

for name in tqdm(protein_name, desc='Calculating RMSD'):
    original = Protein(name,
                       PDB_PATH + name + '.pdb',
                       JSON_PATH + name + '.json')

    protein_df = pd.DataFrame()

    for path in glob.glob('{}{}_*.pdb'.format(PDB_PATH, name)):
        variation = os.path.splitext(os.path.basename(path))[0]
        perturbed = Protein(variation,
                            path,
                            JSON_PATH + name + '.json')
        temp_df = pd.DataFrame()
        temp_df['module'] = original.modules_chain
        temp_df['variation'] = [variation for _ in original.modules_chain]
        temp_df['order'] = [i for i in range(len(original.modules_chain))]
        temp_df['rmsd'] = Protein.kabsch(original, perturbed, range(len(original.modules_chain)))
        protein_df = pd.concat([protein_df, temp_df])
    protein_df.to_csv(CSV_PATH + name + '.csv', index=False)
