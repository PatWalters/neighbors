#!/usr/bin/env python

import sys
import os
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.DataStructs import BulkTanimotoSimilarity
from operator import itemgetter
from docopt import docopt
from tqdm import tqdm
import csv

doc = """Usage:
neighbors.py --query QUERY_FILE --db DATABASE_FILE --out OUTPUT_CSV [--sim SIM_CUTOFF --max MAX_NEIGHBORS]

Options:
--query QUERY_FILE query file name
--db DATABASE_FILE database file name
--out OUTPUT_FILE output csv file
--sim SIMILARITY_CUTOFF minimum similarity cutoff for output (default 0.8) 
--max MAX_NEIGHBORS maximum number of neighbors to output (default 5)
"""


def molecule_supplier_from_name(input_file_name):
    ext = os.path.splitext(input_file_name)[-1]
    if ext == ".smi":
        suppl = Chem.SmilesMolSupplier(input_file_name, titleLine=False)
    elif ext == ".sdf":
        suppl = Chem.SDMolSupplier(input_file_name)
    else:
        print("%s is not a valid molecule extension" % ext)
        sys.exit(0)
    return suppl


def fingerprints_from_file(infile_name):
    suppl = molecule_supplier_from_name(infile_name)
    fp_list = []
    for mol in tqdm(suppl,desc=f"Reading {infile_name}"):
        if mol:
            fp_list.append([Chem.MolToSmiles(mol), mol.GetProp("_Name"), MACCSkeys.GenMACCSKeys(mol)])
    return fp_list


def get_neighbors(query_fp, fp_list, cutoff=0.8, max_nbrs=5):
    sim_list = BulkTanimotoSimilarity(query_fp, [x[2] for x in fp_list])
    nbr_list = [x[0] + [x[1]] for x in zip(fp_list, sim_list) if x[1] >= cutoff]
    nbr_list.sort(key=itemgetter(3),reverse=True)
    nbr_list = nbr_list[0:max_nbrs]
    return nbr_list


if __name__ == "__main__":
    options = docopt(doc)
    query_filename = options.get("--query")
    db_filename = options.get("--db")
    output_filename = options.get("--out")
    max_neighbors = options.get("--max") or 5
    max_neighbors = int(max_neighbors)
    sim_limit = options.get("--sim") or 0.8
    sim_limit = float(sim_limit)

    db_fp_list = fingerprints_from_file(db_filename)
    query_fp_list = fingerprints_from_file(query_filename)

    writer = csv.writer(open(output_filename, "w"))
    writer.writerow(["SMILES", "Name", "SIMILARITY", "QUERY"])
    for query_smi, query_name, query_fp in tqdm(query_fp_list, desc="Finding Neighbors"):
        res = get_neighbors(query_fp, db_fp_list, cutoff=sim_limit,max_nbrs=max_neighbors)
        for row in res:
            smiles, name, _, sim = row
            writer.writerow([smiles, name, f"{sim:.2f}", query_name])
