#!/usr/bin/env python3

import sys
from combinaisons import *
from unique_peptides import *
import argparse

os.chmod("Main.py", 0o755)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("peptidesCSVFile", help="The name of the file containing peptides (CSV File)", type=str)
    parser.add_argument("resultsPath", help="The directory (path) wich will contain all the result files (ex : results/), the path must end with a '/'", type=str)
    parser.add_argument("-I", "--peptideThreshold", help="The minimum length of the peptides we keep (default 1)", type=int, default=1)
    args = parser.parse_args()

    pep = parse_csv(args.peptidesCSVFile, args.peptideThreshold)
    dict_p = compare_peptide(pep)

    # On crée le fichier contenant pour chaque peptides les peptides qui lui sont identiques
    pretty_print_all(dict_p, args.resultsPath)

    # On crée le fichier contenant les peptides uniques
    u = unique_peptide(dict_p)
    pretty_print_unique_peptide(u, args.resultsPath)

    # Pour les séquences sans peptides uniques, on va chercher une combinaison (pour la séquence, ou le genre/la famille)
    mainSearchCombinations(dict_p, u, args.resultsPath)

    print("Done")
    
