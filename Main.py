#!/usr/bin/env python3

import sys
from combinaisons import *
from unique_peptides import *
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--peptidesCSVFile", required=True, help="The name of the file containing peptides (CSV File)", type=str)
    parser.add_argument("-o", "--resultsPath", required=True, help="The directory (path) wich will contain all the result files (ex : results/)", type=str)
    parser.add_argument("-p", "--peptideThreshold", help="The minimum length of the peptides we keep (default 1)", type=int, default=1)
    parser.add_argument("-c", "--combinationThreshold", help="The maximum size of a combination (default no limit)", type=int, default=0)
    args = parser.parse_args()

    if not os.path.exists(args.resultsPath):
        os.makedirs(args.resultsPath)
    else:
        raise Exception("This directory already exists. Please give another directory.")

    if args.resultsPath[-1] != '/':
        args.resultsPath += '/'
    if args.resultsPath[0] == '/':
        args.resultsPath = args.resultsPath[1:]

    pep = parse_csv(args.peptidesCSVFile, args.peptideThreshold)
    dict_p = compare_peptide(pep)

    # On crée le fichier contenant pour chaque peptides les peptides qui lui sont identiques
    pretty_print_all(dict_p, args.resultsPath)

    # On crée le fichier contenant les peptides uniques
    u = unique_peptide(dict_p)
    pretty_print_unique_peptide(u, args.resultsPath)

    # Pour les séquences sans peptides uniques, on va chercher une combinaison (pour la séquence, ou le genre/la famille)
    mainSearchCombinations(dict_p, u, args.resultsPath, args.combinationThreshold)

    print("Done")
    
