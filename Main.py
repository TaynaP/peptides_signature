"""Idée 1 : on prend en entrée des peptides de RPG au format fasta, et on compare pour chaque peptide les peptides identiques
fichier txt en sortie"""
import sys
from combinaisons import *
from unique_peptides import *


def main(argv):
    pep = parse_csv(argv[0], 1)
    dict_p = compare_peptide(pep)

    # On crée le fichier contenant pour chaque peptides les peptides qui lui sont identiques
    pretty_print_all(dict_p, argv[1])

    # On crée le fichier contenant les peptides uniques
    u = unique_peptide(dict_p)
    pretty_print_unique_peptide(u, argv[2], argv[4])

    # Pour les séquences sans peptides uniques, on va chercher une combinaison (pour la séquence, ou le genre/la famille)
    mainSearchCombinations(dict_p, u, argv[3], argv[4])

    print("Done")


if __name__ == "__main__":
    main(sys.argv[1:])
