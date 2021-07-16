import itertools
import combinaisons
import os
import csv

def getAllGenus(peptideToProtein):
    """Crée la liste de toutes les genres

            Args:
                peptideToProtein (dict) : les peptides associés à leur protéine

            Raises:
                /

            Returns:
                list : une liste contenant toutes les genres des séquences dans le fichier fasta de départ
            """
    res = []
    for prot, peps in peptideToProtein.items():
        for elt in peps:
            if elt.get_genus() not in res:
                res.append(elt.get_genus())
    return res


def assignPeptideToGenus(peptideToProtein):
    """Crée un dictionnaire qui associe à chaque genre les peptides qui y apparaissent

            Args:
                peptideToProtein (dict) : les peptides associés à leur protéine

            Raises:
                /

            Returns:
                dict : dictionnaire contenant en clé les genres, et en valeurs les peptides appartenant à ces genres
    """
    res = {}
    for prot, otherPeps in peptideToProtein.items():
        for otherP in otherPeps:
            if otherP.get_genus() not in res:
                res[otherP.get_genus()] = [otherP]
            else:
                if otherP not in res[otherP.get_genus()]:
                    res[otherP.get_genus()].append(otherP)
    return res


def seqToGenus(peptideToProtein):
    """Crée un dictionnaire qui associe à chaque genre les séquences leur appartenant

            Args:
                peptideToProtein (dict) : les peptides associés à leur protéine

            Raises:
                /

            Returns:
                dict : dictionnaire qui associe à chaque genre les séquences leur appartenant
     """
    res = {}
    for prot, peptides in peptideToProtein.items():
        if peptides[0].get_genus() not in res:
            res[peptides[0].get_genus()] = [prot]
        else:
            res[peptides[0].get_genus()].append(prot)
    return res


def getGenusWithoutUnique(peptideToGenus, uniquePeptides):
    """Crée une liste contenant les genre n'ayant pas de peptide unique

    Args:
        peptideToGenus (dict) : le dictionnaire associant les peptides à leur protéine
        uniquePeptides (list) : la liste des peptides uniques

    Raises:
        /

    Returns:
        list : une liste contenant les genres n'ayant pas de peptide unique
    """
    res = list(peptideToGenus.keys())
    for pep in uniquePeptides:
        if pep.get_genus() in res:
            res.remove(pep.get_genus())
    return res


def genusOfpeptides(dico):
    dico_genus = {}
    for pep, pep_list in dico.items():
        genus = []
        genus.append(pep.get_genus())
        for pep_identique in pep_list:
            if pep_identique.get_genus() != pep.get_genus():
                if pep_identique.get_genus() not in genus:
                    genus.append(pep_identique.get_genus())
        dico_genus[pep] = genus
    return dico_genus


def where_pep_present_genre(dico):
    """stock the different species the peptide appears in and make sublist of species that are part of the same genus.

    Args:
        dico (dict): the dictonnary made after compare_peptide()

    Returns:
        dict: key = peptide
              values = list of sublist of species the peptide appears in (sublist = genus)
    """
    dico_genre = {}
    for pep, pep_list in dico.items():  # On parcours le dico
        species = []  # liste pour stocker les espèces où le peptide apparaît
        species.append(pep)
        for pep_identique in pep_list:
            if pep_identique.get_nb_prot() != pep.get_nb_prot():
                species.append(pep_identique)
        genus_list = []  # liste pour pouvoir constituer les sous listes d'espèces
        species = sorted(species, key=lambda
            pep: pep.get_genus())  # on range les espèces en fonction du genre auquelles elles appartiennent
        for k, g in itertools.groupby(species,
                                      lambda pep: pep.get_genus()):  # fonction du module itertools (permet de regrouper en sous listes les elements d'une liste de départ en fonction d'un critère, ici le genre)
            genus_list.append(list(g))
        dico_genre[
            pep] = genus_list  # le dico avec pour chaque peptide, les espèces regroupées en genre auquel il apparaît
    return dico_genre

def unique_pep_genre(dico_pep_genre):
    """Let the user find the peptides that are unique to only one genus (if the peptide appears only in one genus).

    Args:
        dico_pep_genre (dict): the dictonnary made in where_pep_present_genre()

    Returns:
        list: the list containing the peptides that are unique for one genus
    """
    unique_pep_genre_list = []
    for pep, genre_list in dico_pep_genre.items():
        unique = True
        for i in range(len(genre_list)):
            if i > 0:  # si la liste contient plus d'une sous liste, c'est que le peptide apparaît dans plus d'un genre donc il n'est pas unique à un genre
                unique = False
                break
        if unique:
            unique_pep_genre_list.append(pep)
    return unique_pep_genre_list


def pretty_print_unique_peptide_genus(liste, output_dir):
    """Permet le formatage du fichier txt seulement pour les peptides uniques pour chaque genre

    Args:
        liste (List): the output of unique_peptide()
        output_file (str): the name of the output file we want
        allResultsFile (str) : the name of the output file that contains all the results

    Raises:
        TypeError: if the parameters is not a list and a str
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + 'unique_pep_genre.csv', 'w', newline='') as results:
        writer_genus = csv.writer(results)
        writer_genus.writerow(["Family", "Genus", "Name of prot", "Position", "Peptide mass", "Peptide seq"])
        for peptide in liste:
            rowToInsert = [peptide.get_family(), peptide.get_genus(), peptide.get_prot_name(), peptide.get_position(), peptide.get_mass(), peptide.get_seq()]
            writer_genus.writerow(rowToInsert)


def mainGenre(dict_p, output_dir, peptidesToProtein):
    dict_g = where_pep_present_genre(dict_p)
    uniquePepGenre = unique_pep_genre(dict_g)

    # Création du fichier contenant les peptides uniques pour le genre
    pretty_print_unique_peptide_genus(uniquePepGenre, output_dir)

    # Liste des séquences n'ayant pas de peptides uniques pour un genre
    seqWithoutUniqueGenre = combinaisons.getSequencesWithoutUnique(peptidesToProtein, uniquePepGenre)
    return seqWithoutUniqueGenre
