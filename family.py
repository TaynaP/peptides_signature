import itertools
import combinaisons
import os
import csv

def getAllFamilies(peptideToProtein):
    """Crée la liste de toutes les familles

        Args:
            peptideToProtein (dict) : les peptides associés à leur protéine

        Raises:
            /

        Returns:
            list : une liste contenant toutes les familles des séquences dans le fichier fasta de départ
        """
    res = []
    for prot, peps in peptideToProtein.items():
        for elt in peps:
            if elt.get_family() not in res:
                res.append(elt.get_family())
    return res


def assignPeptideToFamily(peptideToProtein):
    """Crée un dictionnaire qui associe à chaque famille les peptides qui y apparaissent

        Args:
            peptideToProtein (dict) : les peptides associés à leur protéine

        Raises:
            /

        Returns:
            dict : dictionnaire contenant en clé les familles, et en valeurs les peptides appartenant à ces familles
        """
    res = {}
    for prot, otherPeps in peptideToProtein.items():
        for otherP in otherPeps:
            if otherP.get_family() not in res:
                res[otherP.get_family()] = [otherP]
            else:
                if otherP not in res[otherP.get_family()]:
                    res[otherP.get_family()].append(otherP)
    return res


def familyOfPeptides(dicoPeptides):
    """Crée un dictionnaire contenant pour chaque peptide les familles auxquelles il appartient

        Args:
            dicoPeptides (dict) : les dictionnaire contenant les peptides, et les peptides identiques à celui-ci en valeur

        Raises:
            /

        Returns:
            dict : dictionnaire avec les peptides en clé, et les familles auxquelles ils appartiennet en valeur
        """
    dico_family = {}
    for pep, pep_list in dicoPeptides.items():
        family = []
        family.append(pep.get_family())
        for pep_identique in pep_list:
            if pep_identique.get_family() != pep.get_family():
                if pep_identique.get_family() not in family:
                    family.append(pep_identique.get_family())
        dico_family[pep] = family
    return dico_family


def seqToFamily(peptideToProtein):
    """Crée un dictionnaire qui associe à chaque famille les séquences leur appartenant

        Args:
            peptideToProtein (dict) : les peptides associés à leur protéine

        Raises:
            /

        Returns:
            dict : dictionnaire qui associe à chaque famille les séquences leur appartenant
        """
    res = {}
    for prot, peptides in peptideToProtein.items():
        if peptides[0].get_family() not in res:
            res[peptides[0].get_family()] = [prot]
        else:
            res[peptides[0].get_family()].append(prot)
    return res


def getFamilyWithoutUnique(peptideToFamily, uniquePeptides):
    """Crée une liste contenant les familles n'ayant pas de peptide unique

    Args:
        peptideToFamily (dict) : le dictionnaire associant les peptides à leur famille
        uniquePeptides (list) : la liste des peptides uniques pour une séquence

    Raises:
        /

    Returns:
        list : une liste contenant les familles n'ayant pas de peptide unique
    """
    res = list(peptideToFamily.keys())
    for pep in uniquePeptides:
        if pep.get_family() in res:
            res.remove(pep.get_family())
    return res


def where_pep_present_family(dico):
    """stock the different species the peptide appears in and make sublist of species that are part of the same family.

    Args:
        dico (dict): the dictonnary made after compare_peptide()

    Returns:
        dict: key = peptide
              values = list of sublist of species the peptide appears in (sublist = family)
    """
    dico_family = {}
    for pep, pep_list in dico.items():  # On parcours le dico
        species = []  # liste pour stocker les espèces où le peptide apparaît
        species.append(pep)
        for pep_identique in pep_list:
            if pep_identique.get_nb_prot() != pep.get_nb_prot():
                species.append(pep_identique)
        family_list = []  # liste pour pouvoir constituer les sous listes d'espèces
        species = sorted(species, key=lambda
            pep: pep.get_family())  # on range les espèces en fonction de la famille auquelles elles appartiennent
        for k, g in itertools.groupby(species,
                                      lambda pep: pep.get_family()):  # fonction du module itertools (permet de regrouper en sous listes les elements d'une liste de départ en fonction d'un critère, ici la famille)
            family_list.append(list(g))
        dico_family[
            pep] = family_list  # le dico avec pour chaque peptide, les espèces regroupées en genre auquel il apparaît

    return dico_family


def unique_pep_family(dico_pep_family):
    """Let the user find the peptides that are unique to only one family (if the peptide appears only in one family).

    Args:
        dico_pep_family (dict): the dictonnary made in where_pep_present_family()

    Returns:
        list: the list containing the peptides that are unique for one family
    """
    unique_pep_family_list = []
    for pep, family_list in dico_pep_family.items():
        unique = True
        for i in range(len(family_list)):
            if i > 0:  # si la liste contient plus d'une sous liste, c'est que le peptide apparaît dans plus d'une famille donc il n'est pas unique à une famille
                unique = False
                break
        if unique:
            unique_pep_family_list.append(pep)
    return sorted(unique_pep_family_list, key=lambda pep: pep.get_family())


def pretty_print_unique_peptide_family(liste, listAllFamily, output_dir):
    """Permet le formatage du fichier txt seulement pour les peptides uniques pour chaque famille

    Args:
        liste (List): the output of unique_pep_family()
        output_file (str): the name of the output file we want
        allResultsFile (str) : the name of the output file that contains all the results

    Raises:
        TypeError: if the parameters is not a list and a str
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + 'unique_pep_family.csv', 'w', newline='') as results:
        writer_family = csv.writer(results)
        writer_family.writerow(["Family", "Genus", "Protein Name", "Position", "Peptide mass", "Peptide seq"])
        for peptide in liste:
            rowToInsert = [peptide.get_family(), peptide.get_genus(), peptide.get_prot_name(), peptide.get_position(), peptide.get_mass(), peptide.get_seq()]
            writer_family.writerow(rowToInsert)
        writer_family.writerow("")
        noUnique = []
        for family in listAllFamily:
            tmp = False
            for elt in liste :
                if family == elt.get_family():
                    tmp = True
            if not tmp:
                noUnique.append(family)
        strNoUnique = ",".join(noUnique)
        writer_family.writerow(["Families that don't have a unique peptide :" + strNoUnique])


def mainFamily(dict_p, output_dir, peptidesToProtein):
    dict_f = where_pep_present_family(dict_p)
    uniquePepFamily = unique_pep_family(dict_f)
    AllFamilies = getAllFamilies(peptidesToProtein)

    # Création du fichier contenant les peptides uniques pour chaque famille
    pretty_print_unique_peptide_family(uniquePepFamily, AllFamilies, output_dir)

    # On cherche et renvoie les familles sans peptides uniques
    seqWithoutUniqueFamily = combinaisons.getSequencesWithoutUnique(peptidesToProtein, uniquePepFamily)
    return seqWithoutUniqueFamily

