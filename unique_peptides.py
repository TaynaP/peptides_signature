from peptide_v2 import Peptides
import itertools
import csv
import os

def parse_csv(file, threshold=1):
    """Lit le fichier multifasta de RPG et crée une liste d'objet Peptides avec tout les peptides issus de la digestion

    Args:
        file (str): the output of RPG in csv

    Raises:
        NameError: if the extension of the file is not "csv"

    Returns:
        list: list of all the peptides from RPG
    """
    _, ext = os.path.splitext(file)
    if ext != ".csv":
        raise NameError("the input file must be a csv file")
    list_pep = []
    nb_prot = 0
    current_prot_name = ""
    f= open (file)
    myReader = csv.reader(f)
    next(myReader)
    for row in myReader:
        if int(row[4]) >= threshold:
            id = row[0].split('_')
            family = id[0]
            genus = id[1]
            prot_name = ' '.join(id[2:])
            if prot_name != current_prot_name :
                nb_prot += 1
                current_prot_name = prot_name
            nb_peptide = int(row[1]) + 1
            list_pep.append(Peptides(row[7], int(row[1]) + 1, prot_name, nb_prot, genus, family, row[3]))
    return list_pep


def compare_peptide(liste):
    """Permet de comparer les séquences des peptides et stocker les peptides dans un dictionnaire.
    keys : les peptides
    values : liste avec les peptides identiques

    Args:
        liste (List): the output of parse_csv()

    Raises:
        TypeError: if the parameter is not a list

    Returns:
        Dict: the dictionary with the peptides as keys the peptides that shares the same sequence as values
    """
    if not isinstance(liste, list):
        raise TypeError("the input file must be a list of peptides")
    pep_dict = {}
    for pep_a, pep_b in itertools.combinations(liste,
                                               2):  # combinations() permet d'avoir des couples de peptides uniques donc si j'ai (peptide1, peptide2) je n'aurais pas (peptide2, peptide1)
        if pep_a.get_seq() == pep_b.get_seq():
            if pep_a not in pep_dict:
                pep_dict[pep_a] = [pep_b]
            else:
                pep_dict[pep_a].append(pep_b)
            if pep_b not in pep_dict:
                pep_dict[pep_b] = [pep_a]
            else:
                pep_dict[pep_b].append(pep_a)
        else:
            if pep_a not in pep_dict:
                pep_dict[pep_a] = []
            if pep_b not in pep_dict:
                pep_dict[pep_b] = []

    return dict(sorted(pep_dict.items(), key=lambda item: (item[0].get_nb_prot(), item[
        0].get_nb_peptide())))  # ici c'est pour avoir un dico dans l'ordre des protéines et par ordre de peptides


def pretty_print_all(dico, output_dir):
    """Permet de formatage du fichier txt pour tous les peptides.

    Args:
        dico (Dict): the output of compare_peptide()
        output_file (str): the name of the output file we want

    Raises:
        TypeError: if the parameters is not a dict and a str
    """
    currentSequence = ""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + 'liste_peptides.txt', 'w') as results:
        peptides = list(dico.keys())
        for pep in peptides:
            if pep.get_prot_name() != currentSequence:
                results.write("Protéine {} : {}\n".format(str(pep.get_nb_prot()), pep.get_prot_name()))
                currentSequence = pep.get_prot_name()
            results.write("\tPeptide {} : {}\n".format(str(pep.get_nb_peptide()), '|'.join(
                map(lambda s: str(s).strip('()').replace(' ', ''),
                    [(elt.get_nb_prot(), elt.get_nb_peptide()) for elt in dico[pep]]))))


def unique_peptide(dico):
    """Retoune les peptides uniques dans une liste

    Args:
        dico (Dict): the output of compare_peptide()

    Raises:
        TypeError: if the parameter is not a dictionary
        Exception: if there are no unique peptides

    Returns:
        List: the list containing only the unique peptides
    """
    if not isinstance(dico, dict):
        raise TypeError("The parameter must be a dictionary")
    unique_peptide_list = []
    for pep, pep_id_list in dico.items():
        if pep_id_list == []:
            unique_peptide_list.append(pep)
        else:
            unique = True
            for peptide in pep_id_list:
                if peptide.get_nb_prot() != pep.get_nb_prot():
                    unique = False
                    break
            if unique:
                unique_peptide_list.append(pep)
    if unique_peptide_list == []:
        raise Exception("There are no unique peptides")
    return unique_peptide_list


def pretty_print_unique_peptide(liste, output_dir):
    """Permet le formatage du fichier txt seulement pour les peptides uniques pour chaque protéine

    Args:
        liste (List): the output of unique_peptide()
        output_file (str): the name of the output file we want
        allResultsFile (str) : the name of the output file regrouping the results

    Raises:
        TypeError: if the parameters is not a list and a str
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + 'uniquePeptides.csv', 'w', newline='') as results:
        with open(output_dir + 'allResults.csv', 'w', newline='') as allRes:
            writer_all = csv.writer(allRes, delimiter='|')
            writer_unique = csv.writer(results, delimiter='|')
            writer_all.writerow(["Unique peptides for each sequence :"])
            writer_all.writerow(["Family", "Genus", "Name", "Position", "Peptide seq"])
            writer_unique.writerow(["Family", "Genus", "Name", "Position", "Peptide seq"])
            for peptide in liste:
                rowToInsert = [peptide.get_family(), peptide.get_genus(), peptide.get_prot_name(), peptide.get_position(), peptide.get_seq()]
                writer_all.writerow(rowToInsert)
                writer_unique.writerow(rowToInsert) 
