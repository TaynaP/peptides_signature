from family import *
from genre import *
import warnings
import os

# *********************************************************************************************************************
# FONCTIONS INTERMÉDIAIRES

def dicoNameNumberSequence(listePep):
    """Lit la liste des peptides et crée un dictionnaire liant un nom de séquence à son numéro

    Args:
        listePep (list) : la liste des peptides

    Raises:
        /

    Returns:
        dict : un dictionnaire liant un nom de séquence peptidique à son numéro
    """
    res = {}
    for pep in listePep:
        if pep.get_prot_name() not in res:
            theName = "{} {}".format(pep.get_genus(), pep.get_prot_name())
            res[theName] = pep.get_nb_prot()
    return res


def searchSpeciesForPeptides(dicoPeptides):
    """Lit le dictionnaire entier des peptides et en crée un nouveau indiquant pour chaque peptide les espèces dans lesquelles il apparaît

    Args:
        dicoPeptides (dict) : le dictionnaire de tous les peptides

    Raises:
        /

    Returns:
        dict : le dictionnaire contenant les clés/valeurs : peptide/espèces dans lesquelles il apparaît
    """
    res = {}
    for pep, listPep in dicoPeptides.items():
        if pep not in res:
            res[pep] = [pep.get_nb_prot()]
            for otherPeptide in listPep:
                if otherPeptide.get_nb_prot() not in res[pep]:
                    res[pep].append(otherPeptide.get_nb_prot())
                if otherPeptide not in res:
                    res[otherPeptide] = [otherPeptide.get_nb_prot(), pep.get_nb_prot()]
                    for otherInList in listPep:
                        if otherPeptide != otherInList:
                            if otherInList.get_nb_prot() not in res[otherPeptide]:
                                res[otherPeptide].append(otherInList.get_nb_prot())
    return dict(sorted(res.items(), key=lambda item: (item[0].get_nb_prot(), item[0].get_nb_peptide())))


def assignPeptideToProtein(speciesOfPeptides):
    """Lit le dictionnaires des peptides et leurs espèces, et en crée un en associant les peptides aux protéines

    Args:
        speciesOfPeptides (dict) : les peptides et les espèces dans lesquelles ils apparaissent

    Raises:
        /

    Returns:
        dict : les peptides et leur protéine associée
    """
    res = {}
    for pep, otherPeps in speciesOfPeptides.items():
        if pep.get_nb_prot() not in res:
            res[pep.get_nb_prot()] = [pep]
        else:
            res[pep.get_nb_prot()].append(pep)
    return res


def peptideInSequence(pep, prot, peptideToProtein):
    """Renvoie vrai le peptide appartient à la séquence en paramètre, faux sinon

    Args:
        pep (Peptide) : le peptide dont nous voulons savoir l'appartenance à une séquence
        prot (int) : le numéro de la séquence protéique dans laquelle nous voulons tester l'appartenance du peptide
        peptideToProtein (dict) : le dictionnaire associant les protéines à leurs peptides

    Raises:
        /

    Returns:
        boolean : vrai si le peptide appartient à la séquence, faux sinon
    """
    for peps in peptideToProtein[prot]:
        if pep.get_seq() == peps.get_seq():
            return True
    return False


def getSequencesForPeptides(peptides, peptideToProtein):
    """Lit le dictionnaires des peptides et leurs espèces, et en crée un en associant les peptides aux séquences dans lesquelles ils apparaissent

    Args:
        peptides (list) : les peptides pour lesquels nous voulons la liste des séquences auxquels ils appartiennent
        peptideToProtein (dict) : le dictionnaire associant les protéines à leurs peptides

    Raises:
        /

    Returns:
        list : liste des séquences contenant les peptides donnés en paramètre
    """
    res = []
    for prot, peps in peptideToProtein.items():
        tmp = True
        for peptide in peptides:
            if not peptideInSequence(peptide, prot, peptideToProtein):
                tmp = False
        if tmp:
            res.append(prot)
    return res


def getSequencesWithoutUnique(peptideToProt, uniquePeptides):
    """Crée une liste contenant les séquences n'ayant pas de peptide unique

    Args:
        peptideToProt (dict) : le dictionnaire associant les peptides à leur protéine
        uniquePeptides (list) : la liste des peptides uniques pour une séquence

    Raises:
        /

    Returns:
        list : une liste contenant les séquences n'ayant pas de peptide unique
    """
    res = [i + 1 for i in range(len(peptideToProt))]
    for pep in uniquePeptides:
        if pep.get_nb_prot() in res:
            res.remove(pep.get_nb_prot())
    return res


# **********************************************************************************************************************


def searchCombinations(peptidesOfProteins, speciesOfPeptides, seqWithoutUnique, dicoNameSeq, seqToGenFam={},
                       threshold=0):
    """Cherche les combinaisons de peptides uniques pour chaque séquence

    Args:
        peptidesOfProteins (dict) : le dictionnaire associant les peptides à leur protéine
        speciesOfPeptides (dict) : le dictionnaire donnant les espèces dans lesquelles apparaissent chaque peptide
        seqWithoutUnique (list) : la liste des séquences n'ayant pas de peptide unique
        threshold (int) : seuil déterminant le nombre maximum de peptides dans une combinaison

    Raises:
        /

    Returns:
        dict : le dictionnaire de la combinaison pour chaque séquence
    """
    combinaisons = {}
    for seq in seqWithoutUnique:  # on a les séquences sans uniques
        allSequences = [i + 1 for i in range(len(peptidesOfProteins))]
        if not seqToGenFam:
            allSequences.remove(seq)
        else:
            for genFam, seqDico in seqToGenFam.items():
                if seq in seqDico:
                    for elt in seqDico:
                        allSequences.remove(elt)
        peptides = peptidesOfProteins[seq]  # on a les peptides des séquences sans uniques
        currentPeptide = peptides[0]
        usedPeptides = []
        while not allSequences == []:
            if len(usedPeptides) == len(peptides):
                if threshold != 0:
                    break
                else:
                    combinaisons[seq] = []
                    allSequences = []
            else:
                mini = 20000000
                for pep in peptides:  # pour chaque peptide
                    if pep not in usedPeptides:
                        species = speciesOfPeptides[pep]
                        if len(species) < mini:  # on prend celui qui a le moins de correspondances avec les autres séquences
                            mini = len(species)
                            currentPeptide = pep
                tmp = []
                for sequence in allSequences:  # on supprime les séquences qui ne contiennent pas le peptide
                    if sequence not in speciesOfPeptides[currentPeptide]:
                        tmp.append(sequence)
                for elt in tmp:
                    allSequences.remove(elt)
                if tmp:
                    if seq not in combinaisons:
                        combinaisons[seq] = [currentPeptide]
                    else:
                        combinaisons[seq].append(currentPeptide)
                for elt in peptides:
                    if elt.get_nb_prot() == currentPeptide.get_nb_prot() and elt.get_nb_peptide() == currentPeptide.get_nb_peptide():
                        usedPeptides.append(elt)
                if threshold != 0:
                    if len(combinaisons[seq]) >= threshold:
                        allSequences = []
    if not seqToGenFam:
        warning = []
        for s, p in combinaisons.items():
            if not p:
                warning.append(s)
        lineWarn = "\n*********\nLes séquences suivantes n\'ont ni peptide unique, ni combinaison de peptides uniques :\n"
        for elt in warning:
            theName = ""
            for name, nb in dicoNameSeq.items():
                if nb == elt:
                    theName = name
            lineWarn += theName
            lineWarn += '\n'
        lineWarn += '**********'
        print(lineWarn)
    return combinaisons


def createFile(output_dir, combinations, combsGenus, combsFamily, dicoNameSeq, peptideForProt):
    """Crée le fichier des résultats des combinaisons

    Args:
        output_file (str) : le nom du fichier de sortie
        allResultsFile (str) : le nom du fichier de sortie contenant tous les résultats
        combinations (dict) : les combinaisons à insérer
        combsWithThreshold (dict) : les combinaisons calculées avec un seuil à insérer
        dicoNameSeq (dict) : le dictionnaire contenant le nom d'une séquence avec son numéro associé
        peptideForProt (dict) : les peptides associés à leur protéine

    Raises:
        /

    Returns:
        list : une liste contenant les séquences n'ayant pas de peptide unique
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + 'combinations.csv', 'w', newline='') as results:
        with open(output_dir + 'allResults.csv', 'a', newline='') as allRes:
            writer_all = csv.writer(allRes, delimiter='|')
            writer_combi = csv.writer(results, delimiter='|')
            writer_combi.writerow(["List of unique combinations (for sequences that don't have unique peptides) :"])
            writer_combi.writerow(["Family", "Genus", "Name", "Peptide nb"])
            writer_all.writerow("")
            writer_all.writerow("")
            writer_all.writerow(["List of unique combinations (for sequences that don't have unique peptides) :"])
            writer_all.writerow(["Family", "Genus", "Name", "Peptide nb"])
            theName = ""
            for seq, peptides in combinations.items():
                if peptides:
                    for name, nb in dicoNameSeq.items():
                        if nb == seq:
                            theName = name
                    liste_pep_combi = []
                    for pep in peptides:
                        liste_pep_combi.append(str(pep.get_nb_peptide()))
                    pep = peptides[0]
                    writer_all.writerow([pep.get_family(), pep.get_genus(), theName, ",".join(liste_pep_combi)])
                    writer_combi.writerow([pep.get_family(), pep.get_genus(), theName, ",".join(liste_pep_combi)])

            writer_all.writerow("")
            writer_all.writerow("")
            writer_all.writerow(["List of unique combinations for each sequence, only taking into consideration sequences that have different genus"])
            writer_combi.writerow("")
            writer_combi.writerow("")
            writer_combi.writerow(["List of unique combinations for each sequence, only taking into consideration sequences that have different genus"])
            theName = ''
            for seqG, peptidesG in combsGenus.items():
                if peptides:
                    for name, nb in dicoNameSeq.items():
                        if nb == seqG:
                            theName = name
                liste_pep_combi_genus = []
                for pep in peptidesG:
                    liste_pep_combi_genus.append(pep)
                    peptoUse = pep
                writer_all.writerow([peptoUse.get_family(), peptoUse.get_genus(), theName, ",".join(map(lambda s: str(s.get_nb_peptide()), liste_pep_combi_genus))])
                writer_combi.writerow([peptoUse.get_family(), peptoUse.get_genus(), theName, ",".join(map(lambda s: str(s.get_nb_peptide()), liste_pep_combi_genus))])

            writer_all.writerow("")
            writer_all.writerow("")
            writer_all.writerow(["List of unique combinations for each sequence, only taking into consideration sequences that have different family"])
            writer_combi.writerow("")
            writer_combi.writerow("")
            writer_combi.writerow(["List of unique combinations for each sequence, only taking into consideration sequences that have different family"])
            theName = ''
            for seqF, peptidesF in combsFamily.items():
                if peptides:
                    for name, nb in dicoNameSeq.items():
                        if nb == seqF:
                            theName = name
                liste_pep_combi_family = []
                for pep in peptidesF:
                    liste_pep_combi_family.append(pep)
                    peptoUse = pep
                writer_all.writerow([peptoUse.get_family(), peptoUse.get_genus(), theName, ",".join(map(lambda s: str(s.get_nb_peptide()), liste_pep_combi_family))])
                writer_combi.writerow([peptoUse.get_family(), peptoUse.get_genus(), theName, ",".join(map(lambda s: str(s.get_nb_peptide()), liste_pep_combi_family))])


def prettyPrint_liste_peptides(output_dir, peptideToProtein, dicoNameSeq):
    """Crée l'affichage de la liste complète des peptides avec les positions sur les séquences'

        Args:
            output_file (str) : le nom du fichier de sortie
            dicoNameSeq (dict) : le dictionnaire contenant le nom d'une séquence avec son numéro associé
            peptideToProtein (dict) : les peptides associés à leur protéine

        Raises:
            /

        Returns:
            None
        """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + 'allResults.csv', 'a') as results:
        results.write("*********************************************\n")
        results.write("\nLISTE DES PEPTIDES POUR CHAQUE SÉQUENCE :")
        for prot, peps in peptideToProtein.items():
            for name, num in dicoNameSeq.items():
                if num == prot:
                    results.write("\n\nProtéine : {}".format(name))
            for peptide in peps:
                results.write("\n\tPeptide {} : position {} ; {}".format(peptide.get_nb_peptide(), peptide.get_position(), peptide.get_seq()))


def mainSearchCombinations(dicoPeptides, uniquesPeptides, output_dir, threshold=0):
    dicoSeq = dicoNameNumberSequence(dicoPeptides)
    speciesForPeptides = searchSpeciesForPeptides(dicoPeptides)
    peptidesToProtein = assignPeptideToProtein(speciesForPeptides)
    seqWithoutUnique = getSequencesWithoutUnique(peptidesToProtein, uniquesPeptides)

    seqToGen = seqToGenus(peptidesToProtein)
    seqToFam = seqToFamily(peptidesToProtein)

    seqWithoutUniqueGenre = mainGenre(dicoPeptides, output_dir, peptidesToProtein)
    seqWithoutUniqueFamily = mainFamily(dicoPeptides, output_dir, peptidesToProtein)

    combinations = searchCombinations(peptidesToProtein, speciesForPeptides, seqWithoutUnique, dicoSeq)
    combsFamily = searchCombinations(peptidesToProtein, speciesForPeptides, seqWithoutUniqueFamily, dicoSeq, seqToFam)
    combsGenus = searchCombinations(peptidesToProtein, speciesForPeptides, seqWithoutUniqueGenre, dicoSeq, seqToGen)

    createFile(output_dir, combinations, combsGenus, combsFamily, dicoSeq, peptidesToProtein)
    prettyPrint_liste_peptides(output_dir, peptidesToProtein, dicoSeq)
