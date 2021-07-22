from family import *
from genus import *
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
        if warning:
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
        writer_combi = csv.writer(results)
        writer_combi.writerow(["List of unique combinations (for sequences that don't have unique peptides) :"])
        writer_combi.writerow(["Family", "Genus", "Protein Name", "Position", "Peptide mass", "Sequence"])
        theName = ""
        nocombiName = []
        for seq, peptides in combinations.items():
            if peptides:
                for name, nb in dicoNameSeq.items():
                    if nb == seq:
                        theName = name
                toAppend = []
                toAppend.append(peptides[0].get_family()) ; toAppend.append(peptides[0].get_genus()); toAppend.append(theName);  toAppend.append(peptides[0].get_position()); toAppend.append(peptides[0].get_mass()); toAppend.append(peptides[0].get_seq())
                for pep in peptides[1:]:
                    toAppend.append(pep.get_position()); toAppend.append(pep.get_mass()); toAppend.append(pep.get_seq())
                writer_combi.writerow(toAppend)
            else:
                nocombi = []
                nocombi.append(seq)
                for elt in nocombi:
                    theName = ""
                    for name, nb in dicoNameSeq.items():
                        if nb == elt:
                            theName = name
                            nocombiName.append(theName)
        strNoCombiName = ",".join(nocombiName)                    
        writer_combi.writerow("")
        writer_combi.writerow(["Sequences that don't have a unique peptide nore a unique combination :" + strNoCombiName])


        writer_combi.writerow("")
        writer_combi.writerow("")
        writer_combi.writerow(["List of unique combinations for each sequence, only taking into consideration sequences that have different genus"])
        theName = ''
        nocombiNameG = []
        for seqG, peptidesG in combsGenus.items():
            if peptidesG:
                for name, nb in dicoNameSeq.items():
                    if nb == seqG:
                        theName = name
                toAppendG = []
                toAppendG.append(peptidesG[0].get_family()) ; toAppendG.append(peptidesG[0].get_genus()); toAppendG.append(theName);  toAppendG.append(peptidesG[0].get_position()); toAppendG.append(peptidesG[0].get_mass()); toAppendG.append(peptidesG[0].get_seq())
                for pep in peptidesG[1:]:
                    toAppendG.append(pep.get_position()); toAppendG.append(pep.get_mass()); toAppendG.append(pep.get_seq())
                writer_combi.writerow(toAppendG)
            else:
                nocombiG = []
                nocombiG.append(seqG)
                for elt in nocombiG:
                    theName = ""
                    for name, nb in dicoNameSeq.items():
                        if nb == elt:
                            theName = name
                            nocombiNameG.append(theName)
        strNoCombiNameG = ",".join(nocombiNameG)                    
        writer_combi.writerow("")
        writer_combi.writerow(["Sequences that don't have a unique peptide nore a unique combination :" + strNoCombiNameG])

        writer_combi.writerow("")
        writer_combi.writerow("")
        writer_combi.writerow(["List of unique combinations for each sequence, only taking into consideration sequences that have different family"])
        theName = ''
        nocombiNameF = []
        for seqF, peptidesF in combsFamily.items():
            if peptidesF:
                for name, nb in dicoNameSeq.items():
                    if nb == seqF:
                        theName = name
                toAppendF = []
                toAppendF.append(peptidesF[0].get_family()) ; toAppendF.append(peptidesF[0].get_genus()); toAppendF.append(theName);  toAppendF.append(peptidesF[0].get_position()); toAppendF.append(peptidesF[0].get_mass()); toAppendF.append(peptidesF[0].get_seq())
                for pep in peptidesF[1:]:
                    toAppendF.append(pep.get_position()); toAppendF.append(pep.get_mass()); toAppendF.append(pep.get_seq())
                writer_combi.writerow(toAppendF)
            else:
                nocombiF = []
                nocombiF.append(seqF)
                for elt in nocombiF:
                    theName = ""
                    for name, nb in dicoNameSeq.items():
                        if nb == elt:
                            theName = name
                            nocombiNameF.append(theName)
        strNoCombiNameF = ",".join(nocombiNameF)                    
        writer_combi.writerow("")
        writer_combi.writerow(["Sequences that don't have a unique peptide nore a unique combination :" + strNoCombiNameF])

'''
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
'''

def mainSearchCombinations(dicoPeptides, uniquesPeptides, output_dir, threshold=0):
    dicoSeq = dicoNameNumberSequence(dicoPeptides)
    speciesForPeptides = searchSpeciesForPeptides(dicoPeptides)
    peptidesToProtein = assignPeptideToProtein(speciesForPeptides)
    seqWithoutUnique = getSequencesWithoutUnique(peptidesToProtein, uniquesPeptides)

    seqToGen = seqToGenus(peptidesToProtein)
    seqToFam = seqToFamily(peptidesToProtein)

    seqWithoutUniqueGenre = mainGenre(dicoPeptides, output_dir, peptidesToProtein)
    seqWithoutUniqueFamily = mainFamily(dicoPeptides, output_dir, peptidesToProtein)

    combinations = searchCombinations(peptidesToProtein, speciesForPeptides, seqWithoutUnique, dicoSeq, {}, threshold)
    combsFamily = searchCombinations(peptidesToProtein, speciesForPeptides, seqWithoutUniqueFamily, dicoSeq, seqToFam, threshold)
    combsGenus = searchCombinations(peptidesToProtein, speciesForPeptides, seqWithoutUniqueGenre, dicoSeq, seqToGen, threshold)

    createFile(output_dir, combinations, combsGenus, combsFamily, dicoSeq, peptidesToProtein)
    #prettyPrint_liste_peptides(output_dir, peptidesToProtein, dicoSeq)
