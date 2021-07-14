def searchCombinationsBis(peptidesOfProteins, speciesOfPeptides, seqWithoutUnique, threshold=0):
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
    return combinaisons


def combination_threshold(combinations, peptidesToProtein, speciesForPeptides, method, threshold):
    """Cherche les combinaisons de peptides uniques pour chaque séquence

    Args:
        combinations (dict) : le dictionnaire des combinaisons exactes trouvées (combinaisons uniques pour une seule séquence)
        peptideToProteins (dict) : le dictionnaire associant les peptides à leur protéine
        speciesForPeptides (dict) : le dictionnaire donnant les espèces dans lesquelles apparaissent chaque peptide
        threshold (int) : seuil déterminant le nombre maximum de peptides dans une combinaison

    Raises:
        /

    Returns:
        dict : le dictionnaire de la combinaison pour chaque séquence
    """
    pass
    # protWithoutCombi = []
    # for numSeq, peptides in combinations.items():
    #     if not peptides:
    #         protWithoutCombi.append(numSeq)
    # return searchCombinations(peptidesToProtein, speciesForPeptides, protWithoutCombi, method, threshold)


def createFileBis(output_file, allResultsFile, combinations, combsWithThreshold, dicoNameSeq, peptideForProt):
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
    pass
    # with open(output_file, 'w') as results:
    #     with open(allResultsFile, 'a') as allRes:
    #         allRes.write("\n\nListe des combinaisons uniques (pour les séquences n'ayant pas de peptide unique) : \n\n")
    #         theName = ""
    #         for seq, peptides in combinations.items():
    #             if peptides:
    #                 for name, nb in dicoNameSeq.items():
    #                     if nb == seq:
    #                         theName = name
    #                 lineToInsert = "Protéine {} ({}) : {}\n".format(seq, theName, '|'.join(
    #                     map(lambda s: str(s.get_nb_peptide()).strip('[]').replace(' ', ''), peptides)))
    #                 results.write(lineToInsert)
    #                 allRes.write(lineToInsert)
    #
    #             else:
    #                 for seqT, peptidesT in combsWithThreshold.items():
    #                     if seq == seqT:
    #                         for name, nb in dicoNameSeq.items():
    #                             if nb == seqT:
    #                                 theName = name
    #                         lineToWrite = "Protéine {} ({}) : {} (combinaison unique pour les séquences {})\n".format(
    #                             seqT, theName, '|'.join(
    #                                 map(lambda s: str(s.get_nb_peptide()).strip('[]').replace(' ', ''), peptidesT)),
    #                             ','.join(map(lambda s: str(s).strip('[]').replace(' ', ''),
    #                                          getSequencesForPeptides(
    #                                              peptidesT,
    #                                              peptideForProt))))
    #                         results.write(lineToWrite)
    #                         allRes.write(lineToWrite)
