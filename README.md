# peptides_signature

## Objectif

Ce programme a pour but de trouver les peptides uniques à une espèce (dans le cas où une protéine = une espèce), à un genre ou à une famille, à partir d'une liste de peptides. Dans le cas où aucun peptide unique ne serait trouvé, une combinaison de peptides unique est recherchée. 

## Description

Ce dossier contient 2 dossiers : results (fichiers générés par le code, nos résultats) et fichier (fichiers utilisés, notamment le fichier csv contenant les peptides), suivi de notre code.

- peptide_v2 est notre classe représentant un peptide
- unique_peptides contient les fonctions utilisées pour chercher les peptides uniques dans les séquences.
- genre et famille contiennent les fonctions utilisées pour chercher les peptides uniques pour des familles ou genres
- combinaisons contient les fonctions pour chercher les combinaisons uniques de peptides
- Utils est un fichier où sont stockées les fonctions dont nous ne nous servons plus pour le moment
- Main contient uniquement la fonction main

## Contraintes techniques

- Python 3
- RPG 1.2.4

## Lancement du code

2 manières :

`
$ chmod u+x Main.py
`

`
$ ./Main.py [-h] [-p PEPTIDETHRESHOLD] [-c COMBINATIONTHRESHOLD] peptidesCSVFile resultsPath
`

Ou :

`
$ python3 Main.py [-h] [-p PEPTIDETHRESHOLD] [-c COMBINATIONTHRESHOLD] peptidesCSVFile resultsPath
`
### Paramètres

Obligatoires :
-  peptidesCSVFile : le nom du fichier contenant la liste des peptides (le fichier csv)
- resultsPath : le dossier (chemin) qui contiendra tous les fichiers résultats (ex : results/), le chemin doit se finir par un '/'.

Facultatifs :
- -h : message d'aide
- -p : la taille minimale des peptides gardés (par défaut 1)
- -c : la taille maximale d'une combinaison (par défaut aucune limite)


Exemple de commande fonctionnelle :

`
$ python3 Main.py fichier/COL1A1_trypsin.csv results/
`

### Syntaxe du fichier d'entrée

C'est un fichier CSV qui doit correspondre à une digestion in silico. Il peut être très facilement réalisable par RPG.

Syntaxe :

Nom_prot,No_peptide,Enzyme,Pos_clivage,Taille_pep,Masse_pep,pI,Sequence

## Fichiers résultats

- uniquePeptides.csv : contient les peptides uniques par protéine
- unique_pep_genre.csv : contient les peptides uniques par genre
- unique_pep_family.csv : contient les peptides uniques par famille
- combinations.csv : contient les combinaisons unique de peptides par protéine, genre et famille
- liste_peptides.txt : contient la liste de tous les peptides associés aux peptides qui leurs sont identiques.  

A noter que nous fournissons un fichier exemple "COL1A1_trypsin.csv" pour faire un test.

## Docker

Il est possible de lancer le programme en passant par Docker, en créant une image (à lancer depuis le répertoire où se trouve le code):

`
docker build --tag [nom de l'image] .
`

Et en lançant un conteneur Docker en liant un répertoire sur la machine hôte (pour permettre de retrouver les fichiers résultats créés) :

`
$ docker run -v chemin/répertoire/sur/la/machine/hôte:/app/results [nom de l'image]
`