# signatureCollagene

### Idée 1 : Sans alignement multiple

Le code à utiliser se situe dans le dossier idee_1.

Ce dossier contient 2 dossiers : results (fichiers générés par le code, nos résultats) et fichier (fichiers utilisés, notamment le fichier fasta contenant les peptides), suivi de notre code.

- peptide_v2 est notre classe représentant un peptide
- unique_peptides contient les fonctions utilisées pour chercher les peptides uniques dans les séquences.
- genre et famille contiennent les fonctions utilisées pour chercher les peptides uniques pour des familles ou genres
- combinaisons contient les fonctions pour chercher les combinaisons uniques de peptides
- Utils est un fichier où sont stockées les fonctions dont nous ne nous servons plus pour le moment
- Main contient uniquement la fonction main

##### Lancement du code :

`
$ python3 Main.py [fasta sorti RPG] [nom du fichier avec tout les peptides] [chemin du fichier où la liste des peptides et les équivalents seront écrits] [chemin du fichier où les peptides uniques seront écrits] [chemin du fichier où les combinaisons seront écrites] [chemin du fichier regroupant tous les résultats]
`

Exemple de commande fonctionnelle :

`
python3 Main.py fichier/rpg_trypsine_COL1A1_famille.fasta results/liste_peptides.txt results/uniquePeptides.txt results/combinations.txt results/allResults.txt
`

- À noter qu'une liste indiquant les numéros de peptides et leur position sur la séquence est disponible dans le fichier allResults.txt.

