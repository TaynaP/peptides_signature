class Peptides(object):
    """Classe pour représenter les peptides donnés par RPG
    """

    def __init__(self, seq, nb_peptide, prot_name, nb_prot, genus, family, position):
        self.__seq = seq
        self.__nb_peptide = nb_peptide
        self.__prot_name = prot_name
        self.__nb_prot = nb_prot
        self.__genus = genus
        self.__family = family
        self.__position = position

    # Getter
    def get_seq(self):
        return self.__seq

    def get_nb_peptide(self):
        return self.__nb_peptide

    def get_prot_name(self):
        return self.__prot_name

    def get_nb_prot(self):
        return self.__nb_prot

    def get_genus(self):
        return self.__genus

    def get_family(self):
        return self.__family

    def get_position(self):
        return self.__position
