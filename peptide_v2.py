class Peptides(object):
    """Classe pour représenter les peptides donnés par RPG
    """

    def __init__(self, seq, nb_peptide, prot_name, nb_prot, genus, family, position, mass):
        self.__seq = seq
        self.__nb_peptide = nb_peptide
        self.__prot_name = prot_name
        self.__nb_prot = nb_prot
        self.__genus = genus
        self.__family = family
        self.__position = position
        self.__mass = mass

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
    
    def get_mass(self):
        return self.__mass

    def __str__(self):
        return f"Peptide's  number : {self.__nb_peptide}; position : {self.__position}; protein name : {self.__prot_name}{self.__nb_prot}; seq : {self.__seq}; family : {self.__family}; genus : {self.__genus}; mass : {self.__mass} "

    def __repr__(self):
        return f'Peptide(seq={self.__seq}, nb={self.__nb_peptide}, protein name={self.__prot_name}, nb protein={self.__nb_prot}, genus={self.__genus}, family={self.__family}, position={self.__position}, mass={self.__mass})'