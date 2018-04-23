from pymatgen import Element

class SHELXElement:
    """
    Define class to hold information about an element in the context of a SHELX File
    """
    def __init__(self, name, stoichiometry, index):
        self.name = name
        self.nominal_stoichiometry = stoichiometry
        self.index = index

    def get_name(self, capitalize=False):
        """
        Get the element name
        :param capitalize: Whether to capitalize conventionally (Helium would be He)
        :return: element name
        """
        if capitalize:
            return self.name.capitalize()
        else:
            return self.name

    def get_pymatgen_element(self):
        """
        Get a pymatget element object that corresponds to this element
        :return:
        """
        return Element(self.get_name(True))