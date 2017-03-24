import re

# Parser for .ins and .res SHELLXTL files
class SHELXTLFile():
    def __init__(self, filetxt):
        pre = True
        post = False
        self.pre_text = ""
        self.post_text = ""
        self.composition = []

        # object to store "commands"?
        self.commands = []
        self.named_params = {}

        self.fvar_vals = []

        # object to store crystal site information
        # ex. SM1   4    0.086929    0.000000    0.881287    10.50000    0.00621    0.00453 =
        #           0.00623    0.00000    0.00199    0.00000
        self.crystal_sites = []

        # crystal site info (for easier distance calculation and translation between atom sites and q sites)
        self.q_peaks = []

        self.r2 = 0.0



        for line in filetxt:
            if pre:
                self.pre_text += line
            if post:
                self.post_text += line



        re.search(r"SFAC ([^\n]+)", filetxt)

    # do actual file parsing
    def read(self, filetxt):
        pass

    # write ins file (with commands?)
    def write_ins(self):
        pass

    # various editing methods ...

