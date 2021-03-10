#!/usr/bin/env python
# encoding: utf-8
"""
    :copyright: (c) 2021 by C. Fufezan, M.Koesters, J. Leufken, S. Schulze
    :licence: MIT, see LISCENSE for more details
"""

from __future__ import absolute_import
import sys
import re
from collections import defaultdict as ddict
import copy
from chemical_composition import chemical_composition_kb
import unimod_mapper


class ChemicalComposition(dict):
    """
    Chemical composition class.

    Keyword Arguments:
        sequence (Optional[str]): Peptide or chemical formula sequence
        aa_compositions (Optional[dict]): amino acid compositions
        isotopic_distributions (Optional[dict]): isotopic distributions
        monosaccharide_compositions (Optional[dict]): compositions of monosaccharides (as Hill notation)

    Keyword argument examples:

        **sequence** - Currently this can for example be::
              [
              'H2O',
              '{peptide}'.format(pepitde='ELVISLIVES'),
              '{peptide}#{unimod}:{pos}'.format(
              peptide = 'ELVISLIVES',
              unimod = 'Oxidation',
              pos = 1
              )
              ]

    Examples::
        >>> c = ursgal.ChemicalComposition()
        >>> c.use("ELVISLIVES#Acetyl:1")
        >>> c.hill_notation()
        'C52H90N10O18'
        >>> c.hill_notation_unimod()
        'C(52)H(90)N(10)O(18)'
        >>> c
        {'O': 18, 'H': 90, 'C': 52, 'N': 10}
        >>> c.composition_of_mod_at_pos[1]
        defaultdict(<class 'int'>, {'O': 1, 'H': 2, 'C': 2})
        >>> c.composition_of_aa_at_pos[1]
        {'O': 3, 'H': 7, 'C': 5, 'N': 1}
        >>> c.composition_at_pos[1]
        defaultdict(<class 'int'>, {'O': 4, 'H': 9, 'C': 7, 'N': 1})

        >>> c = ursgal.ChemicalComposition('H2O2H2')
        >>> c
        {'O': 2, 'H': 4}
        >>> c.subtract_chemical_formula('H3')
        >>> c
        {'O': 2, 'H': 1}

    """

    def __init__(
        self,
        sequence=None,
        aa_compositions=None,
        isotopic_distributions=None,
        monosaccharide_compositions=None,
    ):

        self._unimod_parser = None
        self.composition_of_mod_at_pos = {}
        """dict: chemical composition of unimod modifications at given position
        (if peptide sequence was used as input or using the `use` function)

        Note:
            Numbering starts at position 1, since all PSM search engines
            use this nomenclature.
        """
        self.composition_of_aa_at_pos = {}
        """dict: chemical composition of amino acid at given peptide position
        (if peptide sequence was used as input or using the `use` function)

        Note:
            Numbering starts at position 1, since all PSM search engines
            use this nomenclature.
        """
        #                 Examples::
        #                     c.composition_of_mod_at_pos[1] = {'15N': 2, '13C': 6, 'N': -2, 'C': -6}

        #         """
        self.composition_at_pos = {}
        """dict: chemical composition at given peptide position incl modifications
        (if peptide sequence was used as input or using the `use` function)

        Note:
            Numbering starts at position 1, since all PSM search engines
            use this nomenclature.
        """
        self.peptide = None
        self.addon = None
        self.unimod_at_pos = {}
        # self.regex_patterns = {
        #     ':pos' : re.compile( r''':(?P<pos>[0-9]*)''' ),
        #     'aaN'  : re.compile( r'''(?P<aa>[A-Z]{1})(?P<N>[0-9]*)''' ),
        # }
        if aa_compositions is None:
            self.aa_compositions = chemical_composition_kb.aa_compositions
        else:
            self.aa_compositions = aa_compositions
        if isotopic_distributions is None:
            self.isotopic_distributions = chemical_composition_kb.isotopic_distributions
        else:
            self.isotopic_distributions = isotopic_distributions
        if monosaccharide_compositions is None:
            self.monosaccharide_compositions = (
                chemical_composition_kb.monosaccharide_compositions
            )
        else:
            self.monosaccharide_compositions = monosaccharide_compositions
        if sequence is not None:
            self.use(sequence)

        self.isotope_mass_lookup = {}
        for element, isotope_list in self.isotopic_distributions.items():
            for isotope_mass, abundance in isotope_list:
                isotope_mass_key = "{0}{1}".format(
                    str(round(isotope_mass)).split(".")[0], element
                )
                self.isotope_mass_lookup[isotope_mass_key] = isotope_mass

    def __add__(self, other_cc):
        """
        Add chemical composition (dict) to the class

        Arguments:
            other_cc (dict): dictionary containing the number of atoms of each element that will be added
        """
        tmp = copy.deepcopy(self)
        for other_key, other_value in other_cc.items():
            tmp[other_key] += other_value
        return tmp

    def __missing__(self, key):
        """Set missing elements to 0"""
        if key not in self.keys():
            self[key] = 0
        return self[key]

    def clear(self):
        """Resets all lookup dictionaries and self

        One class instance can be used analysing a series of sequences, thereby
        avoiding class instantiation overhead
        """

        self.composition_of_mod_at_pos.clear()
        self.composition_of_aa_at_pos.clear()
        self.composition_at_pos.clear()
        self.unimod_at_pos.clear()

        self.peptide = None
        self.addon = None
        for k in list(self.keys()):
            del self[k]

    def _parse_sequence_unimod_style(self, sequence):
        """
        Parse a sequence of the format "<PEPTIDE>#<unimod>:<pos>"
        with unimod referrig to any unimod PSI-MS name and
        pos referring to the position of the modification within the peptide sequence.
        Modifications of peptide side chains start at position 1,
        while position 0 refers to the N-terminus and position len(peptide)+1 refers
        to the C-terminus.

        Args:
            sequence(str): sequence in unimod style
        """
        if self._unimod_parser is None:
            self._unimod_parser = unimod_mapper.UnimodMapper()
        minPos = sequence.index("#")
        peptide = sequence[:minPos]
        addon = sequence[minPos + 1 :]
        self.peptide = peptide
        if peptide != "":
            self.add_peptide(peptide)
            self["O"] += 1
            self["H"] += 2
        self.addon = addon
        unimods = self.addon.split(";")
        pattern = re.compile(r""":(?P<pos>[0-9]*$)""")
        for unimod in unimods:
            if unimod == "":
                continue
            unimod = unimod.strip()
            if ":" not in unimod:
                sys.exit(
                    """
                    Error in chemical_composition.py:
                    This unimod: {0} requires positional information
                    """.format(
                        unimod
                    )
                )

            for occ, match in enumerate(pattern.finditer(unimod)):
                unimodcomposition = self._unimod_parser.name2composition(
                    unimod[: match.start()]
                )
                if unimodcomposition is None:
                    sys.exit(
                        """
                        Error in chemical_composition.py:
                        Cannot map unimod {0}
                        """.format(
                            unimod[: match.start()]
                        )
                    )
                if occ >= 1:
                    sys.exit(
                        """
                        Error in chemical_composition.py:
                        The unimod {0} contains multiple ":", preventing to map the position correctly
                        """.format(
                            unimod
                        )
                    )
                position = int(match.group("pos"))
                if position not in self.unimod_at_pos.keys():
                    self.unimod_at_pos[position] = []
                self.unimod_at_pos[position].append(unimod[: match.start()])

            for k, v in unimodcomposition.items():
                self[k] += v
            # Storing position related modifications
            position = int(match.group("pos"))
            if position == 0:
                # E.g. Acetylation at pos 0 indicates N-Term
                # but has to be counted for position 1 in this class
                position = 1

            if position not in self.composition_of_mod_at_pos.keys():
                self.composition_of_mod_at_pos[position] = ddict(int)
            if position not in self.composition_at_pos.keys():
                self.composition_at_pos[position] = ddict(int)
            for k, v in unimodcomposition.items():
                self.composition_of_mod_at_pos[position][k] += v
                self.composition_at_pos[position][k] += v
        return

    def use(self, sequence):
        """Re-initialize the class with a new sequence

        This is helpful if one wants to use the same class instance
        for multiple sequence since it remove class instantiation overhead.

        Args:
            sequence (str): See top for possible input formats.
        """

        self.clear()
        # reset the shiznit
        if "#" in sequence:
            # Unimod Style format
            self._parse_sequence_unimod_style(sequence)
        elif bool(re.search(r"\d", sequence)):
            self.add_chemical_formula(sequence)
        else:
            self.add_amino_acids(sequence)
            self["O"] += 1
            self["H"] += 2
        return

    def add_chemical_formula(self, chemical_formula, factor=1):
        """Adds chemical formula to the instance

        Args:
            chemical_formula (str): chemical composition given as Hill notation

        Keyword Arguments:
            factor (int): multiplication factor to add the same chemical formula
                multiple times
        """
        self._merge(chemical_formula, mode="addition", factor=factor)
        return

    def add_amino_acids(self, aa_sequence):
        """
        Adds a sequence of amino acids to the instance

        Args:
            aa_sequence (str): string of amino acids (only standard amino acids and U allowed so far), e.g. "PEPTIDE"

        Note:
            In contrast to the previous function add_peptide, this function
            does not allow for labeled amino acids to be added, please use unimods
            for that purpose instead.
        """
        if self.peptide is None:
            self.peptide = ""
        for pos, aa in enumerate(aa_sequence):
            aa_pos = len(self.peptide) + pos + 1
            try:
                aa_compo = self.aa_compositions[aa]
            except:
                sys.exit(
                    """
                    Error in chemical_composition.py:
                    Do not know aa composition for {0}
                    in {1}
                    """.format(
                        aa, aa_sequence
                    )
                )
            self.add_chemical_formula(aa_compo)

            composition = self._chemical_formula_to_dict(aa_compo)
            self.composition_of_aa_at_pos[pos] = composition
            if pos not in self.composition_at_pos.keys():
                self.composition_at_pos[pos] = ddict(int)
            for k, v in composition.items():
                self.composition_at_pos[pos][k] += v
        self.peptide += aa_sequence
        return

    def add_peptide(self, peptide):
        """Adds peptide sequence to the instance"""
        # pattern = self.regex_patterns['aaN']
        print(
            """
            [Warning] You are using the function "add_peptide", which has been replaced with "add_amino_acids"
            """
        )
        pattern = re.compile(r"""(?P<aa>[A-Z]{1})(?P<N>[0-9]*)""")
        # pattern = re.compile(r'[A-Z]{1}[0-9]*')
        number_offset = 0
        # this are the count for e.g. SILAC aa, i.e. R0 R1 or C0 and so on ...
        # print( peptide, type( peptide ))
        for aaN_match in pattern.finditer(peptide):
            aa = aaN_match.group("aa")
            N = aaN_match.group("N")
            pos = int(aaN_match.start()) - number_offset + 1
            if N != "":
                number_offset += len(N)
            try:
                aa_compo = self.aa_compositions[aa + N]
            except:
                sys.exit(
                    """
                    Error in chemical_composition.py:
                    Do not know aa composition for {0}
                    in {1}
                    """.format(
                        aa + N, peptide
                    )
                )
            self.add_chemical_formula(aa_compo)

            composition = self._chemical_formula_to_dict(aa_compo)
            self.composition_of_aa_at_pos[pos] = composition
            if pos not in self.composition_at_pos.keys():
                self.composition_at_pos[pos] = ddict(int)
            for k, v in composition.items():
                self.composition_at_pos[pos][k] += v

    def add_glycan(self, glycan):
        """Adds a glycan to the instance.

        Args:
            glycan (str): sequence of monosaccharides given in unimod format,
                e.g.: HexNAc(2)Hex(3)dHex(1)Pent(1),
                available monosaccharides are listed in chemical_composition_kb
        """
        pattern = re.compile(r"""(?P<monosacch>[A-z0-9]*)(?P<count>\([0-9]*\))""")
        for glyc_match in pattern.finditer(glycan):
            monosacch = glyc_match.group("monosacch")
            if glyc_match.group("count") == "":
                count = 1
            else:
                count = int(glyc_match.group("count").strip("(").strip(")"))
            if monosacch in self.monosaccharide_compositions.keys():
                monosacch_compo = self.monosaccharide_compositions[monosacch]
            else:
                sys.exit("Do not know glycan composition for {0}".format(monosacch))
            self.add_chemical_formula(monosacch_compo, factor=count)
        return

    def _chemical_formula_to_dict(self, chemical_formula):
        """
        Converts chemical formula into chemical composition dictionary
        """
        unimod_style = False
        if "(" in chemical_formula:
            unimod_style = True
        chem_dict = {}
        if unimod_style:
            pattern = re.compile(
                r"(?P<isotope>[0-9]*)(?P<element>[A-Z][a-z]*)\((?P<count>[0-9]*)\)"
            )
        else:
            pattern = re.compile(r"(?P<element>[A-Z][a-z]*)(?P<count>[0-9]*)")
        for match in pattern.finditer(chemical_formula):
            if match.group("count") == "":
                count = 1
            else:
                count = int(match.group("count"))
            if unimod_style:
                element_key = "{0}{1}".format(
                    match.group("isotope"), match.group("element")
                )
            else:
                element_key = match.group("element")
            if element_key not in chem_dict.keys():
                chem_dict[element_key] = 0
            chem_dict[element_key] += count
        return chem_dict

    def hill_notation(self, include_ones=False, cc=None):
        """
        Formats chemical composition into `Hill notation`_ string.

        .. _Hill Notation:
            https://en.wikipedia.org/wiki/Hill_system

        Keyword Arguments:
            cc (dict, optional): can format other element dicts as well.
            include_ones (bool): Include "1" for elements with a count of 1 (otherwise the count will be omitted)

        Returns:
            str: Hill notation format of self.
                For examples::
                    C50H88N10O17
        """
        MAJORS = ["C", "H"]
        s = ""
        if cc is None:
            cc_dict = self
        else:
            cc_dict = cc

        for major in MAJORS:
            if major in cc_dict.keys():
                if cc_dict[major] == 0:
                    continue
                s += major
                if include_ones or cc_dict[major] > 1:
                    s += str(cc_dict[major])
        for k in sorted(cc_dict.keys()):
            if k not in MAJORS:
                if cc_dict[k] == 0:
                    continue
                s += k
                if include_ones or cc_dict[k] > 1:
                    s += str(cc_dict[k])
        return s

    def hill_notation_unimod(self, cc=None):
        """
        Formats chemical composition into `Hill notation`_ string
        adding `unimod`_ features.

        .. _Hill Notation:
            https://en.wikipedia.org/wiki/Hill_system

        .. _unimod:
            http://www.unimod.org/fields.html

        Args:
            cc (dict, optional): can format other element dicts as well.

        Returns:
            str: Hill notation format including unimod format rules of self.
                For example::
                    C(50)H(88)N(10)O(17)

        """
        MAJORS = ["C", "H"]
        s = ""
        if cc is None:
            cc_dict = self
        else:
            cc_dict = cc

        for major in MAJORS:
            if major in cc_dict.keys():
                if cc_dict[major] == 0:
                    continue
                s += "{0}({1})".format(
                    major.replace("(", "").replace(")", ""), cc_dict[major]
                )
        for k in sorted(cc_dict.keys()):
            if k not in MAJORS:
                if cc_dict[k] == 0:
                    continue
                s += "{0}({1})".format(k.replace("(", "").replace(")", ""), cc_dict[k])
        return s

    def _mass(self, cc=None):
        """
        Calculate the mass of the chemical composition.
        Optional cc can be specified, i.e. a cc dict in the style of
        { element : count , ... }

        This does however not work with enriched elements, e.g. 15N
        Better use pyQms isotopologue library for more accurate mass
        calculations.
        """
        mass = 0
        if cc is None:
            cc_mass_dict = self
        else:
            cc_mass_dict = cc
        for element, count in cc_mass_dict.items():
            isotope_mass = None
            try:
                isotope_mass = sorted(
                    self.isotopic_distributions[element],
                    key=lambda x: x[1],
                    reverse=True,
                )[0][0]
                ##### add test with Li
            except:
                # we check for 15N or 13C or other isotopes
                assert (
                    element in self.isotope_mass_lookup.keys()
                ), """
                The given element {0}
                is not included yet. Please add it to the isotopic_distributions
                dictionary in ursgal/chemical_composition_kb.py.
                Isotopic masses and distributions can be found at
                http://ciaaw.org/atomic-masses.html
                """.format(
                    element
                )
                isotope_mass = self.isotope_mass_lookup[element]
            mass += count * isotope_mass

        return mass

    def _merge(self, chemical_formula, mode="addition", factor=1):
        """
        Generalized function that allows addition and subtraction
        of chemical formulas to/from the instance

        Args:
            chemical_formula (str|dict): chemical formula that should be merged with the instance

        Keyword Arguments:
            mode (str): "addition" to add the defined chemical formula, "subtraction" to subtract it from the instance
            factor (int): defines how often the chemical formula should be added/subtracted
        """
        if mode == "addition":
            sign = +1
        elif mode == "subtraction":
            sign = -1
        else:
            sys.exit(
                """
                Error in chemical_composition.py:
                Do not know which mode to use for _merge
            """
            )
        if isinstance(chemical_formula, str):
            chemical_formula = self._chemical_formula_to_dict(chemical_formula)
        for element, count in chemical_formula.items():
            self[element] = self[element] + sign * count * factor
        return

    def subtract_amino_acids(self, aa_sequence):
        """Subtract amino acid sequence from instance"""
        for aa in aa_sequence:
            try:
                aa_compo = self.aa_compositions[aa]
            except:
                sys.exit(
                    """
                    Error in chemical_composition.py:
                    Do not know aa composition for {0}
                    in {1}
                    """.format(
                        aa, aa_sequence
                    )
                )
            self.subtract_chemical_formula(aa_compo)

        if self.peptide.endswith(aa_sequence):
            self.peptide = "".join(self.peptide.split(aa_sequence)[:-1])
            for pos, aa in enumerate(aa_sequence):
                aa_pos = len(self.peptide) + pos + 1
                del self.composition_of_aa_at_pos[pos]
                del self.composition_at_pos[pos]
        else:
            print(
                """
                [WARNING] The subtracted amino acid sequence does not match the end of the stored peptide sequence.
                The compositions at specific positions can therefore not be updated!
                """
            )
        return

    def subtract_chemical_formula(self, chemical_formula, factor=1):
        """Subtract chemical formula from instance

        Args:
            chemical_formula (str): chemical composition given as Hill notation

        Keyword Arguments:
            factor (int): multiplication factor to add the same chemical formula
                multiple times
        """
        self._merge(chemical_formula, mode="subtraction", factor=factor)
        return

    def generate_cc_dict(self):
        """
        Generate a dictionary of the instance
        """
        tmp = {}
        tmp.update(self)
        return tmp
