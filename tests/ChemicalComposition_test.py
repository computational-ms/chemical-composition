#!/usr/bin/env python
# encoding: utf-8
import chemical_composition
import unittest


class TestChemicalComposition(unittest.TestCase):
    def setUp(self):
        self.lib = chemical_composition.ChemicalComposition()
        self.water1 = chemical_composition.ChemicalComposition(deprecated_format="+H2O")
        self.water_unimod_style_chemical_formula = (
            chemical_composition.ChemicalComposition(deprecated_format="+H(2)16O(1)")
        )

    def tearDown(self):
        pass

    def generate_cc_using_uni_mod_terminology_test(self):
        cc_returned = self.water_unimod_style_chemical_formula.generate_cc_dict()
        assert len(cc_returned.keys()) == 2
        self.assertEqual(cc_returned["16O"], 1)
        self.assertEqual(cc_returned["H"], 2)
        return

    def test_add_chemical_formulas(self):
        water1 = chemical_composition.ChemicalComposition(deprecated_format="+H2O")
        water2 = chemical_composition.ChemicalComposition(deprecated_format="+H2O")
        two_waters = water1 + water2
        self.assertEqual(two_waters.hill_notation(), "H4O2")

    # def test_representation_is_hill_notatation(self):
    #     self.assertEqual(str(self.water1), "H2O")

    def test_clear_removes_everything(self):
        cc = chemical_composition.ChemicalComposition(sequence="KLEINERTEST")
        cc.clear()
        self.assertEqual(cc.hill_notation(), "")

    def test_mass_even_though_you_should_not_use_it(self):
        cc = chemical_composition.ChemicalComposition(formula="H2O")
        self.assertAlmostEqual(cc.mass(), 18.0105646844)

    def test_mass_with_heavy_isotopes_TMT6Plex_example(self):
        cc = chemical_composition.ChemicalComposition(modifications="TMT6plex:0")
        self.assertAlmostEqual(cc.mass(cc=cc), 229.162932, 6)

    def fail_test(self):
        cc = chemical_composition.ChemicalComposition()
        failing_sequence_list = ["#Oxidation_w/o_pos", "#Oxidation_with_pos:1"]
        for failing_sequence in failing_sequence_list:
            with self.assertRaises(SystemExit) as system_exit_check:
                cc.add_modifications(failing_sequence)

            assert len(system_exit_check.exception.code) != 0
        cc.clear()
        cc.add_modifications(";")
        assert len(cc.keys()) == 0

        return

    def generate_cc_test(self):
        cc = chemical_composition.ChemicalComposition(formula="+H2O")
        cc_returned = cc.generate_cc_dict()
        assert len(cc_returned.keys()) == 2
        return


TESTS = [
    {
        "input": "A#Oxidation:1",
        "keyword": "deprecated_format",
        "output": "C(3)H(7)N(1)O(3)",
        "mod_pos_info": [{"pos": 1, "cc_mods": {"O": 1}}],
    },
    {
        "input": "A#Oxidation:1;Acetyl:0",
        "keyword": "deprecated_format",
        # zero should end up at pos 1 ...
        "output": "C(5)H(9)N(1)O(4)",
        "mod_pos_info": [{"pos": 1, "cc_mods": {"O": 2, "C": 2, "H": 2}}],
    },
    {
        "input": "X",
        "keyword": "deprecated_format",
        "aa_compositions": {"X": "C12"},
        "output": "C(12)H(2)O(1)",
        "mod_pos_info": [{"pos": 1, "cc_mods": None}],
    },
    # {
    #     "input": "RAA",
    #     # 'aa_compositions' :  {'X': 'C12', 'R0': 'C12'},
    #     "output": "C(12)H(24)N(6)O(4)",
    # },
    # {
    #     "input": "R0X",
    #     "aa_compositions": {"X": "C12", "R0": "C12"},
    #     "output": "C(24)H(2)O(1)",
    #     # 'mod_pos_info' : (0, [{'C': 12 }], 1, {'O' : 1, 'C' : 2, 'H' : 2}])
    # },
    {
        "input": "R#Label:13C(6)15N(2):1",
        "keyword": "deprecated_format",
        "output": "H(14)13C(6)15N(2)N(2)O(2)",
    },
    # {
    #     "input": "R0XR1",
    #     "aa_compositions": {"R0": "C12", "X": "C13", "R1": "C14"},
    #     "output": "C(39)H(2)O(1)",
    #     "aa_pos_info": [
    #         {"pos": 1, "cc_mods": {"C": 12}},
    #         {"pos": 3, "cc_mods": {"C": 14}},
    #     ],
    # },
    {
        "input": "RR#Label:13C(6)15N(2):1",
        "keyword": "deprecated_format",
        "output": "C(6)H(26)13C(6)15N(2)N(6)O(3)",
        "aa_pos_info": [
            {"pos": 1, "cc_mods": {"C": 6, "O": 1, "N": 4, "H": 12}},
            {"pos": 2, "cc_mods": {"C": 6, "O": 1, "N": 4, "H": 12}},
        ],
        "mod_pos_info": [{"pos": 1, "cc_mods": {"C": -6, "13C": 6, "N": -2, "15N": 2}}],
    },
    {
        "input": "+H2O2H2",
        "keyword": "deprecated_format",
        "output": "H(4)O(2)",
    },
    # {"input": "H2O2H2-HO", "output": "H(3)O(1)"},
]


# SUBTRACTION_TEST_SET = [
#     # functionality tested yet ...
#     {"in": "+H2O2H2", "sub": "H2", "out": {"H": 2, "O": 2}}
# ]


def pepitde_with_unimod_test():
    for test_id, test_dict in enumerate(TESTS):
        yield check_unimod_pep_input, test_dict


def check_unimod_pep_input(test_dict):
    if "aa_compositions" not in test_dict.keys():
        test_dict[
            "aa_compositions"
        ] = chemical_composition.chemical_composition_kb.aa_compositions
    if test_dict["keyword"] == "deprecated_format":
        cc = chemical_composition.ChemicalComposition(
            deprecated_format=test_dict["input"],
            aa_compositions=test_dict["aa_compositions"],
        )
    else:
        pass
    print("hill:", cc.hill_notation_unimod())
    print(cc.composition_of_aa_at_pos)
    print(cc.composition_of_mod_at_pos)
    assert cc.hill_notation_unimod() == test_dict["output"]
    if "mod_pos_info" in test_dict.keys():
        for mod_dict in test_dict["mod_pos_info"]:
            cc_mods = cc.composition_of_mod_at_pos.get(mod_dict["pos"], None)
            assert cc_mods == mod_dict["cc_mods"]
    if "aa_pos_info" in test_dict.keys():
        for aa_dict in test_dict["aa_pos_info"]:
            cc_aa = cc.composition_of_aa_at_pos.get(aa_dict["pos"], None)
            assert cc_aa == aa_dict["cc_mods"]


def Simple_amino_acid_test():
    for (
        aa,
        chemformula,
    ) in chemical_composition.chemical_composition_kb.aa_compositions.items():
        yield check_hill_notation, aa, chemformula


def check_hill_notation(aa, chemformula):
    cc = chemical_composition.ChemicalComposition(aa)
    cc.subtract_chemical_formula("H2O")
    # print(aa, chemformula)
    # print(cc)
    assert cc.hill_notation() == chemformula


if __name__ == "__main__":
    Simple_amino_acid_test()
