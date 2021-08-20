#!/usr/bin/env python3.2
# encoding: utf-8
import chemical_composition
from pathlib import Path
from chemical_composition import chemical_composition_kb

PEPTIDE_TESTS = [
    {
        "input": "A#Oxidation:1",
        "output": "C(3)H(7)N(1)O(3)",
        "mod_pos_info": [{"pos": 1, "cc_mods": {"O": 1}}],
    },
    {
        "input": "A#Oxidation:1;Acetyl:0",
        # zero should end up at pos 1 ...
        "output": "C(5)H(9)N(1)O(4)",
        "mod_pos_info": [{"pos": 1, "cc_mods": {"O": 2, "C": 2, "H": 2}}],
    },
    {
        "input": "X",
        "aa_compositions": {"X": "C12"},
        "output": "C(12)H(2)O(1)",
        "mod_pos_info": [{"pos": 1, "cc_mods": None}],
    },
    {
        "input": "X#",
        "aa_compositions": {"X": "C12"},
        "output": "C(12)H(2)O(1)",
        "mod_pos_info": [{"pos": 1, "cc_mods": None}],
    },
    {
        "input": "RAA",
        # 'aa_compositions' :  {'X': 'C12', 'R0': 'C12'},
        "output": "C(12)H(24)N(6)O(4)",
    },
    {
        "input": "R0X",
        "aa_compositions": {"X": "C12", "R0": "C12"},
        "output": "C(24)H(2)O(1)",
        # 'mod_pos_info' : (0, [{'C': 12 }], 1, {'O' : 1, 'C' : 2, 'H' : 2}])
    },
    {"input": "R#Label:13C(6)15N(2):1", "output": "H(14)13C(6)15N(2)N(2)O(2)"},
    {
        "input": "R0XR1",
        "aa_compositions": {"R0": "C12", "X": "C13", "R1": "C14"},
        "output": "C(39)H(2)O(1)",
        "aa_pos_info": [
            {"pos": 1, "cc_mods": {"C": 12}},
            {"pos": 3, "cc_mods": {"C": 14}},
        ],
    },
    {
        "input": "RR#Label:13C(6)15N(2):1",
        "output": "C(6)H(26)13C(6)15N(2)N(6)O(3)",
        "aa_pos_info": [
            {"pos": 1, "cc_mods": {"C": 6, "O": 1, "N": 4, "H": 12}},
            {"pos": 2, "cc_mods": {"C": 6, "O": 1, "N": 4, "H": 12}},
        ],
        "mod_pos_info": [{"pos": 1, "cc_mods": {"C": -6, "13C": 6, "N": -2, "15N": 2}}],
    },
]

FORMULA_TESTS = [
    {"input": "+H2O2H2", "output": "H(4)O(2)"},
    {"input": "+H2O2H2-HO", "output": "H(3)O(1)"},
]

SUBTRACTION_TEST_SET = [
    # functionality tested yet ...
    {"in": "+H2O2H2", "sub": "H2", "out": {"H": 2, "O": 2}}
]

UNIMOD_FILE_TEST = [
    {
        "input": "K#SILAC TMT:1",
        "output": "C(8)H(34)13C(10)15N(1)N(3)O(4)",
        "mod_pos_info": [
            {
                "pos": 1,
                "cc_mods": {"13C": 10, "15N": 1, "C": 2, "H": 20, "N": 1, "O": 2},
            }
        ],
        "usermod": Path(__file__).parent.joinpath("usermod.xml"),
    }
]


def pepitde_with_usermod_test():
    for test_id, test_dict in enumerate(UNIMOD_FILE_TEST):
        if "aa_compositions" not in test_dict.keys():
            test_dict[
                "aa_compositions"
            ] = chemical_composition.chemical_composition_kb.aa_compositions

        cc = chemical_composition.ChemicalComposition(
            sequence=test_dict["input"],
            aa_compositions=test_dict["aa_compositions"],
            unimod_file_list=[test_dict["usermod"]],
            add_default_files=False
        )
        yield check_results, cc, test_dict


def pepitde_with_unimod_test():
    for test_id, test_dict in enumerate(PEPTIDE_TESTS):
        if "aa_compositions" not in test_dict.keys():
            test_dict[
                "aa_compositions"
            ] = chemical_composition.chemical_composition_kb.aa_compositions

        cc = chemical_composition.ChemicalComposition(
            sequence=test_dict["input"], aa_compositions=test_dict["aa_compositions"]
        )
        yield check_results, cc, test_dict


def formula_with_unimod_test():
    for test_id, test_dict in enumerate(FORMULA_TESTS):
        if "aa_compositions" not in test_dict.keys():
            test_dict[
                "aa_compositions"
            ] = chemical_composition.chemical_composition_kb.aa_compositions

        cc = chemical_composition.ChemicalComposition(
            formula=test_dict["input"], aa_compositions=test_dict["aa_compositions"]
        )
        yield check_results, cc, test_dict


def old_format_with_unimod_test():
    for test_id, test_dict in enumerate(PEPTIDE_TESTS + FORMULA_TESTS):
        if "aa_compositions" not in test_dict.keys():
            test_dict[
                "aa_compositions"
            ] = chemical_composition.chemical_composition_kb.aa_compositions

        cc = chemical_composition.ChemicalComposition(
            deprecated_format=test_dict["input"],
            aa_compositions=test_dict["aa_compositions"],
        )
        yield check_results, cc, test_dict


def check_results(cc, test_dict):
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
    print(aa, chemformula)
    print(cc)
    assert cc.hill_notation() == chemformula


if __name__ == "__main__":
    Simple_amino_acid_test()
