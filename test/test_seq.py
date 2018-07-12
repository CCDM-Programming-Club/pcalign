""" Unittests for the Seq class. """

from types import GeneratorType

import pytest
from pcalign.seq import Seq

# Define some test data that we'll be using

seqs = [
""">seq1 testing
TTTCCGGGGCACATAATCTTCAGCCGGGCGC""",
""">seq2 testing2
ACTAAGTAGTCTTTTTGAGGTCGTTAACTCTTATAAAGCGGCGCAGCATACCTCCCGAGA
CTATAGTTTTTCTCAATGCTGAACGCCTCATGGCTTGCCGGGCTCAATGCTGTAATCTGT
CTCGGTTCCTGTATACTAGCCGGTACTCCCCAGTTAATTCGACTCGTTGTTTCTCTGTAT
GTCTCCGATACATCCTAATATAATGTCCCCATGCTTACGCCTATAAAATCGCAATACTGT
CTAAGGGAGGTCACTTAATTGTGAAGAGAGCCTAGACAGCGTTCGATTTAGAGCGTCCGT
ACCAGGATCTTCTATCGGGCTCTGTGATGATTATAGCTATCGCTGACCGCCGGCTCGTCC
TAGCGTTTAATACGGCGTACCGACCACTAGGGGGGAGGAAGTAGTTACCATTATCATCCA
TATCGATTAAGGTGCGTTTCAAGCGTTGACATTAAAGCCGAAACGCAAGGGCAATGCAAG
TTCTGGTGTAATCATGGAAGTTAATGTTGCTCGCGTGTGTCAATGCTGGGTACAGGAGAA
TAGTGTGTATGCGTGTCAGATCCCCCAGCCGCAAAGTCCCCTTCAGTCGTGCCAAGGCGG
GAAATTCCAACTCTCGTGTCCCCATTCCCGCGCCTTGCTAAAGACATTACTAGATACGCT
TGCTTACGGAGCTACGAAACATGTGTGGCAACTCTCCAGTGCGCAGCGCCCCATAGGTTA
GGCACGGAGACAGTTCGCGTACCAGGTTCTAAATTGAGTAGGTTCGCCATGAGCAGTTAC
CACATACTACCTTGTCTGACACAGGTGACATACCGGCGGGCTGAGTATTGTGATCATGGT
GCGTATATATTGTTTCCCGTCCGTCCCCCCGGTGCACGAACTATCATCTAGCCGGCTATT
TCGTTCAGTTAGCGTAGCTCGTTGCAGAGAAGTGAATTACGTTAAGGGGATGAGCGCCCA
GTCCTCGCCCTCGCCGCTGCCATGGATATAGCAACGTT""",
""">seq3
ATCCAGCT"""
]

expected = [
    Seq(id="seq1", desc="testing", seq="TTTCCGGGGCACATAATCTTCAGCCGGGCGC"),
    Seq(id="seq2",
        desc="testing2",
        seq=("ACTAAGTAGTCTTTTTGAGGTCGTTAACTCTTATAAAGCGGCGCAGCATACCTCCCGAGA"
             "CTATAGTTTTTCTCAATGCTGAACGCCTCATGGCTTGCCGGGCTCAATGCTGTAATCTGT"
             "CTCGGTTCCTGTATACTAGCCGGTACTCCCCAGTTAATTCGACTCGTTGTTTCTCTGTAT"
             "GTCTCCGATACATCCTAATATAATGTCCCCATGCTTACGCCTATAAAATCGCAATACTGT"
             "CTAAGGGAGGTCACTTAATTGTGAAGAGAGCCTAGACAGCGTTCGATTTAGAGCGTCCGT"
             "ACCAGGATCTTCTATCGGGCTCTGTGATGATTATAGCTATCGCTGACCGCCGGCTCGTCC"
             "TAGCGTTTAATACGGCGTACCGACCACTAGGGGGGAGGAAGTAGTTACCATTATCATCCA"
             "TATCGATTAAGGTGCGTTTCAAGCGTTGACATTAAAGCCGAAACGCAAGGGCAATGCAAG"
             "TTCTGGTGTAATCATGGAAGTTAATGTTGCTCGCGTGTGTCAATGCTGGGTACAGGAGAA"
             "TAGTGTGTATGCGTGTCAGATCCCCCAGCCGCAAAGTCCCCTTCAGTCGTGCCAAGGCGG"
             "GAAATTCCAACTCTCGTGTCCCCATTCCCGCGCCTTGCTAAAGACATTACTAGATACGCT"
             "TGCTTACGGAGCTACGAAACATGTGTGGCAACTCTCCAGTGCGCAGCGCCCCATAGGTTA"
             "GGCACGGAGACAGTTCGCGTACCAGGTTCTAAATTGAGTAGGTTCGCCATGAGCAGTTAC"
             "CACATACTACCTTGTCTGACACAGGTGACATACCGGCGGGCTGAGTATTGTGATCATGGT"
             "GCGTATATATTGTTTCCCGTCCGTCCCCCCGGTGCACGAACTATCATCTAGCCGGCTATT"
             "TCGTTCAGTTAGCGTAGCTCGTTGCAGAGAAGTGAATTACGTTAAGGGGATGAGCGCCCA"
             "GTCCTCGCCCTCGCCGCTGCCATGGATATAGCAACGTT")),
    Seq(id="seq3", desc=None, seq="ATCCAGCT")
    ]

# Define the actual tests

@pytest.mark.parametrize("fasta,expected",
    list(zip([seq.split("\n") for seq in seqs], expected))
)
def test_Seq_read(fasta, expected):
    actual = Seq.read(fasta)
    assert actual.id == expected.id
    assert actual.desc == expected.desc
    assert actual.seq == expected.seq
    return


@pytest.mark.parametrize("fasta,expected", [
    ("\n".join(seqs).split("\n"), expected)
])
def test_Seq_parse(fasta, expected):
    actual = Seq.parse(fasta)

    # We want parse to return a generator over the sequences.
    assert isinstance(actual, GeneratorType)

    # Just this time, get a list so we can compare lenths.
    actual = list(actual)

    assert len(actual) == len(expected)

    for act, exp in zip(actual, expected):
        assert act.id == exp.id
        assert act.desc == exp.desc
        assert act.seq == exp.seq
    return

@pytest.mark.parametrize("expected,obj",
    list(zip(seqs, expected))
)
def test_Seq_str(expected, obj):
    """ Note that this is 'reciprocal' to the read method """

    actual = str(obj)
    assert actual == expected
    return


def test_Seq_len():
    """ Simple wrapper function, so gets a simple test. """
    seq = Seq(id="test", desc=None, seq="four")
    assert len(seq) == 4
    return


def test_Seq_getitem():
    """ Simple wrapper function, so gets a simple test. """
    seq = Seq(id="test", desc=None, seq="aaabbbbbaaa")
    assert seq[3:].seq == "bbbbbaaa"
    assert seq[1].seq == "a"
    return


@pytest.mark.parametrize("line,exp_id,exp_desc", [
    (">id desc", "id", "desc"),
    (">id", "id", None),
    ])
def test_Seq__split_id_line(line, exp_id, exp_desc):
    """ Are we splitting the fasta header lines correctly? """
    act_id, act_desc = Seq._split_id_line(line)

    assert act_id == exp_id
    assert act_desc == exp_desc
    return
