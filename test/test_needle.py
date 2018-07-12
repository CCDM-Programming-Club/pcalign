""" Tests for the needle submodule. """

import numpy as np
import pytest

from pcalign.seq import Seq

from pcalign.needle import construct_matrix
from pcalign.needle import init_matrix
from pcalign.needle import score
from pcalign.needle import movements
from pcalign.needle import backtrack
from pcalign.needle import align

seq1s = [
    "GCATGCU",
    "GCAT",
]

seq2s = [
    "GATTACA",
    "GATTACA",
]


matrices = [
    np.array([[ 0, -1, -2, -3, -4, -5, -6, -7],
              [-1,  1,  0, -1, -2, -3, -4, -5],
              [-2,  0,  0, -1, -2, -3, -2, -3],
              [-3, -1,  1,  0, -1, -1, -2, -1],
              [-4, -2,  0,  2,  1,  0, -1, -2],
              [-5, -3, -1,  1,  1,  0, -1, -2],
              [-6, -4, -2,  0,  0,  0,  1,  0],
              [-7, -5, -3, -1, -1, -1,  0,  0]]),
    np.array([[ 0, -1, -2, -3, -4, -5, -6, -7],
              [-1,  1,  0, -1, -2, -3, -4, -5],
              [-2,  0,  0, -1, -2, -3, -2, -3],
              [-3, -1,  1,  0, -1, -1, -2, -1],
              [-4, -2,  0,  2,  1,  0, -1, -2]])
    ]


def all_zero(matrix):
    """ Utility function to check if all elements in a matrix are zero """
    return (matrix == 0).all()


@pytest.mark.parametrize("seq1,seq2,nrows,ncols", list(
    zip(seq1s, seq2s, [len(i) + 1 for i in seq1s], [len(i) + 1 for i in seq2s])
))
def test_construct_matrix(seq1, seq2, nrows, ncols):
    actual = construct_matrix(seq1, seq2)

    # Correct dimensions
    assert len(actual) == nrows
    assert len(actual[0]) == ncols

    # Is square?
    assert all(len(actual[0]) == len(row) for row in actual)

    # Is all zero?
    assert all_zero(actual)
    return


def test_init_matrix():
    example = np.array([[0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0]])

    expected1 = np.array([[ 0, -1, -2, -3, -4],
                          [-1,  0,  0,  0,  0],
                          [-2,  0,  0,  0,  0],
                          [-3,  0,  0,  0,  0]])

    expected2 = np.array([[ 0, -2, -4, -6, -8],
                          [-2,  0,  0,  0,  0],
                          [-4,  0,  0,  0,  0],
                          [-6,  0,  0,  0,  0]])

    actual1 = init_matrix(example, -1)
    assert (actual1 == expected1).all()

    actual2 = init_matrix(example, -2)
    assert (actual2 == expected2).all()
    return


@pytest.mark.parametrize("seq1,seq2,matrix", list(
    zip(seq1s, seq2s, matrices)
))
def test_score(seq1, seq2, matrix):
    """ This isn't really a great test.

    For starters we don't check whether the best matches are actually the best.
    Other issues, the loop for runnning it is similar to what we'll eventually
    have to implement.
    It's late. Bear with me :)
    """

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            expected = matrix[i, j]

            actual, best = score(i, j, seq1, seq2, matrix, match_reward=1,
                                 mismatch_penalty=-1, indel_penalty=-1)
            print(i, j, actual, expected, best)
            assert actual == expected
    return


@pytest.mark.parametrize("seq1,seq2,matrix,expected", [
    (
        "GA",
        "GG",
        np.array([[ 0, -1, -2],
                  [-1,  0,  0],
                  [-2,  0,  0]]),
        {(1, 1): [(0, 0)], # Score 1
         (2, 1): [(1, 1)], # Score 0
         (1, 2): [(0, 1), (1, 1)], # Score 0
         (2, 2): [(1, 1)]} # Score -1
    ),
])
def test_movements(seq1, seq2, matrix, expected):
    """ NB i should really be testing switching out the scorer """
    actual = movements(seq1, seq2, matrix)

    assert len(actual) == len(expected)

    for ij in actual:
        assert actual[ij] == expected[ij]

    return


@pytest.mark.parametrize("seq1,seq2,movements,expected", [
    (
        "GA",
        "GG",
        {(1, 1): [(0, 0)], # Score 1
         (2, 1): [(1, 1)], # Score 0
         (1, 2): [(0, 1), (1, 1)], # Score 0
         (2, 2): [(1, 1)]}, # Score -1
        [("GA", "GG")]
    ),
    (
        "GATTACA",
        "GCATGCU",
        {(1, 1): [(0, 0)],
         (1, 2): [(1, 1)],
         (2, 3): [(1, 2)],
         (3, 3): [(2, 3)],
         (3, 4): [(2, 3)],
         (4, 4): [(3, 3), (3, 4)],
         (4, 5): [(3, 4)],
         (5, 5): [(4, 4), (4, 5)],
         (6, 6): [(5, 5)],
         (7, 7): [(6, 6)]},
        [("G-ATTACA", "GCATG-CU"),
         ("G-ATTACA", "GCA-TGCU"),
         ("G-ATTACA", "GCAT-GCU")]
    ),
    ])
def test_backtrack(seq1, seq2, movements, expected):
    alignments = backtrack(movements, seq1, seq2)

    assert len(alignments) == len(expected)
    for alignment in alignments:
        assert alignment in expected
    return

@pytest.mark.parametrize("seq1,seq2,expected", [
        (
            Seq(id="one", desc=None, seq="GATTACA"),
            Seq(id="two", desc="desc", seq="GCATGCU"),
            [
                Seq(id="a_0|one", desc=None, seq="G-ATTACA"),
                Seq(id="a_0|two", desc="desc", seq="GCATG-CU"),
                Seq(id="a_1|one", desc=None, seq="G-ATTACA"),
                Seq(id="a_1|two", desc="desc", seq="GCA-TGCU"),
                Seq(id="a_2|one", desc=None, seq="G-ATTACA"),
                Seq(id="a_2|two", desc="desc", seq="GCAT-GCU")
            ],
        )
    ])
def test_align(seq1, seq2, expected):
    actual = align(seq1, seq2, match_reward=1,
                   mismatch_penalty=-1, indel_penalty=-1,
                   seqid_template="a_{}")
    assert len(actual) == len(expected)

    expected_seqs = {s.seq for s in expected}
    expected_ids = {s.id for s in expected}

    for seq in actual:
        assert seq.id in expected_ids
        assert seq.seq in expected_seqs

    return
