""" Code related to the needleman-wunsch alignment algorithm. """

import logging
from collections import defaultdict
from collections import deque
from functools import partial
from copy import deepcopy

import numpy as np

from pcalign import utils
from pcalign import seq

logger = logging.getLogger(__name__)


@utils.log(logger, logging.DEBUG)
def construct_matrix(seq1, seq2):
    """ Construct an x by y matrix filled with zeros.

    x and y are the lengths of the input sequences + 1.
    i.e. len(seq1) + 1

    Keyword arguments:
    seq1 -- A Seq or string. Will form the rows.
    seq2 -- As for seq 1. Will form the columns.

    Examples:

    >>> construct_matrix('ABC', 'DEFG')
    array([[0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0.],
           [0., 0., 0., 0., 0.]])
    """

    return np.zeros((len(seq1) + 1, len(seq2) + 1))


@utils.log(logger, logging.DEBUG)
def init_matrix(matrix, indel_penalty):
    """ The first row and column will decrease from the top left corner.

    e.g.

     0 -1 -2 -3
    -1  0  0  0
    -2  0  0  0
    -3  0  0  0
    -4  0  0  0

    Keyword arguments:
    matrix -- A 2 dimensional numpy array filled with zeros.
    indel_penalty -- The penalty for indels in the scoring system.

    Returns -- A matrix with the first row and column having decreasing numbers.

    Examples:

    >>> init_matrix(np.array([[0, 0, 0, 0, 0],
    ...                       [0, 0, 0, 0, 0],
    ...                       [0, 0, 0, 0, 0],
    ...                       [0, 0, 0, 0, 0]]), -1)
    array([[ 0, -1, -2, -3, -4],
           [-1,  0,  0,  0,  0],
           [-2,  0,  0,  0,  0],
           [-3,  0,  0,  0,  0]])
    """

    # Update the first column
    for i in range(1, matrix.shape[0]):
        matrix[i, 0] = indel_penalty * i

    # Update the first row
    for j in range(1, matrix.shape[1]):
        matrix[0, j] = indel_penalty * j
    return matrix


@utils.log(logger, logging.DEBUG)
def score(i, j, seq1, seq2, matrix, match_reward=1,
          mismatch_penalty=-1, indel_penalty=-1):
    """ Score the possible moves to a point in the matrix.

    Keyword arguments:
    i -- The row index
    j -- The column index
    seq1 -- A string or Seq object in the rows
    seq2 -- A string or Seq object in the columns
    match_reward -- The score to give base pair matches (1).
    mismatch_penalty -- The score to give base pair mismatches (-1).
    indel_penalty -- The score to give indels (-1).

    Returns:
    max_score -- The maximum score from a move to this point.
    max_indices -- A list of tuples containing the indices giving the max score.

    Examples:
    >>> score(1, 1, "G", "G", np.array([[0, -1], [-1, 0]]), 1, -1, -1)
    (1, [(0, 0)])
    >>> score(2, 1, "GA", "G", np.array([[0, -1], [-1, 1], [-2, 0]]), 1, -1, -1)
    (0, [(1, 1)])
    >>> score(1, 1, "A", "G", np.array([[1, 1], [0, 0]]), 1, -1, -1)
    (0, [(0, 0), (0, 1)])
    """

    # To reduce typing we'll write out the indices
    diag_ij = (i - 1, j - 1)
    top_ij = (i - 1, j)
    left_ij = (i, j - 1)

    # First check the diagonal score.
    # If the bases are a match at this point...
    if seq1[i - 1] == seq2[j - 1]:
        diag = matrix[diag_ij] + match_reward
    else:
        diag = matrix[diag_ij] + mismatch_penalty

    # Next get the indel scores (top and left elements of the matrix).
    top = matrix[top_ij] + indel_penalty
    left = matrix[left_ij] + indel_penalty

    # Collate the scores and indices into a list
    scores = [diag, top, left]
    indices = [diag_ij, top_ij, left_ij]

    # Find the max score, and get the indices of each move that gives that score.
    max_score = max(scores)
    max_indices = [ij for ij, s in zip(indices, scores) if s == max_score]

    return max_score, max_indices


@utils.log(logger, logging.DEBUG)
def movements(
        seq1,
        seq2,
        matrix,
        scorer=partial(score,
                       match_reward=1,
                       mismatch_penalty=-1,
                       indel_penalty=-1)
        ):
    """ Score the matrix and find the possible paths through it.

    Keyword arguments:
    seq1 -- A string or Seq object to compare in the rows.
    seq2 -- A string or Seq object to compare in the columns.
    matrix -- An initialised matrix.
    scorer -- A function that takes an index, the seqs and the matrix, and
              returns the max score and the paths that give that score.
              See section below.

    The "scorer" function has the signature:
    scorer(i: int, j: int, seq1: Seq, seq2: Seq, matrix: np.array) ->
                          (int, [(int, int)])

    The method of scoring is independent of the actual traversal of the matrix,
    so we could swap out the method for say a substitution matrix based method.
    It's good practice to keep your code modular, and providing functions as
    arguments can help with that a lot. Note the use of the `partial` function
    from `functools` to supply default arguments (used here for illustrative
    purposes, those were already the defaults).

    Returns:
    A dictionary keyed by tuples containing a list of tuples as the value.
    The key tuples correspond to indices in the matrix, and the list of tuples
    are the indices of previous cells where the movement had the maximum score.
    """

    memo = dict()

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score_ij, movements_ij = scorer(i, j, seq1, seq2, matrix)
            matrix[i, j] = score_ij
            memo[(i, j)] = movements_ij

    return memo


@utils.log(logger, logging.DEBUG)
def direction(child, parent, seq1, seq2):
    """ Find the alignment between seqs in the matrix for the child indices.

    This function takes two points of indices from a matrix, it decides what
    direction the movement between those points is (diagonal, vertical, or
    horizontal) and outputs the corresponding alignment (match, gap in seq2,
    gap in seq1) as a tuple of strings.

    Raises an error if an unsupported movement types is encountered, e.g.
    if i or j is less than x or y.

    Remember that the indices are for the matrix, not the sequences. So the
    bases returned for matrix column 7 will be for seq2 base 6.

    Keyword arguments:
    child -- A matrix index corresponding to the node closer to the bottom
        right of the matrix.
    parent -- A matrix index representing the node closer to the top left of
        matrix, where the alignment path has moved from to the child.
    seq1 -- The sequence in the rows of the alignment matrix.
    seq2 -- The sequence in the columns of the alignment matrix.

    Returns:
    A 2-tuple of aligned bases or gaps.
    """

    i, j = child
    x, y = parent

    if i > x and j > y: # diagonal
        return seq1[i - 1], seq2[j - 1]
    elif i > x and j == y: # vertical
        return seq1[i - 1], "-"
    elif j > y and i == x: # horizontal
        return "-", seq2[j - 1]
    else:
        # We should have a better error message usually, but this should never
        # happen if we used the movements and backtrack functions.
        raise ValueError("Illegal movement")
    return


@utils.log(logger, logging.DEBUG)
def backtrack(mvmt, seq1, seq2):
    """ Take the paths through the matrix and reconstruct all best alignments.

    This function takes a dictionary of lists of child matrix indices
    keyed by parent matrix indices. Starting at the bottom right of the matrix,
    we trace back through all of the highest scoring paths through the matrix.

    Note that a lot of example NW implementations simply choose one of the
    optimal alignments.

    The implementation of the function is based on a non-recursive breadth
    first tree traversal algorithm. We use a deque as the FIFO queue
    because it has similar characteristics to a doubly linked list.

    The recursive algorithm is much simpler (give it a go!), but recursion is
    computationally and memory expensive in python. In fact it performs so
    badly that python limits the recursion depth to 1000 to prevent a memory
    error called a stack overflow. Practically speaking, that means that we
    wouldn't be able to align sequences longer than 1000 bp.

    Keyword arguments:
    mvmt -- A dictionary keyed by 2-tuples containing a list of 2-tuples.
        Output from the `movements` function is intended input.
    seq1 -- The sequence aligned in the matrix rows.
    seq2 -- The sequence aligned in the matrix columns.

    Returns:
    A list of 2-tuples containing aligned sequences.
    """

    # The first index to get will always be the bottom right of the matrix
    idx = (len(seq1), len(seq2))

    # We keep track of the aligned paths in this dictionary.
    memo = defaultdict(list)
    # The last cell of the matrix has an empty sequence.
    memo[idx].append([])

    # The visited set tells us if we've been to a node before.
    # It helps us avoid redoing work, or duplicating results.
    visited = set()

    # We use a deque a bit like a linked list.
    # Essentially, it's a list that is computationally cheap to remove the
    # head from. For a normal list to 'pop' the first element (return
    # the value and remove it from the list) python would have to create a whole
    # new list in the background.
    queue = deque([idx])

    # As long as the queue has some elements we keep doing the loop.
    while len(queue) > 0:
        # Take the first queue element and remove it.
        idx = queue.popleft()

        logger.debug("Running %s", idx)
        logger.debug("Paths at %s %s", idx, memo[idx])

        # Skip the node if we've already processed it.
        if idx in visited:
            continue
        else:
            visited.add(idx)

        # Get the possible next nodes from the movements dict.
        # If the current node isn't in the mvmt dict, then we've either reached
        # (0, 0) or somethings gone horribly wrong. Here we assume it's (0, 0)
        # i.e. we have reached the top left of the matrix and all possible
        # alignments have been found, so we break the loop.
        try:
            parents = mvmt[idx]
        except KeyError:
            logger.debug("Breaking backtrack at %s", idx)
            break

        for parent in parents:
            # We push the parent node onto the end of the queue so that we
            # remember to process them later.
            queue.append(parent)

            # We get the alignment between the two bases.
            base1, base2 = direction(idx, parent, seq1, seq2)
            logger.debug("idx %s parent %s alignment %s-%s",
                         idx, parent, base1, base2)

            # And we append this alignment to all of the paths that end at
            # this cell, and add those updated paths to the parent nodes.
            for seq in memo[idx]:
                # We need to deepcopy because lists are mutable pointers, and
                # we don't want weird side effects.
                this_seq = deepcopy(seq)
                # Sequences are represented as a list of tuples, where the tuple
                # has a single base from seq1, and seq2.
                this_seq.append((base1, base2))
                memo[parent].append(this_seq)

    # Finally we'll restructure the output from a list of list of tuples, to
    # a list of tuples of sequences i.e. [("ATGC", "A-GC"), ...]
    output = []
    for path in memo[(0, 0)]:
        # Note that we have to reverse the sequences because backtracking
        # goes backwards through the matrix.
        s1 = "".join(tup[0] for tup in reversed(path))
        s2 = "".join(tup[1] for tup in reversed(path))
        output.append((s1, s2))

    return output


@utils.log(logger, logging.DEBUG)
def align(seq1, seq2, match_reward, mismatch_penalty,
          indel_penalty, seqid_template="alignment_{:02d}"):
    """ Run the needleman-wunsch algorithm and return all best alignments.

    Keyword arguments:
    seq1 -- A Seq object to align.
    seq2 -- Another Seq object to align.
    match_reward -- The score when bases match.
    mismatch_penalty -- The score when bases don't match, note that negative
        numbers are required to actually penalise.
    indel_penalty -- The score for indels. Again negative numbers required to
        penalise this type of match. This model assumes no difference between
        indel start and extension.
    seqid_template -- A string using pythons formatter syntax to use as names
        for the output sequences. E.g. default "alignment_{:02d}" would give
        "alignment_01", "alignment_02". A single parameter is passed to format
        a number indicating the alignment number.

    Returns:
    A list of 2-tuples of Seq objects, where the tuples represent alignment.
    """

    matrix = init_matrix(construct_matrix(seq1.seq, seq2.seq),
                         indel_penalty)

    # Define a helper function using partial to use custom score penalties
    # and rewards.
    scorer = partial(score,
                     match_reward=match_reward,
                     mismatch_penalty=mismatch_penalty,
                     indel_penalty=indel_penalty)

    mvmts = movements(seq1.seq, seq2.seq, matrix, scorer)
    alignments = backtrack(mvmts, seq1.seq, seq2.seq)

    # Reformat the output alignment strings as Seq objects.
    output = []
    for i, (ali1, ali2) in enumerate(alignments):
        name = seqid_template.format(i)

        out_seq1 = seq.Seq(id="{}|{}".format(name, seq1.id),
                           desc=seq1.desc,
                           seq=ali1)

        out_seq2 = seq.Seq(id="{}|{}".format(name, seq2.id),
                           desc=seq2.desc,
                           seq=ali2)

        output.extend([out_seq1, out_seq2])

    return output
