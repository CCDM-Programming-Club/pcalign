"""
"""

__version__ = "0.0.0"

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s %(name)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

import sys
import argparse

from pcalign import utils
from pcalign.seq import Seq
from pcalign.needle import align

logger = logging.getLogger("pcalign")

@utils.log(logger, logging.DEBUG)
def cli(prog, args):
    """ Process command line arguments.

    Often this is simply put in the main function or the
    if __name__ == "__main__" block, but keeping it as a function allows
    testing if we had more complex command line interfaces.

    Keyword arguments:
    prog -- The name of the program. Usually this will be the first element of
        the sys.argv list.
    args -- The command line arguments for the program. Usually this will be
        a slice from i.e sys.argv[1:].

    Returns:
    An argparse Args object containing the parsed args.

    Raises:
    Standard argparse exceptions if the CLI conditions are not met.
    Required parameters, type constraints etc.
    """

    parser = argparse.ArgumentParser(
        prog=prog,
        description="Project description."
        )

    # For file IO we're using some special syntax from argparse.
    # Rather than having to open and close the files, handle exceptions, or
    # check that the files exist, we just let argparse do all of that for us.
    # It also enables the standard alias for stdin/stdout '-' by default.
    parser.add_argument(
        "-1", "--seq1",
        required=True,
        type=argparse.FileType('r'),
        help="FASTA file 1. Use '-' for stdin.",
        )

    parser.add_argument(
        "-2", "--seq2",
        required=True,
        type=argparse.FileType('r'),
        help="FASTA file 2. Use '-' for stdin.",
        )

    parser.add_argument(
        "-o", "--output",
        type=argparse.FileType('w'),
        default=sys.stdout,
        help="Output fasta file path. Default stdout.",
        )

    parser.add_argument(
        "-m", "--match-score",
        dest="match_score",
        type=float,
        default=1,
        help="The score to assign for matching base pairs. Default = 1."
        )

    parser.add_argument(
        "-n", "--mismatch-score",
        dest="mismatch_score",
        type=float,
        default=-1,
        help="The score to assign for mis-matching base pairs. Default = -1."
        )

    parser.add_argument(
        "-i", "--indel-score",
        dest="indel_score",
        type=float,
        default=-1,
        help="The score to assign for a gap insertion. Default = -1."
        )

    # This one is interesting, a user can provide this flag multiple times,
    # and the number of times it's specified is counted, giving the value.
    # I've only ever seen it used to control the "verbosity" of running updates
    # or logging.
    parser.add_argument(
        "-v", "--verbose",
        help=("Print progress updates to stdout. Use twice (i.e. -vv or -v -v)"
              "to show debug output."),
        action="count",
        default=0
        )

    # This is also a special action, that just prints the software name and
    # version, then exits.
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__
        )
    return parser.parse_args(args)


def main():
    """ The command line version of pcalign.

    This is the entrypoint specified in setup.py as `pcalign`.
    """

    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    if args.verbose > 1:
        log_level = logging.DEBUG
    elif args.verbose > 0:
        log_level = logging.INFO
    else:
        log_level = logging.WARNING

    logger.setLevel(log_level)
    logger.info("Running pcalign")
    logger.info("Using parameters:")
    logger.info("- seq1 = %s", args.seq1.name)
    logger.info("- seq2 = %s", args.seq2.name)
    logger.info("- output = %s", args.output.name)
    logger.info("- match score = %d", args.match_score)
    logger.info("- mismatch score = %d", args.mismatch_score)
    logger.info("- indel_score = %d", args.indel_score)


    logger.info("Parsing sequences")
    try:
        # Use this seq1_loaded thing to see which file fails, if one does.
        seq1_loaded = False

        # next() will return the next element of an iterator.
        # I'm using it to pop off the first sequence.
        # If there were more sequences, we simply ignore them
        seq1 = Seq.read(args.seq1)
        seq1_loaded = True

        seq2 = Seq.read(args.seq2)
    except ValueError as e:
        # next will yield this exception if the iterator is empty.
        # This will happen if there are no fasta headers in the file.
        # To handle this user error, we log the problem and stop the program.
        # It's common to handle errors like this when we're dealing with user
        # input.
        if seq1_loaded:
            fname = args.seq2.name
        else:
            fname = args.seq1.name

        logger.error("Parsing fasta files failed.")
        logger.error("Offending file is '%s'", fname)
        logger.error(e)
        sys.exit(1)


    logger.info("Aligning sequences")
    aligned = align(seq1,
                    seq2,
                    match_reward=args.match_score,
                    mismatch_penalty=args.mismatch_score,
                    indel_penalty=args.indel_score)

    logger.info("Writing alignments")
    for seq in aligned:
        print(seq, file=args.output)

    return
