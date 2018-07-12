# PCalign

This is an implementation of the Needleman-Wunsch alignment algorithm based on 
the description in this wikipedia article <https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm>.

It is mostly intended as a learning exercise for some ways of structuring
python programs and some of the tooling that goes along with that.


Briefly, a NW aligner takes two sequences and finds all optimal alignments that
cover all of both sequences (a global alignment). Optimality is decided based
on some scoring system, and it is possible to have multiple alignments with the
best score. This script implements the simple scoring method, and returns all
alignments with the best score.

Quickstart:

```bash
# Clone this repo
# Make sure you have a github account and have your ssh keys setup first!
git clone git@github.com:CCDM-Programming-Club/pcalign.git && cd pcalign

# Setup a virtual environment
python3 -m venv env
source env/bin/activate

# On Windows the second line may be:
# $ source env/Scripts/activate

# Install the package
pip install -e .

# Run the program
pcalign -1 example/example1.fasta -2 example/example2.fasta
```

You should also be able to import the modules from within python, e.g.

```python
from pcalign.seq import Seq
from pcalign.needle import align

with open("example/example1.fasta") as handle:
    seq1 = Seq.read(handle)

with open("example/example2.fasta") as handle:
    seq2 = Seq.read(handle)

alignments = align(seq1, seq2, 1, -1, -1)
```

Have a look at the source code for more on what these do.

To run the tests, do:

```bash
# within the virtualenv from earlier

python setup.py test
```

## Packaging python projects.

Packaging python projects up basically boils down to some folder structure,
`__init__.py` files, and a special file called `setup.py`.

Creating your own module is as simple as putting your script in a folder, and
putting an `__init__.py` file in the same folder. The `__init__.py` file can
just be empty because it just tells python that this folder is something you
want to be able to import. In this code the `__init__.py` file contains the
code for running our command line version of the progam.

Mimicking this project should get you pretty far, but for more information
check out the official [python packaging user guide](https://packaging.python.org)
which contains, a great basic tutorial <https://packaging.python.org/tutorials/packaging-projects/>.

One nice aspect of pip is that you can install the package and still edit it, 
and use those changes.

The command:

```python
pip install -e .
```

Tells python to install the current directory (i.e. `.`) as a package, but
allow us to use the edits in the installed version `-e`.
So you could open a jupyter notebook and play with the package, then edit the
package, and to use those new edits, you just need to reload the package or
restart the notebook. Same goes for scripts.
N.B. pip requires the `setup.py` file to work properly.


## Packaging command line scripts.

You might have noticed in the quick start that we didn't have to specify a path
to the `pcalign` script. This is because pip put the script inside `./env/bin/` for us!
We're using a parameter called `entrypoints` in the `setup.py` to tell pip to do that for us.

To learn more about this, check out <http://python-packaging.readthedocs.io/en/latest/command-line-scripts.html>.


## Unit-tests

Testing your code seems like a waste of time, until you come to use it again in
six months, or you realise that your output was all wrong.

Save yourself and the community from this pain with unit-tests.
Testing is a big concept and I wouldn't claim to be an expert at it, so
you might need to do some googling. The concepty are not python specific so
pick something that explains the ideas well for you. There are a lot of good
presentations on youtube too.

The basic idea of unit-tests is that you test the small functions and methods
that you use to create the big functions or your overall script. If you've
got those core bits right, then the big bits that piece those small bits
together are more likely to work. It's much easier to test the small functions
because they don't take long to run, and it's easier to figure out how those
functions might come up with the right answer and write a test to make sure
you handle that in your code.


For python, I find [pytest](https://docs.pytest.org/en/latest/) to be the easiest testing library.
You can have a look at the test files in the `test` folder.

Note that I'm also testing the documentation within my code!
There is a standard called doctest, which specifies how code in docstrings should
be written and tested <https://docs.python.org/3.7/library/doctest.html>.
Pytest can then use a plugin to run doctest for you.
Now both your code and documentation is going to be checked!


## Logging

Sometimes things aren't quite going right and you can't figure out where it's going wrong.
Rather than putting `print` calls everywhere or trying to use a
[debugger](https://docs.python.org/3/library/pdb.html), logging can be nice
because you can leave the code there and just switch that verbose logging on when
you need it.

The python logging module does this for us <https://docs.python.org/3/howto/logging.html#logging-basic-tutorial>.
It's a fairly old module, so the interface is slightly different to what we might be used to (e.g. string formatting).

In our script I used logging to write output to the screen instead of print statements.
We use INFO level logging for regular information (e.g. "aligning sequences", "finished", etc), 
and we use DEBUG level logging to help us understand what's happening where.

We control the level of the logging output (e.g. do we want to show DEBUG, INFO, or neither?),
using a command line flag `-v`.

I've also used a python decorator function to automatically log the input and output
of functions.
This is a bit complicated, but just know that the funny "@"s all over the place just do that.

You can learn more about decorators here: <https://realpython.com/primer-on-python-decorators/>, and <https://www.learnpython.org/en/Decorators>
My logging decorator is in the `utils.py` file.
