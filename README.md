snoop
==

A collection of Python scripts for exploratory analyses of next-generation sequencing data stored in a [multi-string BWT (msBWT)](https://pypi.python.org/pypi/msbwt) data structure.

Dependencies
--
* Python >= 2.7
* Cython (http://cython.org/)
* numpy (http://www.numpy.org/)
* pysam (http://code.google.com/p/pysam/)
* bedtools (http://bedtools.readthedocs.org/en/latest/)
* pybedtools (http://pythonhosted.org/pybedtools/)
* pyfasta (https://pypi.python.org/pypi/pyfasta/)
* msbwt (https://pypi.python.org/pypi/msbwt)

Scripts
--
* `fa2bwt.py`: build a msBWT from one or more fasta files (ie. a reference genome)
* `bwt_query.py`: query a msBWT with an arbitrary string; return matching strings in a pretty alignment
* `kmer_profile.py`: compute an approximation to the Jensen-Shannon divergence between a pair of msBWTs, possibly with downsampling to save time
* `count_reference_kmers.py`: given one or more msBWTs, a reference sequence, and a BED file specifying intervals on that sequence, query the msBWT(s) with (arbitrarily-spaced) *k*-mers from the reference and return the counts
* `assemble_me_this.py`: given a msBWT and a seed sequence, perform targeted de novo assembly, trying to extend seed as far as possible before things get messy
* `find_insertions.py`: given a bed file of k-mer sequences and counts along some reference (probably from `count_reference_kmers.py`), attempt to find transposable-element insertions which are absent in that reference
