import numpy
import itertools as it

import HTSeq as hts

# read data
c01 = hts.FastqReader("C2680001.fq", "phred")
c02 = hts.FastqReader("C2680002.fq", "phred")

#

for read in it.islice(c01, 10):
	print read

# 
