#!/bin/bash --login

#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --account=
#SBATCH --export=none

module load mpibash
module load blast+

# MY_BLAST contains a standard Blast+ application command line except:
#
# -query       should not be specified
#              instead use a list of files specified by MY_BLAST_QUERY_FILES
#              which is the list of FASTA inputs to be processed; this allows
#              standard shell file name expansion.
#
# -out         should not be specified
#              for each input file name, the output file name is constructed
#              by replacing the .fasta (or .fa) extension with .blast
#
# -num_threads should not be specified
#              instead specify MY_BLAST_THREADS

export MY_BLAST="tblastn -task tblastn -num_descriptions 16 -num_alignments 1 -db my-data-base"
export MY_BLAST_QUERY_FILES="*.fasta"

# The following can be used unchanged.

export MY_BLAST_THREADS=2

aprun -n 24 -N 24 -S 12 -d ${MY_BLAST_THREADS} -j 2 ./my-query-split.sh \
    ${MY_BLAST_QUERY_FILES}
