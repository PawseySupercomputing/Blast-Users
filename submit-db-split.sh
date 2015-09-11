#!/bin/bash --login

#SBATCH --nodes=2
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

# Set as required...

export MY_BLAST="blastn -task blastn -num_descriptions 8 -num_alignments 0 -db nt"
export MY_BLASTDB="/group/data/blast"
export MY_BLAST_QUERY_FILES="*.fasta"

# Do not change...

export MY_NODES=${SLURM_NNODES}
export MY_BLAST_THREADS=48

aprun -n ${MY_NODES} -N 1 -d ${MY_BLAST_THREADS} -j 2 ./my-db-split.sh ${MY_BLAST_QUERY_FILES}

