#!/usr/bin/env mpibash

#  This script performs a BLAST search specifed by the environment
#  variable MY_BLAST. The number of threads should be specified by
#  MY_BLAST_THREADS. This is described in the accompanying
#  submit-query-split.sh
#
#  The command line argument is the list of query files (.fasta or
#  .fa extension expected).
#
#  The implementation uses mpibash to control parallelism, and the
#  associated circle queue API to provide work-sharing for the
#  different files.

enable -f mpibash.so mpi_init
enable -f circlebash.so circle_init

function main() {

    mpi_init
    mpi_comm_rank rank
    mpi_comm_size size

    circle_init
    circle_set_options split_equal
    circle_enable_logging error

    # Register the circle queue functions and begin

    circle_cb_create  files_add
    circle_cb_process file_blast

    circle_begin

    circle_finalize
    mpi_finalize

}

###############################################################################
#
#  files_add
#
#  This function adds all the files on the command line to the queue
#  for processing.
#
#  As this occurs only on rank 0, we take the opportunity to provide
#  some information to stdout.
#
#  The function is registered as a callback by circle_cb_create
#
###############################################################################

function files_add() {

    local filename

    echo ""
    echo "Blast command: ${MY_BLAST}"

    # Rank 0 enqueues all of the files on the command line

    echo ""
    echo "Starting the circle queue with ${size} MPI tasks"
    echo ""

    for filename in "${BASH_ARGV[@]}"; do
	circle_enqueue "${filename%/}"
    done

    return
}

###############################################################################
#
#  file_blast
#
#  This function processes a file from the queue by running the requested
#  Blast search
#
#  Registered as a callback by circle_cb_process
#
###############################################################################

function file_blast() {

    local stub
    local filename_input
    local filename_output

    circle_dequeue filename_input
    
    # Generate an output file name (replace .fasta or .fa by ".blast")

    stub=`echo $filename_input | sed 's/.f\(ast\)*a$//'`
    filename_output=${stub}.blast

    STARTTIME=$(date +%s)

    ${MY_BLAST} -num_threads ${MY_BLAST_THREADS} -query ${filename_input} -out ${filename_output}

    ENDTIME=$(date +%s)

    # Report

    printf "Rank %2d: %s -> %s   \t%4d s\n" "${rank}" "${filename_input}" "${filename_output}" "$(($ENDTIME-$STARTTIME))"

    return
}

# Execute and finish

main "$@"
