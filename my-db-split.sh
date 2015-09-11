#!/usr/bin/env mpibash

enable -f mpibash.so mpi_init

#  This script performs a simple Blast search on the input file
#  (or files) specified on the command line.
#
#  It is required that the script be executed on one MPI task per
#  node. Parallelism at the node level is provided by the threaded
#  model implemented in blast applications.
#
#  The database is staged to $TMP (/tmp is in memory on Magnus).
#  This means the database is limited to a size of 32GB in total.
#
#  For larger databases, more than one node must be used. In this
#  case, the results from partial searches must be recombined.

#  Search parameters originating in the accompanying submission
#  script are indicated in capital letters.

MY_BLAST_FORMATTER="blast_formatter"
MY_TMP="/tmp"

# For additional diagnostic output
my_verbose=""

function main() {

    local f

    mpi_init
    mpi_comm_rank rank
    mpi_comm_size size

    files="${BASH_ARGV[@]}"
    nfiles=${#BASH_ARGV[@]}

    parse_options

    export BLASTDB=${MY_BLASTDB}
    db=${MY_BLASTDB}/${my_db}

    if [[ $rank == 0 ]]; then
	echo ""
	echo "Blast with               ${MY_BLAST}"
	echo "Database name is         ${my_db}" 
	echo "Blast+ BLASTDB           ${BLASTDB}"
	echo "Path to database is      ${db}"
	echo "Number of query files    ${nfiles}"
	echo "Staging location         ${MY_TMP}"
	echo "Using max_target_seqs    ${my_max_target_seqs}"
	echo "Number of threads        ${MY_BLAST_THREADS}"
    fi

    db_stage

    while [ $nfiles -gt 0 ]; do
	((nfiles-=1))
	f="${BASH_ARGV[$nfiles]}"
	blast_process_input_file $f
    done

    mpi_barrier

    # Now everyone has finished remove the database
    # ...and finalise.

    db_unstage

    mpi_finalize

    return 0
}

###############################################################################
#
#  verbose
#
#  Utility to print a message
#
###############################################################################

function verbose() {

    if [[ ${my_verbose} ]]; then
	echo $1
    fi

    return
}

###############################################################################
#
#  parse_options
#
#  We need to parse MY_BLAST to separate out some of the options.
#
#  This will not parse e.g., '-outfmt "6 tag1 tag2"' correctly.
#  Instead use -outfmt 11 and reformat.
#
###############################################################################

function parse_options() {

    read -a my_blast_args <<< "${MY_BLAST}"
    nargs=${#my_blast_args[@]}

    # Defaults
    my_num_descriptions=500
    my_num_alignments=500
    my_max_target_seqs=500
    my_outfmt=0
    my_db="null"

    n=0

    verbose "parse_options: ${MY_BLAST}"
    verbose "parse options: argc = ${nargs}"

    while (( $n < $nargs )); do

	key="${my_blast_args[$n]}"
	((n+=1))

	case $key in
	    -num_descriptions)
		my_num_descriptions="${my_blast_args[$n]}"
		((n+=1))
		verbose "option: -num_descriptions = ${my_num_descriptions}"
		;;
	    -num_alignments)
		my_num_alignments="${my_blast_args[$n]}"
		((n+=1))
		verbose "option: -num_alignments = ${my_num_alignments}"
		;;
	    -max_target_seqs)
		my_max_target_seqs="${my_blast_args[$n]}"
		((n+=1))
		verbose "option: -max_target_seqs ${my_max_target_seqs}"
		;;
	    -outfmt)
		my_outfmt="${my_blast_args[$n]}"
		((n+=1))
		verbose "option: -outfmt = ${my_outfmt}"
		;;
	    -db)
		my_db="${my_blast_args[$n]}"
		((n+=1))
		verbose "option: -db = ${my_db}"
		;;
	    *)
		# Unknown option
		;;
	esac
    done

    if (( $my_outfmt <= 4 )); then
	# We may need ASN output with max(num_alignments, num_descriptions)
	my_max_target_seqs="${my_num_descriptions}"
	if (( $my_num_alignments > $my_max_target_seqs )); then
	    my_max_target_seqs="${my_num_alignments}"
	fi
    fi

    return
}

###############################################################################
#
#  function db_stage
#
#  is responsible for copying the database (all the files associated with
#  the indexed BLAST database) to ${MY_TMP}
#
#  It is assumed that one MPI task is running per node, so that one task
#  has available /tmp on the node.
#
#  If more than one MPI task is present, this function will attempt to
#  distribute the database volumes between the available nodes. If
#  there are more MPI tasks than data base volumes, this will cause
#  a failure in tasks that have no database volumes.
#
#  Note only the database files required on a particular node are copied;
#  we do not want to copy all the database files on all nodes.
#
#  A more official way to create the new alias file would be to run
#  blastdb_aliastool -db list -dbtype ... -out ... -title
#
###############################################################################

function db_stage() {

    local alias_ext;

    # Look for the alias file, either .nal or .pal
    # If no alias file is found, one could issue a warning, or try
    # to count the volumes in some other way, but here we just
    # default to trying to copy all the files with base name

    alias_ext=""
    [[ -e ${db}.nal ]] && alias_ext="nal"
    [[ -e ${db}.pal ]] && alias_ext="pal"

    # Stage database

    STARTTIME=$(date +%s)

    if [[ $size == 1 || $alias_ext == "" ]]; then
	# Single node; just copy all files
	cp ${db}.* ${MY_TMP}
    else

	# Multiple nodes: distribute round-robin

        # Count the total number of database fragments "ndb"
	# and record the size of the data base.
        # Make up the list of fragments to go in the .nal or .pal

	var=`grep DBLIST ${db}.${alias_ext}`
	ndb=`grep -o "${my_db}" <<<"${var}" | grep -c .`

	var=`grep LENGTH ${db}.${alias_ext}`
	my_dbsize=`echo $var | sed 's/LENGTH//'`

	list="DBLIST"

	for (( index=$rank; index<${ndb}; index+=$size )); do
	    file=`printf "%2.2d" ${index}`
	    list=`printf "%s \"%s.%s\"" "${list}" "${my_db}" "${file}"`
            # Note this is .${file}. to obtain exact matches with
	    # numerical part of the file name
            cp ${db}.${file}.* ${MY_TMP}
	done

	sed "s/DBLIST.*/${list}/" "${db}.${alias_ext}" > ${MY_TMP}/${my_db}.${alias_ext}

	printf -v var "Rank %2.2d %s\n" "${rank}" "${list}"
	verbose "$var"
    fi

    mpi_barrier
    ENDTIME=$(date +%s)

    if [[ $rank == 0 ]]; then
	echo
	echo "-> `date`"
	echo "-> Time to stage ${my_db} to ${MY_TMP} $(($ENDTIME-$STARTTIME)) s"
    fi

    return 0
}

###############################################################################
#
#  db_unstage
#
#  At the end, completely clear ${MY_TMP}. Caller must arrange
#  synchronistation, if required beforehand.
#
###############################################################################

function db_unstage() {

    if [[ $rank == 0 ]]; then
	echo ""
	echo "Unstage database..."
    fi

    rm -f ${MY_TMP}/*

    if [[ $rank == 0 ]]; then
	echo "Finished normally."
    fi

    return 0
}

###############################################################################
#
#  blast_process_input_file
#
#  The first argument is the file to blast; this is a switch for
#  single task (no data-base decomposition) or multiple task versions.
#
###############################################################################

function blast_process_input_file() {

    if [[ $size == 1 ]]; then
	blast_process_input_file_serial $1
    else
	blast_process_input_file_parallel $1
    fi

    return 0
}

###############################################################################
#
#  blast_process_input_file_serial
#
#  The entire database is accommodated on one node.
#
#  Just run the blast command using the staged data base.
#
###############################################################################

function blast_process_input_file_serial() {

    local f
    local stub
    local output

    STARTTIME=$(date +%s)

    f="$1"
    stub=`echo $f | sed 's/.f\(ast\)*a$//'` 
    output=${stub}.blast

    echo ""
    echo "-> Input  ${f}"
    echo "-> Output ${output}"

    export BLASTDB=${MY_TMP}

    ${MY_BLAST} -num_threads ${MY_BLAST_THREADS} -query ${f} -out ${output}

    ENDTIME=$(date +%s)

    echo "-> Time to blast $f   $(($ENDTIME-$STARTTIME)) s"

    return 0
}

###############################################################################
#
#  blast_process_input_file_parallel
#
#  This is the data-base splitting version.
#
#
#  To avoid significant, and error-prone, manipulation of parallel output,
#  I use the following scheme:
#
#  1. An initial, parallel, scan of the decomposed data base to identify
#     relevant gis for the query
#  2. Synchronise and produce one global gilist file on root
#  3. re-blast the original database using -gilist
#
###############################################################################

function blast_process_input_file_parallel() {

    local f
    local stub
    local gilist
    local gilist_all

    STARTTIME=$(date +%s)

    f="$1"
    stub=`echo $f | sed 's/.f\(ast\)*a$//'` 

    if [[ ${rank} == 0 ]]; then
	echo ""
	echo "-> Input  ${f}"
	echo "-> Output ${stub}.blast"
    fi

    # Run a reduced command in parallel to identify gis

    gilist="./${f}-gilist-${rank}.msk"
    export BLASTDB=${MY_TMP}

    my_blast=`echo ${MY_BLAST} | sed 's/-num_descriptions [0-9][0-9]*//'`
    my_blast=`echo ${my_blast} | sed 's/-num_alignments [0-9][0-9]*//'`
    my_blast=`echo ${my_blast} | sed 's/-max_target_seqs [0-9][0-9]*//'`
    my_blast=`echo ${my_blast} | sed 's/-outfmt [0-9][0-9]*//'`

    verbose "blast: ${my_blast}"

    ${my_blast} -num_threads ${MY_BLAST_THREADS} -query ${f} \
	-max_target_seqs ${my_max_target_seqs} -outfmt "6 sgi" -out ${gilist}

    # synchronise

    mpi_barrier

    if [[ ${rank} == 0 ]]; then
	ENDTIME=$(date +%s)
	echo "-> Time to blast $f   $(($ENDTIME-$STARTTIME)) s"
    fi

    # Now combine the results.

    if [[ ${rank} == 0 ]]; then

	STARTTIME=$(date +%s)

	gilist_all="./${f}-gilist-global.msk"
	sort -u ${f}-gilist-*.msk > ${gilist_all}

	# Note that this search (from original on disk) is limited
	# to one thread irrespective of what is used above. This
	# can prevent contention and is usually significantly quicker
	# than using all the threads available.

	export BLASTDB=${MY_BLASTDB}

	${MY_BLAST} -num_threads 1 -query ${f} \
	    -dbsize ${my_dbsize} -gilist ${gilist_all} -out ${stub}.blast

	ENDTIME=$(date +%s)
	echo "-> Time to merge $f   $(($ENDTIME-$STARTTIME)) s"
    fi

    # Clean up

    mpi_barrier

    rm -f ${gilist}
    rm -f ${gilist_all}

    return 0
}

# Execute and finish

main "$@"

