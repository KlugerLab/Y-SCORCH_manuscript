#!/bin/bash

# Ask user for number of cores
read -p "Enter the number of available cores: " CORES

# Ask user for viral index location
read -p "Enter the viral index location: " VIRAL_INDEX
VIRAL_INDEX_FOLDER=$(dirname $VIRAL_INDEX)

# Ask user for human index location
read -p "Enter the human index location: " HUMAN_INDEX
HUMAN_INDEX_FOLDER=$(dirname $HUMAN_INDEX)

# Check if viral index folder exists
if [[ ! -d "$VIRAL_INDEX_FOLDER" ]]; then
    mkdir $VIRAL_INDEX_FOLDER
    read -p "Enter the viral fasta file: " VIRAL_FASTA
    bowtie2-build $VIRAL_FASTA $VIRAL_INDEX
fi

# Check if viral index folder exists
if [[ ! -d "$HUMAN_INDEX_FOLDER" ]]; then
    mkdir $HUMAN_INDEX_FOLDER
    read -p "Enter the human fasta file: " HUMAN_FASTA
    bowtie2-build $HUMAN_FASTA $HUMAN_INDEX
fi

# Ask for output folder
read -p "Enter the output folder for bowtie2 run: " OUTPUT_FOLDER
if [[ ! -d "$OUTPUT_FOLDER" ]]; then
    mkdir $OUTPUT_FOLDER
fi

# Ask for fastq files location
read -p "Enter the location of the fastq files: " FASTQ_FOLDER

# Ask for gex whitelist location
read -p "Enter the location of the gex whitelist file: " GEX_WHITE

# Ask for atac whitelist location
read -p "Enter the location of the atac whitelist file: " ATAC_WHITE

# Create an array to store file names
declare -a FILE_NAMES

# Find files matching the naming convention and store them in array
while IFS= read -r -d $'\0' file; do
    FILE_NAME=$(basename "$file")

    # Extract the name part
    NAME_PART=$(echo "$FILE_NAME" | awk -F '_S[0-9]*' '{print $1}')

    # Check if name part contains "ATAC"
    if [[ $NAME_PART =~ "ATAC" ]]; then
#    if [[ $NAME_PART =~ "A-" ]]; then
        FILE_NAMES+=("$NAME_PART|$FILE_NAME")
    fi
done < <(find -L $FASTQ_FOLDER -type f -name '*_S*_L*_R3_*.fastq.gz' -print0)
#done < <(find $FASTQ_FOLDER -type l -name '*_S*_L*_R3_*.fastq.gz' -print0)

# Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
max_samples=$((CORES / 8))
max_samples=$((max_samples < ${#FILE_NAMES[@]} ? max_samples : ${#FILE_NAMES[@]}))

# Ask user for amount of memory
read -p "Enter the amount of available memory in GB: " MAX_MEMORY
max_runs_mem=$((MAX_MEMORY / 100))
max_runs=$((max_runs_mem < max_samples ? max_runs_mem : max_samples))

# Determine number of cores and memory per run
cores_per_sample=$((CORES / max_samples))

running_jobs=0

# Run bowtie2 for each file
for file_info in "${FILE_NAMES[@]}"; do
    IFS="|" read -r NAME_PART FILE_NAME <<< "$file_info"
    if [[ ! -d "$OUTPUT_FOLDER/$NAME_PART" ]]; then
        mkdir $OUTPUT_FOLDER/$NAME_PART
    fi
#    if [ ! -e "$OUTPUT_FOLDER/$NAME_PART/output.bam" ]; then
	# Find related FASTQ file
        RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R1_*.fastq.gz" -print -quit)
#	bowtie2 --no-mixed --no-discordant --local -x $VIRAL_INDEX -1 $RELATED_FILE -2 $FASTQ_FOLDER/$FILE_NAME -N 1 -L 10 --very-sensitive -p $cores_per_sample --no-unal | samtools view -bS -o "$OUTPUT_FOLDER/$NAME_PART/output.bam" &
	bowtie2 --no-discordant --local -x $VIRAL_INDEX -1 $RELATED_FILE -2 $FASTQ_FOLDER/$FILE_NAME -N 1 -L 10 --very-sensitive -p $cores_per_sample --no-unal | samtools view -bS -o "$OUTPUT_FOLDER/$NAME_PART/output.bam" &
#    else
#	continue
#    fi
    ((running_jobs++))

    # Check to see if as many jobs as possible are running
    if [ "$running_jobs" -eq "$max_samples" ]; then
        wait -n
        ((running_jobs--))
    fi
done

wait

# Function to process each file
process_file() {
    file_info=$1
    IFS="|" read -r NAME_PART _ <<< "$file_info"

    # Find related FASTQ file
    RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R2_*.fastq.gz" -print -quit)

    # Extract unique read names, sequences, count the occurrences of the first 16 characters, and format the output
    samtools view -h $OUTPUT_FOLDER/${NAME_PART}/output.bam | \
    awk '$0 !~ /^@/ && $0 !~ /GGGGGGGGGGGGGGG/ && $0 !~ /TTTTTTTTTTTTTTT/ && $0 !~ /AAAAAAAAAAAAAAA/ && $0 !~ /CCCCCCCCCCCCCCC/ {if (!seen[$1]++) print $1}' > $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt
#    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt | \
    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R2_viral.txt
    RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R3_*.fastq.gz" -print -quit)
    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R3_viral.txt
    RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R1_*.fastq.gz" -print -quit)
    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R1_viral.txt
    cat $OUTPUT_FOLDER/${NAME_PART}/match_R2_viral.txt | \

    # Main awk script
    awk -v atac_file="$ATAC_WHITE" -v gex_file="$GEX_WHITE" '
    function reverse_complement(seq) {
        n = length(seq)
        rc = ""
        for (i = n; i > 0; i--) {
            base = substr(seq, i, 1)
            if (base == "A") rc = rc "T"
            else if (base == "T") rc = rc "A"
            else if (base == "C") rc = rc "G"
            else if (base == "G") rc = rc "C"
        }
        return rc
    }

    NR % 4 == 2 {
        key = substr($0, length($0) - 15, 16)
#        rest = substr($0, 1, length($0) - 16)
#	combined_key = key "|" rest  # Create a compound key
#        seen[combined_key] = 1
	seen[key]++
    }

    END {
#        for (combined_key in seen) {
#            split(combined_key, parts, "|")
#            key = parts[1]
#	    unique_count[key]++
#	}

#	for (key in unique_count) {
#	    count = unique_count[key]
	for (key in seen) {
	    count = seen[key]

            # Calculate reverse complement of the key
            rc_key = reverse_complement(key)

            # Find the line number in ATAC_WHITE file
            while ((getline line < atac_file) > 0 && (getline gex_line < gex_file) > 0) {
                if (line ~ rc_key) {
                    # Get the corresponding line from GEX_WHITE file
                    print gex_line ", " count
                    break
                }
            }

            # Close the files for the next iteration
            close(atac_file)
            close(gex_file)
        }
    }' | \
    sort -t, -k2,2nr > $OUTPUT_FOLDER/${NAME_PART}/barcode_counts.txt
}

export -f process_file
export OUTPUT_FOLDER FASTQ_FOLDER ATAC_WHITE GEX_WHITE

# Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
parallel -j $max_runs process_file ::: "${FILE_NAMES[@]}"

# Function to process each file
get_seqs() {
    file_info=$1
    IFS="|" read -r NAME_PART FILE_NAME <<< "$file_info"

    RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R1_*.fastq.gz" -print -quit)
    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt | gzip > $OUTPUT_FOLDER/${NAME_PART}/match_R1.fastq.gz
    seqtk subseq $FASTQ_FOLDER/$FILE_NAME $OUTPUT_FOLDER/${NAME_PART}/match_labels.txt | gzip > $OUTPUT_FOLDER/${NAME_PART}/match_R3.fastq.gz
}

export -f get_seqs
export OUTPUT_FOLDER FASTQ_FOLDER

# Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
parallel -j $max_runs get_seqs ::: "${FILE_NAMES[@]}"

running_jobs=0

# Run bowtie2 for each file
for file_info in "${FILE_NAMES[@]}"; do
    IFS="|" read -r NAME_PART FILE_NAME <<< "$file_info"
    if [[ ! -d "$OUTPUT_FOLDER/$NAME_PART" ]]; then
        mkdir $OUTPUT_FOLDER/$NAME_PART
    fi
#    if [ ! -e "$OUTPUT_FOLDER/$NAME_PART/human_output.bam" ]; then
        bowtie2 --local -x $HUMAN_INDEX -1 $OUTPUT_FOLDER/${NAME_PART}/match_R1.fastq.gz -2 $OUTPUT_FOLDER/${NAME_PART}/match_R3.fastq.gz -N 1 -L 10 --very-sensitive -p $cores_per_sample --no-unal | samtools view -bS -o "$OUTPUT_FOLDER/$NAME_PART/human_output.bam" &
#    else
#       continue
#    fi
    ((running_jobs++))

    # Check to see if as many jobs as possible are running
    if [ "$running_jobs" -eq "$max_samples" ]; then
        wait -n
        ((running_jobs--))
    fi
done

wait

# Function to process each file
process_file_human() {
    file_info=$1
    IFS="|" read -r NAME_PART _ <<< "$file_info"

    # Find related FASTQ file
    RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R2_*.fastq.gz" -print -quit)

    # Extract unique read names, sequences, count the occurrences of the first 16 characters, and format the output
    samtools view -h $OUTPUT_FOLDER/${NAME_PART}/human_output.bam | \
    awk '$0 !~ /^@/ && $0 !~ /GGGGGGGGGGGGGGG/ && $0 !~ /AAAAAAAAAAAAAAA/ && $0 !~ /CCCCCCCCCCCCCCC/ && $0 !~ /TTTTTTTTTTTTTTT/ {if (!seen[$1]++) print $1}' | \
#    seqtk subseq $RELATED_FILE - | \
    seqtk subseq $RELATED_FILE - > $OUTPUT_FOLDER/${NAME_PART}/match_R2.txt
    cat $OUTPUT_FOLDER/${NAME_PART}/match_R2.txt | \

    # Main awk script
    awk -v atac_file="$ATAC_WHITE" -v gex_file="$GEX_WHITE" '
    function reverse_complement(seq) {
        n = length(seq)
        rc = ""
        for (i = n; i > 0; i--) {
            base = substr(seq, i, 1)
            if (base == "A") rc = rc "T"
            else if (base == "T") rc = rc "A"
            else if (base == "C") rc = rc "G"
            else if (base == "G") rc = rc "C"
        }
        return rc
    }

    NR % 4 == 2 {
        key = substr($0, length($0) - 15, 16)
        seen[key]++
    }

    END {
        for (key in seen) {
            count = seen[key]

            # Calculate reverse complement of the key
            rc_key = reverse_complement(key)

            # Find the line number in ATAC_WHITE file
            while ((getline line < atac_file) > 0 && (getline gex_line < gex_file) > 0) {
                if (line ~ rc_key) {
                    # Get the corresponding line from GEX_WHITE file
                    print gex_line ", " count
                    break
                }
            }

            # Close the files for the next iteration
            close(atac_file)
            close(gex_file)
        }
    }' | \
    sort -t, -k2,2nr > $OUTPUT_FOLDER/${NAME_PART}/barcode_human.txt
}

export -f process_file_human
export OUTPUT_FOLDER FASTQ_FOLDER ATAC_WHITE GEX_WHITE

# Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
parallel -j $max_runs process_file_human ::: "${FILE_NAMES[@]}"

# Function to process each file
process_file_unique() {
    file_info=$1
    IFS="|" read -r NAME_PART _ <<< "$file_info"

    # Find unique sets of 4 lines
    awk 'NR==FNR {
        block = block $0;
        if (FNR % 4 == 0) {
            a[block];
            block = ""
        } else {
            block = block "\n"
        }
        next
    }
    {
        block = block $0;
        if (FNR % 4 == 0) {
            if (!(block in a)) {
                print block
            }
            block = ""
        } else {
            block = block "\n"
        }
    }' $OUTPUT_FOLDER/${NAME_PART}/match_R2.txt $OUTPUT_FOLDER/${NAME_PART}/match_R2_viral.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R2_viral_unique.txt

    # Process the unique sets
    awk 'NR%4==1 {print substr($0, 2, index($0, " ") - 2)}' $OUTPUT_FOLDER/${NAME_PART}/match_R2_viral_unique.txt > $OUTPUT_FOLDER/${NAME_PART}/unique_labels.txt
    RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R1_*.fastq.gz" -print -quit)
    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/unique_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R1_viral_unique.txt
    RELATED_FILE=$(find -L $FASTQ_FOLDER -type f -name "${NAME_PART}_S*_L*_R3_*.fastq.gz" -print -quit)
    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}/unique_labels.txt > $OUTPUT_FOLDER/${NAME_PART}/match_R3_viral_unique.txt
    cat $OUTPUT_FOLDER/${NAME_PART}/match_R2_viral_unique.txt | \
    # Main awk script
    awk -v atac_file="$ATAC_WHITE" -v gex_file="$GEX_WHITE" '
    function reverse_complement(seq) {
        n = length(seq)
        rc = ""
        for (i = n; i > 0; i--) {
            base = substr(seq, i, 1)
            if (base == "A") rc = rc "T"
            else if (base == "T") rc = rc "A"
            else if (base == "C") rc = rc "G"
            else if (base == "G") rc = rc "C"
        }
        return rc
    }

    NR % 4 == 2 {
        key = substr($0, length($0) - 15, 16)
        seen[key]++
    }

    END {
        for (key in seen) {
            count = seen[key]

            # Calculate reverse complement of the key
            rc_key = reverse_complement(key)

            # Find the line number in ATAC_WHITE file
            while ((getline line < atac_file) > 0 && (getline gex_line < gex_file) > 0) {
                if (line ~ rc_key) {
                    # Get the corresponding line from GEX_WHITE file
                    print gex_line ", " count
                    break
                }
            }

            # Close the files for the next iteration
            close(atac_file)
            close(gex_file)
        }
    }' | \
    sort -t, -k2,2nr > $OUTPUT_FOLDER/${NAME_PART}/barcode_viral_unique.txt
}

export -f process_file_unique
export OUTPUT_FOLDER FASTQ_FOLDER ATAC_WHITE GEX_WHITE

# Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
parallel -j 4 process_file_unique ::: "${FILE_NAMES[@]}"

chmod 777 $OUTPUT_FOLDER
chmod 777 $OUTPUT_FOLDER/*
