#!/bin/bash

# Check if the input file is provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

# Read the input file
INPUT_FILE="$1"

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: File '$INPUT_FILE' not found."
    exit 1
fi

# Parse the input file
while IFS="=" read -r key value; do
    case "$key" in
        OUTPUT_FOLDER)
            OUTPUT_FOLDER="$value"
            ;;
        FASTQ_FOLDER)
            FASTQ_FOLDER="$value"
            ;;
        GEX_WHITE)
            GEX_WHITE="$value"
            ;;
        CORES)
            CORES="$value"
            ;;
        VIRAL_INDEX_FOLDER)
            VIRAL_INDEX_FOLDER="$value"
            ;;
        VIRAL_FASTA)
            VIRAL_FASTA="$value"
            ;;
        ADJUSTMENT_FACTOR)
            ADJUSTMENT_FACTOR="$value"
            ;;
	READ_GREP)
	    READ_GREP="$value"
	    ;;
	CB_GREP)
	    CB_GREP="$value"
	    ;;
        *)
            echo "Warning: Unknown parameter '$key' in the input file."
            ;;
    esac
done < "$INPUT_FILE"

# Check if viral index folder exists
if [[ ! -d "$VIRAL_INDEX_FOLDER" ]]; then
    STAR --runThreadN $CORES --runMode genomeGenerate --genomeDir $VIRAL_INDEX_FOLDER --genomeFastaFiles $VIRAL_FASTA --genomeSAindexNbases $ADJUSTMENT_FACTOR
fi

# Create an array to store file names
declare -a FILE_NAMES

# Find files matching the naming convention and store them in array
while IFS= read -r -d $'\0' file; do
    FILE_NAME=$(basename "$file")

    # Extract the name part
    NAME_PART=$(echo "$FILE_NAME" | awk -F '_S[0-9]*' '{print $1}')

    # Check if name part does not contain "ATAC"
    if [[ ! $NAME_PART =~ "ATAC" ]]; then
        FILE_NAMES+=("$NAME_PART|$FILE_NAME")
    fi
done < <(find -L $FASTQ_FOLDER -type f -name "$READ_GREP" -print0)
#done < <(find $FASTQ_FOLDER -type l -name '*_S*_L*_R3_*.fastq.gz' -print0)

# Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
max_samples=$((CORES / 8))
max_samples=$((max_samples < ${#FILE_NAMES[@]} ? max_samples : ${#FILE_NAMES[@]}))

# Determine number of cores and memory per run
cores_per_sample=$((CORES / max_samples))

running_jobs=0

# Run STAR for each file
for file_info in "${FILE_NAMES[@]}"; do
    IFS="|" read -r NAME_PART FILE_NAME <<< "$file_info"
    if [ ! -e "$OUTPUT_FOLDER/${NAME_PART}Log.out" ]; then
        STAR --genomeDir $VIRAL_INDEX_FOLDER --runThreadN $cores_per_sample --readFilesIn $FASTQ_FOLDER/$FILE_NAME --readFilesCommand zcat --outFilterMultimapNmax 1000 --outFileNamePrefix $OUTPUT_FOLDER/$NAME_PART --genomeLoad LoadAndKeep --outSAMtype BAM Unsorted &
    else
	continue
    fi
    ((running_jobs++))

    # Check to see if as many jobs as possible are running
    if [ "$running_jobs" -eq "$max_samples" ]; then
        wait -n
        ((running_jobs--))
    fi
done

wait

# Run STAR again for each file with additional parameters
for file_info in "${FILE_NAMES[@]}"; do
    IFS="|" read -r NAME_PART FILE_NAME <<< "$file_info"
    if [ ! -e "$OUTPUT_FOLDER/${NAME_PART}_SJLog.out" ]; then
        STAR --genomeDir $VIRAL_INDEX_FOLDER --runThreadN $cores_per_sample --readFilesIn $FASTQ_FOLDER/$FILE_NAME --readFilesCommand zcat --sjdbFileChrStartEnd $OUTPUT_FOLDER/${NAME_PART}SJ.out.tab --outFilterMultimapNmax 1000 --outFileNamePrefix $OUTPUT_FOLDER/${NAME_PART}_SJ --genomeLoad NoSharedMemory --outSAMtype BAM Unsorted &
    else
       continue
    fi
    ((running_jobs++))

    # Check to see if as many jobs as possible are running
    if [ "$running_jobs" -eq "$max_samples" ]; then
        wait -n
        ((running_jobs--))
    fi
done

wait

# Process SAM files
for file_info in "${FILE_NAMES[@]}"; do
    IFS="|" read -r NAME_PART _ <<< "$file_info"
    samtools view $OUTPUT_FOLDER/${NAME_PART}_SJAligned.out.bam | grep -v "GGGGGGGGGGGGGGGGGGGG" > $OUTPUT_FOLDER/${NAME_PART}_cleanfile.sam
    awk 'NR>4{print $1}' $OUTPUT_FOLDER/${NAME_PART}_cleanfile.sam > $OUTPUT_FOLDER/${NAME_PART}_positivereads.txt
    sort $OUTPUT_FOLDER/${NAME_PART}_positivereads.txt | uniq > $OUTPUT_FOLDER/${NAME_PART}_uniquepos.txt
    mv $OUTPUT_FOLDER/${NAME_PART}_uniquepos.txt $OUTPUT_FOLDER/${NAME_PART}_positivereads.txt
    RELATED_FILE=$(find $FASTQ_FOLDER -type f -name "${NAME_PART}${CB_GREP}" -print -quit)
    seqtk subseq $RELATED_FILE $OUTPUT_FOLDER/${NAME_PART}_positivereads.txt > $OUTPUT_FOLDER/${NAME_PART}_positivereads.fq.gz

    cat "$OUTPUT_FOLDER/${NAME_PART}_positivereads.fq.gz" | \
    # Main awk script
    awk -v gex_white="$GEX_WHITE" '
    function hamming_distance(str1, str2) {
        # Calculate Hamming distance between two strings
        dist = 0
        for (i = 1; i <= length(str1); i++) {
            if (substr(str1, i, 1) != substr(str2, i, 1)) {
                dist++
                if (dist > 1) return dist # Early exit if distance > 1
            }
        }
        return dist
    }

    BEGIN {
        # Read GEX_WHITE into an array
        while ((getline line < gex_white) > 0) {
            gex_map[line] = 1
        }
        close(gex_white)
    }

    NR % 4 == 2 {
        # Check if the whole line has been seen before
        if (!line_seen[$0]++) {
            # Extract the key from the first 16 characters
            key = substr($0, 1, 16)
            seen[key]++
        }
    }

    END {
        # Sort keys in descending order of counts
        m = asorti(seen, sorted_keys, "@val_num_desc")

        for (j = 1; j <= m; j++) {
            key = sorted_keys[j]
            count = seen[key]

            # Check for an exact match in preloaded atac_map
            if (key in gex_map) {
                final_map[key] += count
                continue
            }

            approx_found = 0
            match_key = ""

            # Check previously processed keys in final_map for approximate matches
            for (prev_key in final_map) {
                if (final_map[prev_key] > 1 && hamming_distance(key, prev_key) == 1) {
                    approx_found++
                    if (approx_found > 1) break
                    match_key = prev_key
                }
            }

            if (approx_found > 1) continue

            if (approx_found == 1) {
                final_map[match_key] += count
                continue
            }

            # Check for approximate matches in atac_map
            for (gex_key in gex_map) {
                if (hamming_distance(gex_key, key) == 1) {
                    approx_found++
                    if (approx_found > 1) break
                    match_key = gex_key
                }
            }

            if (approx_found == 1) {
                final_map[match_key] += count
            }
        }

        # Print results
        for (key in final_map) {
            print key ", " final_map[key]
        }
    }' | \
    sort -t, -k2,2nr > $OUTPUT_FOLDER/${NAME_PART}_uniqbarcode_count.txt
done
chmod 777 $OUTPUT_FOLDER
chmod 777 $OUTPUT_FOLDER/*
