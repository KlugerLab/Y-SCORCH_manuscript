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
        ATAC_WHITE)
            ATAC_WHITE="$value"
            ;;
        *)
            echo "Warning: Unknown parameter '$key' in the input file."
            ;;
    esac
done < "$INPUT_FILE"

# Create the output folder if it doesn't exist
if [[ ! -d "$OUTPUT_FOLDER" ]]; then
    mkdir -p "$OUTPUT_FOLDER"
fi

# Create an array to store file names
declare -a FILE_NAMES

# Find files matching the naming convention and store them in array
while IFS= read -r -d $'\0' file; do
    FILE_NAME=$(basename "$file")

    # Extract the name part
    NAME_PART=$(echo "$FILE_NAME" | awk -F '_S[0-9]*' '{print $1}')

    # Check if name part contains "ATAC"
    if [[ $NAME_PART =~ "ATAC" ]]; then
        FILE_NAMES+=("$NAME_PART|$FILE_NAME")
    fi
done < <(find -L $FASTQ_FOLDER -type f -name '*_S*_L*_R3_*.fastq.gz' -print0)

process_file() {
    file_info=$1
    IFS="|" read -r NAME_PART FILE_NAME <<< "$file_info"
    if [[ ! -d "$OUTPUT_FOLDER/$NAME_PART" ]]; then
        mkdir $OUTPUT_FOLDER/$NAME_PART
    fi

    # Input files
    files=("$OUTPUT_FOLDER/${NAME_PART}/match_R3_viral_unique.txt" "$OUTPUT_FOLDER/${NAME_PART}/match_R2_viral_unique.txt" "$OUTPUT_FOLDER/${NAME_PART}/match_R1_viral_unique.txt")

    # Output files
    output_files=("$OUTPUT_FOLDER/${NAME_PART}/filtered_match_R3_viral_unique.txt" "$OUTPUT_FOLDER/${NAME_PART}/filtered_match_R2_viral_unique.txt" "$OUTPUT_FOLDER/${NAME_PART}/filtered_match_R1_viral_unique.txt")

    # Temporary arrays to hold triplet lines
    declare -a lines1 lines2 lines3

    # Set to store unique triplets
    declare -A seen_triplets

    # Ensure output files are empty
    for output_file in "${output_files[@]}"; do
        > "$output_file"
    done

    # Open all files in parallel and process them
    exec 3<"${files[0]}" 4<"${files[1]}" 5<"${files[2]}"

    while true; do
        # Read 4-line blocks from each file
        read -r line1_1 <&3 || break
        read -r line2_1 <&3
        read -r line3_1 <&3
        read -r line4_1 <&3

        read -r line1_2 <&4 || break
        read -r line2_2 <&4
        read -r line3_2 <&4
        read -r line4_2 <&4

        read -r line1_3 <&5 || break
        read -r line2_3 <&5
        read -r line3_3 <&5
        read -r line4_3 <&5

	# Create two partial triplets
	partial_triplet1="$line2_1|${line2_2: -16}"
	partial_triplet2="${line2_2: -16}|$line2_3"

	# Check if either partial triplet has been seen before
	if [[ -z "${seen_triplets[$partial_triplet1]}" && -z "${seen_triplets[$partial_triplet2]}" && ! $partial_triplet1 =~ (G{15}|C{15}|T{15}|A{15}) && ! $partial_triplet2 =~ (G{15}|C{15}|T{15}|A{15}) ]]; then
	    # Mark both partial triplets as seen
	    seen_triplets["$partial_triplet1"]=1
	    seen_triplets["$partial_triplet2"]=1

        # Combine the second lines into a triplet
#	triplet="$line2_1|${line2_2: -16}|$line2_3"

        # Check if the triplet has been seen before
#        if [[ -z "${seen_triplets[$triplet]}" && ! $triplet =~ (G{15}|C{15}|T{15}|A{15}) ]]; then
            # Mark the triplet as seen
#            seen_triplets["$triplet"]=1

            # Append the blocks to the output files
            {
                echo "$line1_1"
                echo "$line2_1"
                echo "$line3_1"
                echo "$line4_1"
            } >> "${output_files[0]}"

            {
                echo "$line1_2"
                echo "$line2_2"
                echo "$line3_2"
                echo "$line4_2"
            } >> "${output_files[1]}"

            {
                echo "$line1_3"
                echo "$line2_3"
                echo "$line3_3"
                echo "$line4_3"
            } >> "${output_files[2]}"
        fi
    done

    # Close file descriptors
    exec 3<&- 4<&- 5<&-

    cat "${output_files[1]}" | \
    # Main awk script
    awk -v atac_white="$ATAC_WHITE" -v gex_white="$GEX_WHITE" '
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

    function reverse_complement(seq) {
        # Calculate the reverse complement of a DNA sequence
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

    BEGIN {
        # Read ATAC_WHITE into an array
	atac_counter = 0
        while ((getline line < atac_white) > 0) {
            atac_array[++atac_counter] = line
        }
	close(atac_white)

	gex_counter = 0
        # Read GEX_WHITE into an array
        while ((getline line < gex_white) > 0) {
            gex_map[atac_array[++gex_counter]] = line
        }
	close(gex_white)
    }

    NR % 4 == 2 {
        # Extract the key from the sequence
        key = substr($0, length($0) - 15, 16)
        seen[key]++
    }

    END {
        # Sort keys in descending order of counts
        m = asorti(seen, sorted_keys, "@val_num_desc")

        for (j = 1; j <= m; j++) {
            key = sorted_keys[j]
            count = seen[key]
            rc_key = reverse_complement(key)

            # Check for an exact match in preloaded atac_map
            if (rc_key in gex_map) {
                final_map[rc_key] += count
                continue
            }

	    approx_found = 0
	    match_key = ""

            # Check previously processed keys in final_map for approximate matches
            for (prev_key in final_map) {
                if (final_map[prev_key] > 1 && hamming_distance(rc_key, prev_key) == 1) {
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
            for (atac_key in gex_map) {
                if (hamming_distance(atac_key, rc_key) == 1) {
                    approx_found++
                    if (approx_found > 1) break
		    match_key = atac_key
                }
            }

            if (approx_found == 1) {
                final_map[match_key] += count
            }
        }

        # Print results
        for (key in final_map) {
            print gex_map[key] ", " final_map[key]
        }
    }' | \
    sort -t, -k2,2nr > $OUTPUT_FOLDER/${NAME_PART}/filtered_barcode_viral_unique.txt
}

export -f process_file
export OUTPUT_FOLDER FASTQ_FOLDER ATAC_WHITE GEX_WHITE

# Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
parallel -j 4 process_file ::: "${FILE_NAMES[@]}"
