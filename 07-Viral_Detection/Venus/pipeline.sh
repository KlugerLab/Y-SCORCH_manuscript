#!/bin/bash

source activate venus

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
        VENUS_INDEX)
            VENUS_INDEX="$value"
            ;;
        HUMAN_FASTA)
            HUMAN_FASTA="$value"
            ;;
        HUMAN_GTF)
            HUMAN_GTF="$value"
            ;;
        CORES)
            CORES="$value"
            ;;
        VIRAL_FASTA)
            VIRAL_FASTA="$value"
            ;;
        VIRAL_GTF)
            VIRAL_GTF="$value"
            ;;
        VENUS_OUTPUT)
            VENUS_OUTPUT="$value"
            ;;
        FASTQ)
            FASTQ="$value"
            ;;
        VIRUS_THRESHOLD)
            VIRUS_THRESHOLD="$value"
            ;;
        BARCODE_LENGTH)
            BARCODE_LENGTH="$value"
            ;;
        UMI_LENGTH)
            UMI_LENGTH="$value"
            ;;
        WHITELIST)
            WHITELIST="$value"
            ;;
        CHR_REF)
            CHR_REF="$value"
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

process_file() {
    file_info=$1
    IFS="|" read -r NAME_PART FILE_NAME <<< "$file_info"
    # Use awk to filter and count unique fields
    awk '
        BEGIN {
            OFS = "\t";
        }
        !/^@/ && NF > 8 {
            if ($10 !~ /(A{15}|C{15}|G{15}|T{15})/) {
                # Find CR:Z and UR:Z fields
                for (i = 1; i <= NF; i++) {
                    if ($i ~ /^CR:Z:/) {
                        barcode = substr($i, 6);
                    }
                    if ($i ~ /^UR:Z:/) {
                        umi = substr($i, 6);
                    }
                }
                if (barcode && umi) {
                    key = barcode OFS umi;
                    if (!seen[key]++) {
                        unique_first_fields[barcode]++;
                    }
                }
            }
        }
        END {
            for (barcode in unique_first_fields) {
                print barcode ", " unique_first_fields[barcode];
            }
        }
    ' "$VENUS_OUTPUT/detection/$NAME_PART/virus/Aligned.out.sam" | sort -t, -k2,2nr > "$VENUS_OUTPUT/detection/$NAME_PART/barcode_count.txt"
}

if [ ! -d "$VENUS_INDEX/indices" ]; then
    python3 /opt/Venus/src/module-index/module-index.py \
        --humanFASTA "$HUMAN_FASTA" \
        --humanGTF "$HUMAN_GTF" \
        --virusFASTA "$VIRAL_FASTA" \
        --module detection \
        --thread "$CORES" \
        --out "$VENUS_INDEX/indices" \
        --hGenome "$VENUS_INDEX/indices/human.genomeDir" \
        --virusGenome "$VENUS_INDEX/indices/mega_virus.genomeDir"

#    python3 /opt/Venus/src/module-index/module-index.py \
#        --humanFASTA "$HUMAN_FASTA" \
#        --humanGTF "$HUMAN_GTF" \
#        --virusFASTA "$VIRAL_FASTA" \
#        --virusGTF "$VIRAL_GTF" \
#        --module integration \
#        --thread "$CORES" \
#        --out "$VENUS_INDEX/indices" \
#        --hGenome "$VENUS_INDEX/indices/hybrid.genomeDir" \
#        --virusGenome "$VENUS_INDEX/indices/virus.genomeDir"
fi

# Find files matching the naming convention and store them in array
# Create an array to store file names
declare -a FILE_NAMES

while IFS= read -r -d $'\0' file; do
    FILE_NAME=$(basename "$file")

    # Extract the name part
    NAME_PART=$(echo "$FILE_NAME" | awk -F '_S[0-9]*' '{print $1}')

    # Check if name part does not contain "ATAC"
    if [[ ! $NAME_PART =~ "ATAC" ]]; then
        FILE_NAMES+=("$NAME_PART|$FILE_NAME")
    fi
done < <(find -L $FASTQ -type f -name "$READ_GREP" -print0)

# Calculate the largest possible number of parallel runs given a minimum of 8 cores per run
max_samples_cores=$((CORES / 8))
max_samples_cores=$((max_samples_cores < ${#FILE_NAMES[@]} ? max_samples_cores : ${#FILE_NAMES[@]}))

# Determine number of cores and memory per run
cores_per_sample=$((CORES / max_samples_cores))

running_jobs=0

if [ ! -d "$VENUS_OUTPUT" ]; then
    # Loop through samples
    for index in "${FILE_NAMES[@]}"; do
	IFS="|" read -r NAME_PART FILE_NAME <<< "$index"
	cDNA_FILE=$(find -L $FASTQ -type f -name "${NAME_PART}${READ_GREP}" -print -quit)
	CB_UMI_FILE=$(find -L $FASTQ -type f -name "${NAME_PART}${CB_GREP}" -print -quit)
        python3 /opt/Venus/src/module-detection/module-detection.py \
            --read "$cDNA_FILE" "$CB_UMI_FILE" \
            --virusThreshold "$VIRUS_THRESHOLD" \
            --virusChrRef "$CHR_REF" \
            --virusGenome "$VENUS_INDEX/indices/mega_virus.genomeDir" \
            --humanGenome "$VENUS_INDEX/indices/human.genomeDir" \
            --readFilesCommand zcat \
            --thread "$cores_per_sample" \
            --singleCellBarcode 1 "$BARCODE_LENGTH" \
            --singleUniqueMolIdent "$((BARCODE_LENGTH + 1))" "$UMI_LENGTH" \
            --singleWhitelist "$WHITELIST" \
            --out "$VENUS_OUTPUT/detection/${NAME_PART}" &
        ((running_jobs++))

        # Check to see if as many jobs as possible are running
        if [ "$running_jobs" -eq "$max_samples_cores" ]; then
            wait -n
	    ((running_jobs--))
        fi
    done

    # Wait for any remaining jobs to finish
    wait
fi

all_exist=true
for dir in "$VENUS_OUTPUT/detection/"*/; do
    if [ ! -f "${dir}output.txt" ]; then
        all_exist=false
        break
    fi
done

if ! $all_exist; then
    export VENUS_OUTPUT
    export -f process_file

    # Use parallel to process files in parallel, limiting to MAX_PARALLEL_JOBS
    parallel -j "$max_samples_cores" process_file ::: "${FILE_NAMES[@]}"

fi

    #running_jobs=0

    #read -p "Enter the assembly name of the viruses: " virus_name

    # Loop through samples
    #for index in "${!cDNAs[@]}"; do
    #    python3 /opt/Venus/src/module-integration/module-integration.py \
    #        --read "${cDNAs[$index]}" \
    #        --virusGenome "$output_location/indices/mega_virus.genomeDir" \
    #        --hybridGenome "$output_location/indices/hybrid.genomeDir" \
    #        --guideFASTA "$output_location/reference_files/integrSeq.fna" \
    #        --readFilesCommand zcat \
    #        --virusChr "$virus_name" \
    #        --thread "$cores_per_sample" \
    #        --geneBed /opt/Venus/reference_files/genes.bed \
    #        --out "$output_location/integration/single-end" &
    #    ((running_jobs++))

        # Check to see if as many jobs as possible are running
    #    if [ "$running_jobs" -eq "$max_concurrent_jobs" ]; then
    #        wait -n
    #        ((running_jobs--))
    #    fi
    #done

    # Wait for any remaining jobs to finish
    #wait

conda deactivate
chmod 777 $VENUS_OUTPUT/detection/*/virus
chmod 777 $VENUS_OUTPUT/detection/*/virus/*
