#!/bin/bash

if [ $# -lt 2 ] ; then
    echo ""
    echo "usage: count_fastq.sh [output_file] [fastq_file1] <fastq_file2> ..|| *.fastq.gz"
    echo "counts the number of reads in a fastq file"
    echo ""
    exit 0
fi

output_file="$1"
shift

# Initialiser le premier nom de fichier
previous_file=""

# Imprimer l'en-tête initial
printf "%-50s %-10s\n" "File_Name" "Reads_Count" > "$output_file"
printf "%-50s %-10s\n" "----------------" "-----------" >> "$output_file"

for i in "$@"; do
    # Extraire le nom de fichier sans extension
    filename=$(basename "$i" | cut -d. -f1)

    # Si le nom de fichier actuel est différent du nom de fichier précédent, imprimer les en-têtes pour le nouveau groupe de fichiers
    if [ "$filename" != "$previous_file" ]; then
        if [ "$previous_file" != "" ]; then
            printf "%-50s %-10s\n" "----------------" "-----------" >> "$output_file"
        fi
        previous_file="$filename"
    fi

    if [[ "$i" == *.gz ]]; then
        lines=$(zcat "$i" | wc -l)
    else
        lines=$(wc -l < "$i")
    fi

    count=$((lines / 4))

    printf "%-50s %-10d\n" "$(basename "$i")" "$count" >> "$output_file"  # Écrire le nom du fichier et le nombre de lectures dans le fichier de sortie
done

# Imprimer une ligne de séparation à la fin
printf "%-50s %-10s\n" "----------------" "-----------" >> "$output_file"
