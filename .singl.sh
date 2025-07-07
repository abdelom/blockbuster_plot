#!/bin/bash

# Dossier de sortie
output_dir="resultats_blks"
mkdir -p "$output_dir"

# Fichier contenant les données à tracer
plot_data="plot_data.txt"
> "$plot_data"  # Vide le fichier s'il existe

# Boucle sur les tailles d'échantillon de 2 à 1000 par pas de 20
for n in $(seq 2 20 1000); do
    outfile="$output_dir/sfs_n${n}.blk"
    ./bin/blockbuster_simulator -n "$n" -t 1e9 -f 1 -e 1 -p 1e-4,0.001 -o "$outfile"

    # Extraire les deux premiers nombres de la première ligne
    if [ -s "$outfile" ]; then
        read v1 v2 <<< $(head -n 1 "$outfile")
        echo "$n $v1 $v2" >> "$plot_data"
    fi
done

# Tracer avec gnuplot
gnuplot -persist <<EOF
set title "Valeurs extraites vs Taille d'échantillon"
set xlabel "Taille d'échantillon (n)"
set ylabel "Valeurs (1er et 2e nombres)"
set grid
plot \
    "$plot_data" using 1:2 with linespoints title "1er nombre", \
    "$plot_data" using 1:3 with linespoints title "2e nombre"
EOF
