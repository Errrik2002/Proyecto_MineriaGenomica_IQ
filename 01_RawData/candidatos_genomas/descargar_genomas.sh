#!/bin/bash

# metadata
metadata="candidatos.csv"

Red="$(tput setaf 1)"
Green="$(tput setaf 2)"
Bold=$(tput bold)
reset=$(tput sgr0)

while IFS=, read -r field1; do
    # Omitir líneas vacías o espacios en blanco
    [[ -z "$field1" ]] && continue

    # Limpiar posibles retornos de carro (\r) si el CSV se hizo en Windows
    field1=$(echo "$field1" | tr -d '\r' | xargs)

    echo ""
    echo "${Red}${Bold}Downloading ${reset}: ${field1}"
    
    # 1. Limpiar restos de ejecuciones anteriores para evitar "ncbi_dataset (1)"
    rm -rf "${field1}.zip" ncbi_dataset* README.md md5sum.txt
    
    # 2. Descargar
    datasets download genome accession "${field1}" --filename "${field1}.zip"
    
    if [ ! -f "${field1}.zip" ]; then
        echo "${Red}Error: No se pudo descargar ${field1}. Pasando al siguiente...${reset}"
        continue
    fi

    echo "${Bold}Extracting ${field1}.zip ${reset}"
    unzip -q "${field1}.zip"
    
    # 3. Entrar a la carpeta de manera segura con subshell (paréntesis)
    (
        cd "ncbi_dataset/data/${field1}" || exit 1
        echo "${Bold}Moving ${field1} fasta file into home directory${reset}"
        mv *.fna ../../../
    )

    # 4. Limpieza final de la iteración
    rm -rf "${field1}.zip" ncbi_dataset* README.md md5sum.txt
    
    echo "${Green}${Bold}Download_completed ${reset}: ${field1}"
    echo ""

done < "${metadata}"
