#!/bin/bash
# Para ejecutarlo en la terminal:
# El 'caffeinate' es para que la terminal de la mac no pause el proceso cuando entra en reposo
# caffeinate ./run_antismash_clados.sh
#!/bin/bash

# Parar ejecución en caso de error crítico
set -e

# Número de hilos CPU a utilizar por genoma
THREADS=4

echo "Iniciando análisis masivo con antiSMASH..."

# Recorrer la carpeta del clado
for dir in cladoA_completo_genomas; do
    if [ -d "$dir" ]; then
        clado_name="clado_A"
        echo ""
        echo "--------------------------------------------------"
        echo ">>> CLADO ACTUAL: $clado_name"
        echo "--------------------------------------------------"

        for fasta in "$dir"/*.fna; do
            [ -e "$fasta" ] || continue

            sample=$(basename "$fasta" .fna)
            output_path="resultados_antismash/${clado_name}/${sample}"

            # Crear carpeta de salida completa
            mkdir -p "$output_path"

            if [ -f "${output_path}/index.html" ]; then
                echo "[SALTAR] $sample ya fue procesado correctamente."
                continue
            fi

            echo "[EJECUTANDO] Procesando $sample ..."

            antismash \
                --clusterhmmer \
                --tigrfam \
                --asf \
                --cc-mibig \
                --cb-general \
                --cb-subclusters \
                --cb-knownclusters \
                --pfam2go \
                --tfbs \
                --genefinding-tool prodigal \
                --html-title "${clado_name}_${sample}" \
                --output-dir "$output_path" \
                -c "$THREADS" \
                "$fasta"

            echo "[OK] $sample finalizado."
        done
    fi
done

echo ""
echo "=================================================="
echo "¡Análisis finalizado para todos los clados!"
echo "=================================================="





