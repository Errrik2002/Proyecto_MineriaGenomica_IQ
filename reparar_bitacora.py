archivo = "bitacora.ipynb"  # Sustituye por el nombre real de tu archivo

with open(archivo, "r", encoding="utf-8") as f:
    lineas = f.readlines()

# Filtra las líneas con marcas de conflicto de Git
lineas_limpias = [
    linea for linea in lineas 
    if not (linea.startswith("=======") or linea.startswith(">>>>>>>") or linea.startswith("<<<<<<<"))
]

with open("archivo_reparado.ipynb", "w", encoding="utf-8") as f:
    f.writelines(lineas_limpias)

print("¡Listo! Intenta abrir 'archivo_reparado.ipynb'.")