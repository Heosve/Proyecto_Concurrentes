# Análisis de rendimiento de Dotplots: Una visión en profundidad

## Descripción
Este proyecto realiza un análisis de dot plot entre dos secuencias de ADN utilizando Python.

## Requisitos

- Python 3.11
- Librerías requeridas: numpy, multiprocessing, time, Bio, mpi4py.

## Instalación

1. Clona el repositorio en tu máquina local.
2. Asegúrate de tener Python 3.11 instalado.
3. Instala las librerías necesarias ejecutando el siguiente comando:
pip install numpy multiprocessing time biopython

## Ejecucion secuencial

Para ejecutar el codigo secuencial solo es necesario la ejecucion del comando :

"python secuencial.py datos/E_coli.fna datos/Salmonella.fna 10000"

Dentro de la carpeta principal de los archivos, debe asegurarse se tener todas las dependencias instaladas en su kernel.

## Ejecucion Multiprocessing
Para ejecutar el código, utiliza el siguiente comando:
python archivo.py archivo1.fasta archivo2.fasta --maxlength

## Ejecucion MPI4PY

mpiexec -n 10 python .\archivo.py

Reemplaza `archivo.py` con el nombre del archivo que contiene el código.
Reemplaza `archivo1.fasta` y `archivo2.fasta` con las rutas de los archivos FASTA que contienen las secuencias de ADN que deseas comparar.
Opcionalmente, puedes especificar la opción `--maxlength` seguida de un número entero para limitar la longitud máxima de las secuencias.

## Resultados

El programa generará una matriz de dot plot y la visualizará en una imagen guardada en el archivo `filtered_dotplot.svg`. Además, mostrará por consola el tamaño de la matriz y el tiempo de ejecución.

¡Disfruta del análisis de dot plot!
