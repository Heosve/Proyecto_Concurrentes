#!~\anaconda3\python.exe

import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import argparse
import time

def create_dot_plot(seq1, seq2, window_size):
    length1 = len(seq1)  # Obtener la longitud de la secuencia 1
    length2 = len(seq2)  # Obtener la longitud de la secuencia 2

    tiempos =[]
    dot_plot = np.zeros((length1, length2), dtype=np.int8)  # Crear una matriz de ceros del tamaño adecuado para el dot plot

    for i in range(0, length1, window_size):
        for j in range(0, length2, window_size):
            start_time = time.time()
            start_i = i  # Inicio de la sección de la secuencia 1
            end_i = min(i + window_size, length1)  # Fin de la sección de la secuencia 1
            start_j = j  # Inicio de la sección de la secuencia 2
            end_j = min(j + window_size, length2)  # Fin de la sección de la secuencia 2

            section1 = np.array(list(seq1[start_i:end_i]))  # Obtener la sección de la secuencia 1 como un arreglo NumPy
            section2 = np.array(list(seq2[start_j:end_j]))  # Obtener la sección de la secuencia 2 como un arreglo NumPy

            match_matrix = section1[:, np.newaxis] == section2  # Comparar los elementos de las secciones y obtener una matriz de booleanos
            dot_plot[start_i:end_i, start_j:end_j] = match_matrix.astype(np.int8)  # Asignar la matriz de coincidencias al dot plot
            
            plt.imsave(f'img/dot_plot{start_i}-{end_i}vs{start_j}-{end_j}.png', dot_plot[start_i:end_i, start_j:end_j], cmap='gray')
            
            elapsed_times_bloques = time.time() - start_time
            tiempos.append(elapsed_times_bloques)
           

    return tiempos

if __name__ == '__main__':
    # Configurar la interfaz de línea de comandos
    parser = argparse.ArgumentParser(description='Crear un dot plot a partir de dos secuencias FASTA.')
    parser.add_argument('file1', type=str, help='Archivo FASTA de la secuencia 1')
    parser.add_argument('file2', type=str, help='Archivo FASTA de la secuencia 2')
    parser.add_argument('window_size', type=int, help='Tamaño de la ventana deslizante')
    args = parser.parse_args()

    # Lee las secuencias de los archivos fasta y limita su tamaño
    sequences1 = [str(record.seq[:40000]) for record in SeqIO.parse(args.file1, "fasta")]
    sequences2 = [str(record.seq[:40000]) for record in SeqIO.parse(args.file2, "fasta")]

    # Medir el tiempo de ejecución
    times_bloques = []


    # Itera sobre las secuencias
    for sequence1 in sequences1:
        for sequence2 in sequences2:
            start_time = time.time()

            # Calcula el dot plot por secciones para cada par de secuencias
            tiempos = create_dot_plot(sequence1, sequence2, args.window_size)
            # Realiza las operaciones necesarias con el dot plot obtenido
            
            plt.plot(range(1, len(tiempos) + 1), tiempos)
            plt.xlabel("Iteración de bloque")
            plt.ylabel("Tiempo transcurrido (segundos)")
            plt.title("Tiempos de ejecución por iteración del bucle sequence2")
            plt.savefig("img/tiemposBloque.png")
            break
        print(f"tiempo total : {time.time() - start_time}")
        break