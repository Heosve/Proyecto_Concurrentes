import time
import sys
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from Bio import SeqIO
from mpi4py import MPI


def merge_sequences_from_fasta(file_path):
    sequences = []  # List to store all sequences
    for record in SeqIO.parse(file_path, "fasta"):
        # record.seq gives the sequence
        sequences.append(str(record.seq))
    return "".join(sequences)


def crear_dotplot(args):
    secuencia1, secuencia2, indice = args
    codigos_secuencia1 = np.frombuffer(secuencia1.encode(), dtype=np.uint8)
    codigos_secuencia2 = np.frombuffer(secuencia2.encode(), dtype=np.uint8)

    dotplot = np.zeros((len(secuencia1), len(secuencia2)), dtype=np.uint8)
    for i in range(len(secuencia1)):
        matches = codigos_secuencia1[i] == codigos_secuencia2
        dotplot[i, matches] = 1
    return (indice, dotplot)


def procesar_comparacion(secuencia1, secuencia2):
    comm = MPI.COMM_WORLD
    num_procesos = comm.Get_size()
    rank = comm.Get_rank()

    if rank == 0:
        # Dividir secuencia1 en partes y enviar a los otros procesos
        subsecuencias1 = dividir_secuencia(secuencia1, num_procesos)
    else:
        subsecuencias1 = None

    # Cada proceso recibe su subsecuencia1 correspondiente
    subsecuencia1 = comm.scatter(subsecuencias1, root=0)

    # Cada proceso calcula su parte del dotplot
    resultado_parcial = crear_dotplot((subsecuencia1, secuencia2, rank))

    # Recopilar los resultados parciales en el proceso 0
    resultados = comm.gather(resultado_parcial, root=0)

    if rank == 0:
        dotplot = np.zeros((len(secuencia1), len(secuencia2)), dtype=np.uint8)

        for i, resultado in tqdm(enumerate(resultados), total=num_procesos):
            indice, resultado_parcial = resultado
            dotplot[indice * len(secuencia1) // num_procesos: (indice + 1) * len(secuencia1) // num_procesos] = resultado_parcial

        return dotplot
    else:
        return None


def dividir_secuencia(secuencia, num_partes):
    longitud_subsecuencia = len(secuencia) // num_partes
    subsecuencias = []
    for i in range(num_partes):
        inicio = i * longitud_subsecuencia
        fin = (i + 1) * longitud_subsecuencia
        subsecuencia = secuencia[inicio:fin]
        subsecuencias.append(subsecuencia)
    return subsecuencias


def calcular_peso_matriz(matriz):
    # Obtener el tamaño en bytes de la matriz
    bytes_matriz = sys.getsizeof(matriz)

    # Convertir a megabytes
    megabytes_matriz = bytes_matriz / (1024 ** 2)

    return megabytes_matriz


def draw_dotplot(matrix, fig_name='dotplot.svg'):
    plt.figure(figsize=(500, 500))
    plt.imshow(matrix, cmap='gray', aspect='auto')
    plt.ylabel("Secuencia 1")
    plt.xlabel("Secuencia 2")
    plt.savefig(fig_name)


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    comm_size = comm.Get_size()

    secuencia1 = merge_sequences_from_fasta('datos/E_coli.fna')
    secuencia2 = merge_sequences_from_fasta('datos/Salmonella.fna')
    seccion_matriz = 46000

    # Listas para almacenar los tiempos de ejecución y los números de procesos
    tiempos = []
    procesos = []

    for num_procesos in range(1, comm_size + 1):
        comm.Barrier()  # Sincronizar todos los procesos

        if rank == 0:
            inicio_tiempo = time.time()

        if num_procesos == 1:
            dotplot = procesar_comparacion(secuencia1[0:seccion_matriz], secuencia2[0:seccion_matriz])
        else:
            dotplot = None

        comm.Barrier()  # Sincronizar todos los procesos

        if rank == 0:
            fin_tiempo = time.time()
            tiempo_ejecucion = fin_tiempo - inicio_tiempo

            tiempos.append(tiempo_ejecucion)
            procesos.append(num_procesos)

            # Imprimir el tiempo de ejecución para el número actual de procesos
            print(f"Tiempo de ejecución con {num_procesos} proceso(s): {tiempo_ejecucion} segundos")

    if rank == 0:
        # Graficar los tiempos de ejecución
        plt.plot(procesos, tiempos)
        plt.xlabel("Número de procesos")
        plt.ylabel("Tiempo de ejecución (segundos)")
        plt.title("Rendimiento con diferentes números de procesos")
        plt.show()