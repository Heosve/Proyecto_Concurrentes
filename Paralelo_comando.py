import sys
import argparse
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import time
from tqdm import tqdm
from Bio import SeqIO


def draw_dotplot(matrix, fig_name='dotplot.svg'):
    plt.figure(figsize=(5, 5))
    plt.imshow(matrix, cmap='Greys', aspect='auto')
    plt.ylabel("Secuencia 1")
    plt.xlabel("Secuencia 2")
    plt.savefig(fig_name)


def merge_sequences_from_fasta(file_path, max_length=None):
    sequences = []  # List to store all sequences
    for record in SeqIO.parse(file_path, "fasta"):
        seq = str(record.seq)
        if max_length is not None:
            seq = seq[:max_length]  # Limitar la longitud de la secuencia
        sequences.append(seq)
    return "".join(sequences)


def worker(args):
    i, Secuencia1, Secuencia2 = args
    return [Secuencia1[i] == Secuencia2[j] for j in range(len(Secuencia2))]


def parallel_dotplot(Secuencia1, Secuencia2, threads=mp.cpu_count()):
    with mp.Pool(processes=threads) as pool:
        results = []
        total_tasks = len(Secuencia1) * len(Secuencia2)
        with tqdm(total=total_tasks, desc="Progress", unit="task") as pbar:
            for i, result in enumerate(pool.imap_unordered(worker, [(i, Secuencia1, Secuencia2) for i in range(len(Secuencia1))])):
                results.append(result)
                pbar.update(len(Secuencia2))
    return results


if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate dot plot of two DNA sequences.')
    parser.add_argument('file1', type=str, help='Path to the first sequence file in FASTA format')
    parser.add_argument('file2', type=str, help='Path to the second sequence file in FASTA format')
    parser.add_argument('--max_length', type=int, help='Maximum length of the sequences')
    args = parser.parse_args()

    file_path_1 = args.file1
    file_path_2 = args.file2

    merged_sequence_1 = merge_sequences_from_fasta(file_path_1, max_length=args.max_length)
    merged_sequence_2 = merge_sequences_from_fasta(file_path_2, max_length=args.max_length)

    begin = time.time()
    dotplot = np.array(parallel_dotplot(merged_sequence_1, merged_sequence_2, threads=4))

    # Redimensionar la matriz de resultados
    matrix_size = (len(merged_sequence_1), len(merged_sequence_2))
    dotplot = np.array(dotplot).reshape(matrix_size)

    print("La matriz de resultado tiene tamaño:", dotplot.shape)
    print(f"\nEl código se ejecutó en: {time.time() - begin} segundos")

    ## Vamos a visualizar el dotplot
    draw_dotplot(dotplot, fig_name='dotplot_paralelizado.svg')
    draw_dotplot(dotplot[:500, :500])