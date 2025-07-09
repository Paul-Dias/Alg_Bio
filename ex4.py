from Bio.SeqIO.FastaIO import FastaIterator
import os
import math
from Bio.Seq import Seq

# File paths
file_path = os.getcwd() + "/dengueIN.fasta"
output_file_path = os.getcwd() + "/fastaIN.txt"

# Dicionário de codificação de aminoácidos baseado em RNA
rna_codon_table = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W',
}

# Função para traduzir RNA em proteína considerando todos os frames
def translate_rna_to_protein_all_frames(rna_sequence):
    protein_sequences = []
    reverse_rna_sequence = rna_sequence.reverse_complement()
    
    for frame in range(3):
        protein_sequence = translate_rna_to_protein(rna_sequence[frame:])
        protein_sequences.append(protein_sequence)
        
        reverse_protein_sequence = translate_rna_to_protein(reverse_rna_sequence[frame:])
        protein_sequences.append(reverse_protein_sequence)
        
    return protein_sequences



# Função para traduzir RNA em proteína
def translate_rna_to_protein(rna_sequence):
    protein_sequence = ""
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        if len(codon) == 3:
            protein_sequence += rna_codon_table.get(codon.upper(), 'X')
    return protein_sequence


# Listas para armazenar sequências de aminoácidos em diferentes frames
protein_sequences_frames = [[] for _ in range(6)]  # Expanda o tamanho da lista para 18

# Lists to store GC content percentages and melting temperatures
gc_content_list = []
melting_temp_list = []

# Calculate melting temperature
def calculate_melting_temp(gc_content, dna_length):
    # Constants
    Na_conc = 100  # mM
    tm = 81.5 + 16.6 * math.log10(Na_conc) + 0.41 * gc_content - (500 / dna_length)
    return tm

# Read fasta file and calculate GC content and melting temperature
with open(output_file_path, 'w', encoding='utf-8') as output_file:
    with open(file_path) as handle:
        for values in FastaIterator(handle):
            fasta_sequence = values.seq
            fasta_dna = Seq(values.seq).transcribe()
            quant_a = fasta_sequence.count("A")
            quant_c = fasta_sequence.count("C")
            quant_t = fasta_sequence.count("T")
            quant_g = fasta_sequence.count("G")
            c_g = quant_g + quant_c
            tot = quant_g + quant_c + quant_a + quant_t
            gc_percentage = (c_g / tot) * 100
            gc_content_list.append(gc_percentage)
            melting_temp_list.append(calculate_melting_temp(gc_percentage, len(fasta_sequence)))
            protein_sequences = translate_rna_to_protein_all_frames(fasta_dna)
            output_file.writelines("\n \nArquivo fasta:\n")
            output_file.writelines(f"Sequência em DNA: \n")
            output_file.writelines(f"{fasta_sequence}\n")
            output_file.writelines(f"Sequencia em RNA traduzida da sequência de DNA: \n")
            output_file.writelines(f"{fasta_dna} \n")
            output_file.writelines(f"Sequências de proteína traduzidas do RNA: \n")
            for frame, seq in enumerate(protein_sequences):
                output_file.writelines(f"{seq}\n")
                protein_sequences_frames[frame].append(seq)  # Append the sequence to the corresponding frame


# Salva as sequências de proteína em arquivos diferentes
for frame, sequences in enumerate(protein_sequences_frames):
    with open(f"frame_{frame+1}_protein_sequence.txt", 'w', encoding='utf-8') as frame_file:
        frame_file.writelines(sequences)
