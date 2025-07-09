from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Função para ler o arquivo fasta
def read_fasta(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
    seq1 = lines[1].strip()
    seq2 = lines[3].strip()
    return seq1, seq2

# Função para calcular o escore usando uma matriz de substituição
def calculate_score(seq1, seq2, matrix):
    total_score = 0
    match_count = 0
    scores = []
    for a, b in zip(seq1, seq2):
        score = matrix.get((a, b), 0)
        scores.append(score)
        total_score += score
        if score > 0:
            match_count += 1
    return total_score, scores, match_count

# Definindo as sequências
seq1, seq2 = read_fasta('Seq1Seq2.txt')

# Matriz Identidade
identity_matrix = {
    ('A', 'A'): 1, ('A', 'T'): 0, ('A', 'C'): 0, ('A', 'G'): 0, ('A', '-'): 0,
    ('T', 'A'): 0, ('T', 'T'): 1, ('T', 'C'): 0, ('T', 'G'): 0, ('T', '-'): 0,
    ('C', 'A'): 0, ('C', 'T'): 0, ('C', 'C'): 1, ('C', 'G'): 0, ('C', '-'): 0,
    ('G', 'A'): 0, ('G', 'T'): 0, ('G', 'C'): 0, ('G', 'G'): 1, ('G', '-'): 0,
    ('-', 'A'): 0, ('-', 'T'): 0, ('-', 'C'): 0, ('-', 'G'): 0, ('-', '-'): 0,
}

# Matriz de Escore
score_matrix = {
    ('A', 'A'): 2, ('A', 'T'): -1, ('A', 'C'): -1, ('A', 'G'): -2, ('A', '-'): -2,
    ('T', 'A'): -1, ('T', 'T'): 2, ('T', 'C'): -2, ('T', 'G'): -1, ('T', '-'): -2,
    ('C', 'A'): -1, ('C', 'T'): -2, ('C', 'C'): 2, ('C', 'G'): -1, ('C', '-'): -2,
    ('G', 'A'): -2, ('G', 'T'): -1, ('G', 'C'): -1, ('G', 'G'): 2, ('G', '-'): -2,
    ('-', 'A'): -2, ('-', 'T'): -2, ('-', 'C'): -2, ('-', 'G'): -2, ('-', '-'): 0,
}

# Cálculo do escore para a Matriz Identidade
identity_total_score, identity_scores, identity_matches = calculate_score(seq1, seq2, identity_matrix)
identity_percentage = (identity_matches / len(seq1)) * 100

# Cálculo do escore para a Matriz de Escore
score_total_score, score_scores, _ = calculate_score(seq1, seq2, score_matrix)
score_average = score_total_score / len(seq1)

# Alinhamento utilizando o módulo Bio.pairwise2
alignments_identity = pairwise2.align.globalds(seq1, seq2, identity_matrix, -2, -2)
alignments_score = pairwise2.align.globalds(seq1, seq2, score_matrix, -2, -2)

# Exibindo resultados
print("1) Cálculos utilizando as Matrizes de Identidade e Escore:")

print("\nMatriz Identidade:")
print("Escore Total:", identity_total_score)
print("Identidade (%):", identity_percentage)

print("\nMatriz de Escore:")
print("Escore Total:", score_total_score)
print("Escore Médio por Nucleotídeo:", score_average)

print("\n2) Alinhamento utilizando Bio.pairwise2:")
print("\nAlinhamento com Matriz Identidade:")
for alignment in alignments_identity:
    print(format_alignment(*alignment))

print("\nAlinhamento com Matriz de Escore:")
for alignment in alignments_score:
    print(format_alignment(*alignment))

# Comparação dos escores
bio_identity_score = alignments_identity[0][2]
bio_score_matrix_score = alignments_score[0][2]

print("\nComparação dos Escores:")
print(f"Escore Matriz Identidade (Manual): {identity_total_score}")
print(f"Escore Matriz Identidade (Bio.pairwise2): {bio_identity_score}")

print(f"Escore Matriz de Escore (Manual): {score_total_score}")
print(f"Escore Matriz de Escore (Bio.pairwise2): {bio_score_matrix_score}")
