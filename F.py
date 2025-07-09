def smith_waterman_score_matrix(seq1, seq2, match=2, mismatch=-1, gap=-2):
    # Inicialização da matriz de escore
    m, n = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    # Preenchimento da matriz de escore
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                score = match
            else:
                score = mismatch
            diag = score_matrix[i - 1][j - 1] + score
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(0, diag, up, left)

    return score_matrix

# Função para exibir a matriz de escore
def print_score_matrix(score_matrix, seq1, seq2):
    m, n = len(seq1), len(seq2)

    # Imprimir a primeira linha (seq2)
    print("      ", "  ".join(seq2))
    print("   " + " ".join(["{:3}".format(score_matrix[0][j]) for j in range(n + 1)]))
    
    # Imprimir as demais linhas (seq1)
    for i in range(1, m + 1):
        print(seq1[i - 1] + " " + " ".join(["{:3}".format(score_matrix[i][j]) for j in range(n + 1)]))

# Função para o alinhamento das sequências
def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    score_matrix = smith_waterman_score_matrix(seq1, seq2, match, mismatch, gap)
    max_score = 0
    max_pos = None

    # Encontrar a posição do máximo valor na matriz
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if score_matrix[i][j] >= max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Rastreamento para encontrar o alinhamento ótimo
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = max_pos

    while score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        diag = score_matrix[i - 1][j - 1]
        up = score_matrix[i - 1][j]
        left = score_matrix[i][j - 1]

        if current_score == diag + (match if seq1[i - 1] == seq2[j - 1] else mismatch):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif current_score == up + gap:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:  # current_score == left + gap
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1.reverse()
    aligned_seq2.reverse()

    return ''.join(aligned_seq1), ''.join(aligned_seq2), max_score

# Sequências fornecidas
seq1 = "GAATTCAGTTA"
seq2 = "GGATCGA"

# Calculando a matriz de escore
score_matrix = smith_waterman_score_matrix(seq1, seq2)

# Exibindo a matriz de escore
print("Matriz de Escore do Alinhamento Local:")
print_score_matrix(score_matrix, seq1, seq2)

# Alinhamento das sequências
aligned_seq1, aligned_seq2, score = smith_waterman(seq1, seq2)

print("\nSequências Alinhadas:")
print(f"Sequência 1: {aligned_seq1}")
print(f"Sequência 2: {aligned_seq2}")
print(f"Pontuação de alinhamento: {score}")
