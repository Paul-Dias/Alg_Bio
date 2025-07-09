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

# Exemplo de uso
seq1 = "ACACACTA"
seq2 = "AGCACACA"

score_matrix = smith_waterman_score_matrix(seq1, seq2)

print("Matriz de Escore do Alinhamento Local:")
print_score_matrix(score_matrix, seq1, seq2)
