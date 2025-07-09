def smith_waterman(seq1, seq2, match=2, mismatch=-1, gap=-2):
    # Inicialização da matriz de pontuação
    m, n = len(seq1), len(seq2)
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    max_score = 0
    max_pos = None

    # Preenchimento da matriz de pontuação
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

# Exemplo de uso
seq1 = "ACACACTA"
seq2 = "AGCACACA"

aligned_seq1, aligned_seq2, score = smith_waterman(seq1, seq2)

print(f"Sequência 1 alinhada: {aligned_seq1}")
print(f"Sequência 2 alinhada: {aligned_seq2}")
print(f"Pontuação de alinhamento: {score}")
