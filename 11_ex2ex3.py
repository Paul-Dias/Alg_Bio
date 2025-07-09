def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    # Comprimento das sequências
    m, n = len(seq1), len(seq2)
    
    # Inicialização da matriz de escores
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Preenchimento da primeira linha e primeira coluna com penalidades de gap
    for i in range(1, m + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap
    for j in range(1, n + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap
    
    # Preenchimento da matriz de escores
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # Pontuação para match/mismatch
            if seq1[i - 1] == seq2[j - 1]:
                score = match
            else:
                score = mismatch
            # Calcula o valor da célula considerando diagonal, cima e esquerda
            diagonal = score_matrix[i - 1][j - 1] + score
            up = score_matrix[i - 1][j] + gap
            left = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(diagonal, up, left)
    
    # Função para exibir a matriz de escores
    def print_score_matrix():
        print("     " + " ".join(seq2))
        for i in range(m + 1):
            if i == 0:
                print(" ", end=" ")
            else:
                print(seq1[i - 1], end=" ")
            print(" ".join(f"{score_matrix[i][j]:2}" for j in range(n + 1)))
    
    # Alinhamento usando traceback
    def traceback():
        align1, align2 = [], []
        i, j = m, n
        while i > 0 or j > 0:
            current_score = score_matrix[i][j]
            if i > 0 and j > 0 and score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch) == current_score:
                align1.append(seq1[i - 1])
                align2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif i > 0 and score_matrix[i - 1][j] + gap == current_score:
                align1.append(seq1[i - 1])
                align2.append('-')
                i -= 1
            else:
                align1.append('-')
                align2.append(seq2[j - 1])
                j -= 1
        align1.reverse()
        align2.reverse()
        return ''.join(align1), ''.join(align2)
    
    # Exibindo a matriz de escores
    print("Matriz de Escore:")
    print_score_matrix()
    
    # Calculando e exibindo o alinhamento
    aligned_seq1, aligned_seq2 = traceback()
    print("\nSequências Alinhadas:")
    print(f"Sequência 1: {aligned_seq1}")
    print(f"Sequência 2: {aligned_seq2}")
    print(f"Pontuação de alinhamento: {score_matrix[m][n]}")

# Sequências fornecidas
seq1_1 = "AGCT"
seq2_1 = "ATGCT"

seq1_2 = "GATTACA"
seq2_2 = "GATACCA"

# Executando o alinhamento global para o primeiro par de sequências
print("Alinhamento para seq1: AGCT e seq2: ATGCT")
needleman_wunsch(seq1_1, seq2_1)

# Linha de separação
print("\n" + "="*40 + "\n")

# Executando o alinhamento global para o segundo par de sequências
print("Alinhamento para seq1: GATTACA e seq2: GATACCA")
needleman_wunsch(seq1_2, seq2_2)

