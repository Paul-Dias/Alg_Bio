import re

def extrair_numeros(texto):
    # Expressão regular para encontrar números
    padrao = r'\d+'
    numeros = re.findall(padrao, texto)
    return numeros

def encontrar_region_numeros(arquivo):
    intervalos = []
    with open(arquivo, 'r') as f:
        for linha in f:
            if 'Region' in linha:
                numeros = extrair_numeros(linha)
                if len(numeros) == 2:
                    inicio, fim = map(int, numeros)
                    intervalos.append((inicio, fim))
                    print(f'Region encontrado: {inicio}:{fim}')
    return intervalos

def extrair_sequencia_intervalo(arquivo_sequencia, intervalos):
    sequencia = ''
    with open(arquivo_sequencia, 'r') as f:
        next(f)  # Ignora a primeira linha
        for linha in f:
            sequencia += linha.strip()

    sequencias = []
    for inicio, fim in intervalos:
        sequencia_intervalo = sequencia[inicio - 1:fim]
        sequencias.append(sequencia_intervalo)
    
    return sequencias

def gerar_arquivo_fasta(sequencias, arquivo_saida):
    with open(arquivo_saida, 'w') as f:
        for i, seq in enumerate(sequencias, start=1):
            f.write(f'>Seq_{i}\n')
            f.write(f'{seq}\n')

# Nome do arquivo contendo as regiões
arquivo = 'sequence.gp'

# Extrair os intervalos das regiões
intervalos = encontrar_region_numeros(arquivo)

# Arquivo contendo a sequência de letras
arquivo_sequencia = 'seqdump.txt'

# Extrair as sequências dentro dos intervalos
sequencias = extrair_sequencia_intervalo(arquivo_sequencia, intervalos)

# Gerar o arquivo FASTA com as sequências
arquivo_saida = 'sequencias_region.fasta'
gerar_arquivo_fasta(sequencias, arquivo_saida)
