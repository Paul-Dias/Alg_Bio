# A) Algoritmo

# declarar 4 variáveis correspondentes a A,C,G,T
# Receber uma sequncia de strings
# a verificação e contagem se da utilizando o metodo count, o método lê toda a sequência e conta apenas as STRING selecionadas
# vericar se A,C,G,T estão na sequencia e conta-los individualmente
# somar o resultado individual e exibir o resultado total
# exibir o resultado de cada uma das letras 


#B) Programa/Código Python

def count_elemnt(string: str):
    if string == '':
        return "String Vazia"
    
    quant_a,quant_c,quant_g,quant_t = (0,0,0,0)

    for letra in string:
        if letra == 'A':
            quant_a +=1
        elif letra == 'C':
            quant_c +=1
        elif letra == 'T':
            quant_t +=1
        elif letra == 'G':
            quant_g+=1

    tot= quant_g+quant_a+quant_c+quant_t

    return f'A: {quant_a} \n C: {quant_c} \n T: {quant_t} \n G: {quant_g} \n total nucleotídeos: {tot}'


resp = input('Digite a sequência: ')
print(count_elemnt(resp))





# C) Resultado




