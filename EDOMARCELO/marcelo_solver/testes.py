from .coef_pol import *
import matplotlib.pyplot as plt

def teste_condterm(teste):
    """
        Realiza testes de funcionamento na função de determinação do 
        coeficiente de condutividade térmica.

        1° - Teste de operação sobre vetor de temperaturas válidas:
        Recebe:
        teste = "VECTEMP-VAL"
            Este método cria um vetor no intervalo de temperaturas válido para
            para o modelo de coeficiente de condutividade térmica.
            Espera-se o cálculo correto dos valores de h expressos de forma gráfica.

        2° - Teste de operação sobre vetor de temperaturas com valores inválidos:
        teste = "VECTEMP-INVAL"
            Este método cria um vetor com valores fora do intervalo de 
            temperaturas válido para para o modelo de coeficiente de 
            condutividade térmica.
            Espera-se o levantamento de erro durante a tentativa de execuação.

        3° - Teste de operação sobre escalares randômicos válidos:
        teste = "RANDTEMP-VAL"
            Neste método de teste, um gerador de números aleatórios é 
            recorrentemente chamado para gerar N valores no intervalo de 
            operação da função.
            Espera-se o cálculo correto dos valores de h expressos de forma gráfica.
    """
    print("Testes na função cfcond_termica().")

    if teste == "VECTEMP-VAL":
        temp = np.arange(20.,1515.,0.1)
        matriz = cfcond_termica(temp)
        fig, ax = plt.subplots()
        ax.plot(matriz[0], matriz[1])
        return fig, ax
    elif teste == "VECTEMP-PARC-INVAL":
        temp = np.arange(0.,2000.,0.1)
        print(f"Testando um vetor de temperaturas entre {temp[0]} e {temp[-1]}")
        matriz = cfcond_termica(temp)
        fig, ax = plt.subplots()
        ax.plot(matriz[0], matriz[1])
        return fig, ax
    elif teste == "RANDTEMP-VAL":
        temp = np.random.randint(20,1515, size=100)
        hs = np.empty_like(temp)
        for i,t in enumerate(temp):
            hs[i] = cfcond_termica(t, 'c')

        fig, ax  = plt.subplots()
        ax.plot(temp, hs, ".b")
        return fig, ax

def teste_cv(teste):
    """
        Realiza testes de funcionamento na função cálculo do calor específico
        volumar.

        1° - Teste de operação sobre vetor de temperaturas válidas:
        Recebe:
        teste = "VECTEMP-VAL"
            Este método cria um vetor no intervalo de temperaturas válido para
            para o modelo de coeficiente de condutividade térmica.
            Espera-se o cálculo correto dos valores de h expressos de forma gráfica.

        2° - Teste de operação sobre vetor de temperaturas com valores inválidos:
        teste = "VECTEMP-INVAL"
            Este método cria um vetor com valores fora do intervalo de 
            temperaturas válido para para o modelo de coeficiente de 
            condutividade térmica.
            Espera-se o levantamento de erro durante a tentativa de execuação.

        3° - Teste de operação sobre escalares randômicos válidos:
        teste = "RANDTEMP-VAL"
            Neste método de teste, um gerador de números aleatórios é 
            recorrentemente chamado para gerar N valores no intervalo de 
            operação da função.
            Espera-se o cálculo correto dos valores de h expressos de forma gráfica.
    """
    print("Testes na função teste_cv().")

    if teste == "VECTEMP-VAL":
        temp = np.arange(20.,1515.,0.1)
        matriz = calesp_vol(temp)
        fig, ax = plt.subplots()
        ax.plot(matriz[0], matriz[1])
        return fig, ax
    elif teste == "VECTEMP-PARC-INVAL":
        temp = np.arange(0.,2000.,0.1)
        print(f"Testando um vetor de temperaturas entre {temp[0]} e {temp[-1]}")
        matriz = calesp_vol(temp)
        fig, ax = plt.subplots()
        ax.plot(matriz[0], matriz[1])
        return fig, ax
    elif teste == "RANDTEMP-VAL":
        temp = np.random.randint(20,1515, size=200)
        cv = np.empty_like(temp)
        for i,t in enumerate(temp):
            cv[i] = calesp_vol(t)

        fig, ax  = plt.subplots()
        ax.plot(temp, cv, ".b")
        return fig, ax

def teste_dens(teste):
    """
        Realiza testes de funcionamento na função cálculo da densidade.

        1° - Teste de operação sobre vetor de temperaturas válidas:
        Recebe:
        teste = "VECTEMP-VAL"
            Este método cria um vetor no intervalo de temperaturas válido para
            para o modelo de coeficiente de condutividade térmica.
            Espera-se o cálculo correto dos valores de h expressos de forma gráfica.

        2° - Teste de operação sobre vetor de temperaturas com valores inválidos:
        teste = "VECTEMP-INVAL"
            Este método cria um vetor com valores fora do intervalo de 
            temperaturas válido para para o modelo de coeficiente de 
            condutividade térmica.
            Espera-se o levantamento de erro durante a tentativa de execuação.

        3° - Teste de operação sobre escalares randômicos válidos:
        teste = "RANDTEMP-VAL"
            Neste método de teste, um gerador de números aleatórios é 
            recorrentemente chamado para gerar N valores no intervalo de 
            operação da função.
            Espera-se o cálculo correto dos valores de h expressos de forma gráfica.
    """
    print("Testes na função densidade().")

    if teste == "VECTEMP-VAL":
        temp = np.arange(20.,1515.,0.1)
        matriz = densidade(temp)
        fig, ax = plt.subplots()
        ax.plot(matriz[0], matriz[1])
        return fig, ax
    elif teste == "VECTEMP-PARC-INVAL":
        temp = np.arange(0.,2000.,0.1)
        print(f"Testando um vetor de temperaturas entre {temp[0]} e {temp[-1]}")
        matriz = densidade(temp)
        fig, ax = plt.subplots()
        ax.plot(matriz[0], matriz[1])
        return fig, ax
    elif teste == "RANDTEMP-VAL":
        temp = np.random.randint(25,1515, size=200)
        rho = np.empty_like(temp)
        for i,t in enumerate(temp):
            rho[i] = densidade(t)
        fig, ax  = plt.subplots()
        ax.plot(temp, rho, ".b")
        return fig, ax
