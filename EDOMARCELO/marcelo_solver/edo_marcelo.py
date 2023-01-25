from .coef_pol import *
import numpy as np
from scipy.integrate import solve_ivp

def edo_marcelo(
    parm,
    cond_cont,
    interv,
    method = "RK45", # Runge-Kutta de ordem 5(4)
    f_k = cfcond_termica,
    f_cv = calesp_vol,
    f_rho = densidade,
    f_h = tc_convec
):
    """
        Recebe:
        
        · parm:tuple (
            coef_w:float,
            dic_h:dict, # ATENÇÃO AO TIPO DISTOANTE
            coef_Tinf:float, 
            coef_j:float, 
            coef_r:float,
            coef_delh:float,
            coef_V: float,
        ) - uma tupla de 7 parâmetros do tipo float para o modelo.
        
        · cond_cont:tuple (
            [cd00:float,cd01:float]:list, 
            [cd10:float,cd11:float]:list
        ) - uma tupla de duas listas com pares ordenados descritores das 
        condições de contorno. 
        
        · interv: tuple (if:float, sup:float, disc:float) - tupla contendo, 
        os limites inferior e superior de integração e o a discretização do 
        vetor.

        Por padrão a função usa Runge-Kutta de ordem 5(4) como método numérico.

        Retorna:
        
        Dois objetos scipy de grupo com os seguintes campos definidos em cada
        um.
        ·t:ndarray, shape (n_points,)
            Pontos de temperatura.
        ·y:ndarray, shape (n, n_points)
            Valores das soluções em cada t.
        ·sol:OdeSolution ou None
            Solução como um objeto OdeSolution; None se dense_output foi 
            passado como False.
        ·t_events:list de ndarray ou None
            Contém para cada tipo de evento uma lista de arrays nos quais um 
            evento desse tipo foi detectado. None se nenhum evento foi 
            encontrado.
        ·y_events:list de ndarray ou None
            Para cada valor de t_events, a solução correspondente. None se 
            nenhum evento foi encontrado.
        ·nfev:int
            Número de avaliações do lado direito.
        ·njev:int
            Número de avaliações do Jacobiano.
        ·nlu:int
            Número de decomposições LU.
        ·status:int
            Causas de término do algoritmo:
            -1: Falha no passo de integração.
            0: O solver alcançou com sucesso o final do intervalo de integração.
            1: Ocorrência de evento de término. 
        ·message:string
            Descrição legível para humanos do motivo da término.
        ·success:bool
            True se o solucionador atingiu o fim do intervalo ou ocorreu um evento de término (status >= 0).


    """
    status = None
    def modelo(v_z, v_temp):
        temp = v_temp[0]
        dtdz = v_temp[1]
        status = f"::-> Integrando a temperatura {temp:.3f}°C em {v_z:.3f} mm"
        print(status, end='\r')
        inv_k = 1/f_k(temp)
        parc1 = - parm[0] * f_rho(temp) * f_cv(temp) * dtdz
        parc2 = f_h(v_z, parm[1]) * (temp - parm[2])
        parc3 = - np.power(parm[3], 2) * parm[4]

        dvdz = inv_k*(parc1 + parc2 + parc3)

        return [temp, dvdz]

    dtdz0 = 1/cfcond_termica(cond_cont[0][1]) * (
        parm[0] * densidade(cond_cont[0][1]) * parm[5] -\
        parm[3] * parm[6]
    )
    #try:
    sol = solve_ivp(
        fun = modelo, 
        t_span = [interv[0], interv[1]],
        y0 = [cond_cont[0][1], dtdz0],
        t_eval= np.arange(interv[0], interv[1]+interv[2], interv[2]),
        method = "RK45",
        #max_step=interv[2]
    )
    return sol
    #except Exception as a:
    #    with open("ult.log", 'w') as arquivo:
    #        arquivo.write(
    #            "Erro:\n" + str(a)
    #        )
    #        arquivo.write(str(status)+'\n')


    