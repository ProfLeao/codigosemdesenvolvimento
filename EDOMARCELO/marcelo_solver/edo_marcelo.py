from .coef_pol import *
import numpy as np
from scipy.integrate import solve_ivp

def edo_marcelo(
    tup_parm,
    par_cond_cont,
    interv,
    method = "RK45", # Runge-Kutta de ordem 5(4)
    f_k = cfcond_termica,
    f_cv = calesp_vol,
    f_rho = densidade
):
    """
        Recebe:
        
        · parm:tuple (
            coef_w:float,
            coef_h:float, 
            coef_Tinf:float, 
            coef_j:float, 
            coef_r:float
        ) - uma tupla de 5 parâmetros do tipo float para o modelo.
        
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

    #def modelo(v_z, v_temp):
    #    # Primeira parcela
    #    parc1 = cfcond_termica(v_temp[0]) *\
    #         v_temp[2] + v_temp[1]*np.gradient(cfcond_termica(v_temp[0]))
    #
    #    # Segunda parcela
    #    parc2 = tup_parm[0] * f_rho(v_temp[0]) *\
    #        (f_cv(v_temp[0]) * v_temp[1] +\
    #            v_temp[0] * np.gradient(f_cv(v_temp[0])))
    #
    #    # Terceira parcela
    #    parc3 = - tup_parm[1] * (v_temp[0] - tup_parm[2])
    #
    #    # Quarta parcela
    #    parc4 = np.power(tup_parm[3], 2) * tup_parm[4]
    #
    #    return parc1 + parc2 + parc3 + parc4
    #
    #sol0 = solve_ivp(
    #    fun = modelo, 
    #    t_span = interv[0:2],
    #    y0 = np.array([par_cond_cont[0][1]]),
    #    max_step = interv[2],
    #    method = "RK45"
    #)
    #sol1 = solve_ivp(
    #    fun = modelo, 
    #    t_span = np.flip(interv[0:2]),
    #    y0 = np.array([par_cond_cont[1][1]]),
    #    max_step = -interv[2],
    #    method = "RK45"
    #)

    #return sol0, sol1
    
    def modelo(v_z, v_temp):
        temp = v_temp[0]
        dtdz = v_temp[1]
        inv_k = 1/cfcond_termica(temp)
        parc1 = -tup_parm[0] * densidade(temp) * calesp_vol(temp) * dtdz
        parc2 = tup_parm[1] * (temp- tup_parm[2])
        parc3 = - np.power(tup_parm[3], 2) * tup_parm[4]

        dvdz = parc1 + parc2 + parc3

        return [temp, dvdz]

    sol0 = solve_ivp(
        fun = modelo, 
        t_span = [interv[0], interv[1]],
        y0 = [par_cond_cont[1][1], par_cond_cont[0][1]],
        t_eval= np.arange(interv[0], interv[1], interv[2]),
        method = "RK45"
    )
    sol1 = solve_ivp(
        fun = modelo, 
        t_span = [interv[1], interv[0]],
        y0 = [par_cond_cont[0][1], par_cond_cont[1][1]],
        t_eval= np.arange(interv[1], interv[0], -interv[2]),
        method = "RK45"
    )

    return sol0, sol1