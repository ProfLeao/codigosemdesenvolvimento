from .coef_pol import *
import numpy as np
from scipy.integrate import solve_ivp

def edo_marcelo(
    tup_parm,
    par_cond_cont,
    iterv,
    method = "RK45", # Runge-Kutta de ordem 5(4)
    f_k = cfcond_termica,
    f_cv = calesp_vol,
    f_rho = densidade
):
    """
        Recebe:
        
        · tup_parm:tuple (
            coef_w:float,
            coef_h:float, 
            coef_Tinf:float, 
            coef_j:float, 
            coef_r:float
        ) - uma tupla de 5 parâmetros do tipo float para o modelo.
        
        · par_cond_cont:tuple (
            [cd00:float,cd01:float]:list, 
            [cd10:float,cd11:float]:list
        ) - uma tupla de duas listas com pares ordenados descritores das 
        condições de contorno. 
        
        · iterv: tuple (if:float, sup:float, disc:float) - tupla contendo, 
        os limites inferior e superior de integração e o a discretização do 
        vetor.

        Por padrão a função usa Runge-Kutta de ordem 5(4) como método numérico.
    """

    def modelo(v_temp, v_z):
        # Primeira parcela
        parc1 = cfcond_termica(v_temp[0]) *\
             v_temp[2] + v_temp[1]*np.gradient(cfcond_termica(v_temp[0]))

        # Segunda parcela
        parc2 = tup_parm[0] * f_rho(v_temp[0]) *\
            (f_cv(v_temp[0]) * v_temp[1] +\
                v_temp[0] * np.gradient(f_cv(v_temp[0])))

        # Terceira parcela
        parc3 = - tup_parm[1] * (v_temp[0] - tup_parm[2])

        # Quarta parcela
        parc4 = np.power(tup_parm[3], 2) * tup_parm[4]

        return parc1 + parc2 + parc3 + parc4
    
    comprimento = np.arange(iterv[0], iterv[1]+iterv[2], iterv[2])
    dom = np.array([
        temperature,
        np.linspace(
            par_cond_cont[0][1], 
            par_cond_cont[1][1], 
            temperature.size
        )
    ])

    sol = solve_ivp(
        func = modelo, 
        t_span = interv[0:2],
        y0 = dom,
        method = "RK45"
    )

    return sol 
    
