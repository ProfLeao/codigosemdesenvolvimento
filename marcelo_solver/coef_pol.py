import numpy as np

def cfcond_termica(
    temp, 
    unid = 'K',
    temps_inf_interval = [20.,800.], # intervalo de baixas temperaturas
    temps_sup_interval = [800.,1515.] # intervalo de temperaturas superiores
):
    """
        Determina o coeficiente de condutividade térmica para temperaturas 
        ou vetores de temperaturas.

        Recebe:
        · temp(float ou ndarray) - temperaturas de cálculo. 
        · [unid(str)] - caractere 'K' ou 'C' para determinar a unidade de 
        temperatura 
        · [temps_inf_interval(list)] - lista de dois elementos com as 
        temperaturas inferior e superior para os polinômios de condutividade 
        térmica.
        · [temps_sup_interval(list)] - lista de dois elementos com as 
        temperaturas inferior e superior para os polinômios de condutividade 
        térmica. 

        Retorna:
        · ndarray ou float  de condutividades térmicas em W/m.K
    """
    if unid.lower() == 'c':
        temp += 273.15
    elif unid.lower() == 'k':
        pass
    else:
        raise ValueError(
            f"ERRO!\nA unidade de medida {unid} é desconhecida."
        )

    if str(type(temp)) == "<class 'numpy.ndarray'>":
        any_interval = [False, False]
        try:
            # intevalo 1: 20 <= temp <= 800
            temps_bx = temp[
                np.greater_equal(
                    temp, [temps_inf_interval[0] + 273.15]
                )
            ]
            temps_bx = temps_bx[
                np.less_equal(
                    temp, [temps_inf_interval[1] + 273.15]
                )
            ]
            coef_k = 54. - 3.33 * np.power(temps_bx, [2])
            any_interval[0] = True
        except:
            pass
        try:
            # intevalo 1: 800 < temp <= 1515
            temps_alt = temp[
                np.greater(
                    temp, [temps_sup_interval[0] + 273.15]
                )
            ]
            temps_alt = temps_bx[
                np.less_equal(
                    temp, [temps_sup_interval[1] + 273.15]
                )
            ]
            coef_k = 27.3
            any_interval[1] = True             
        except:
            pass

        if not True in any_interval:
            raise ValueError(
                "ERRO!\nIntervalo de temperatura indevido."
                )
    elif type(temp) in [int, float, complex]:
        if temp >=  temps_inf_interval[0] and temp <=  temps_inf_interval[1]:
            coef_k = 54. - 3.33 * np.power(temp, [2])
        elif temp > temps_sup_interval[0] and temp <=  temps_sup_interval[1]:
            coef_k = 27.3
    
    return coef_k