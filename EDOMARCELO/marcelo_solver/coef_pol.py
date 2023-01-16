import numpy as np

def cfcond_termica(
    temp, 
    unid = 'K',
    temps_inf_interval = [293.15, 1073.15], # intervalo de baixas temperaturas
    temps_sup_interval = [1073.15, 1788.15] # intervalo de temperaturas superiores
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
    # Type casting para as diferentes entradas de tipo. 
    try:
        temp = float(temp)
    except TypeError:
        temp = np.array(temp, dtype=np.float64)

    # Verifica mudança nos paraâmetros default
    change_def1 = not (
        temps_inf_interval == cfcond_termica.__defaults__[1]
    )
    change_def2 = not (
        temps_sup_interval == cfcond_termica.__defaults__[2]
    )
    change_def = change_def1 and change_def2

    temps_inf_interval =  np.array(temps_inf_interval, dtype=np.float64)
    temps_sup_interval =  np.array(temps_sup_interval, dtype=np.float64)

    if unid.lower() == 'c' and change_def:
        temp += 273.15
        temps_inf_interval += 273.15
        temps_sup_interval += 273.15
    elif unid.lower() == 'c' and change_def1:
        temp += 273.15
        temps_inf_interval += 273.15
    elif unid.lower() == 'c' and change_def2:
        temp += 273.15
        temps_sup_interval += 273.15
    elif unid.lower() == 'c':
        temp += 273.15
    elif unid.lower() == 'k':
        pass
    else:
        raise ValueError(
            f"ERRO!\nA unidade de medida {unid} é desconhecida."
        )

    # Aplica as equações segundo o tipo e o intervalo de temperatura
    if str(type(temp)) == "<class 'numpy.ndarray'>":
        any_interval = [False, False]

        # intevalo 1: 20 <= temp <= 800
        temps_bx = temp[
            np.greater_equal(
                temp, [temps_inf_interval[0]]
            )
        ]
        temps_bx = temps_bx[
            np.less_equal(
                temps_bx, [temps_inf_interval[1]]
            )
        ]
        coef_k1 = 54. - 3.33e-2 * (temps_bx-273.15)
        if coef_k1.size != 0: any_interval[0] = True

        # intevalo 1: 800 < temp <= 1515
        temps_alt = temp[
            np.greater(
                temp, [temps_sup_interval[0]]
            )
        ]
        temps_alt = temps_alt[
            np.less_equal(
                temps_alt, [temps_sup_interval[1]]
            )
        ]
        coef_k2 = np.ones_like(temps_alt)
        coef_k2 = coef_k2 * 27.3 
        if coef_k2.size != 0: any_interval[1] = True            

        if not True in any_interval:
            raise ValueError(
                "ERRO!\nIntervalo de temperatura indevido."
                )
        
        temps_coefs = np.array(
            [
                np.concatenate([temps_bx, temps_alt]),
                np.concatenate([coef_k1, coef_k2])
            ]
        )

        return temps_coefs
            
    elif type(temp) in [int, float, complex]:
        if temp >=  temps_inf_interval[0] and temp <=  temps_inf_interval[1]:
            coef_k = 54. - 3.33e-2 * (temp-273.15)
        elif temp > temps_sup_interval[0] and temp <=  temps_sup_interval[1]:
            coef_k = 27.3
        else:
            raise ValueError(
                "ERRO!\n"+\
                f"A temperatura {temp} não é definida no intervalo de"+\
                "temperaturas do modelo:"+\
                f"\n{temps_inf_interval}\n{temps_sup_interval}"
            )
    return coef_k

def calesp_vol(
    temp, 
    unid = 'K',
    temps_itv_1 = [20. + 273.15, 600. + 273.15], 
    temps_itv_2 = [600. + 273.15, 735. + 273.15],
    temps_itv_3 = [735. + 273.15, 900. + 273.15],
    temps_itv_4 = [900. + 273.15, 1515. + 273.15]
):
    """
        Determina o calor específico volumar para temperaturas (escalar)
        ou vetores de temperaturas.

        Recebe:
        · temp(float ou ndarray) - temperaturas de cálculo. 
        · [unid(str)] - caractere 'K' ou 'C' para determinar a unidade de 
        temperatura 
        · [temps_itv_1(list)] - lista de dois elementos com as 
        temperaturas inferior e superior para o primeiro polinômio de cv.
        · [temps_itv_2(list)] - lista de dois elementos com as 
        temperaturas inferior e superior para o segundo polinômio de cv.
        · [temps_itv_3(list)] - lista de dois elementos com as 
        temperaturas inferior e superior para o terceiro polinômio de cv.
        · [temps_itv_4(list)] - lista de dois elementos com as 
        temperaturas inferior e superior para o quarto polinômio de cv.

        Retorna:
        · ndarray ou float  de calores específicos volumares em K/kg.K
    """
    # Type casting para as diferentes entradas de tipo.
    try:
        temp = float(temp)
    except TypeError:
        temp = np.array(temp, dtype=np.float64)
    
    # Verifica mudança nos parâmetros default
    change_def = []
    temps_itrvs = [
        temps_itv_1, temps_itv_2,
        temps_itv_3, temps_itv_4
    ]

    for tp, pd in zip(temps_itrvs, cfcond_termica.__defaults__[1:]):
        change_def.append(not (tp == pd))

    change_def.append(
        (change_def[0] and change_def[1]) and 
        (change_def[2] and change_def[3])
    ) # O último registro da lista marca a 
      # mudança de todos os parâmetros

    if unid.lower() == 'c' and change_def[-1]:
        temp += 273.15
        for tidx in range(len(temps_itrvs)):
            temps_itrvs[tidx] = temps_itrvs[tidx] + 273.15
    elif unid.lower() == 'c' and change_def[0]:
        temp += 273.15
        temps_itrvs[0] = temps_itrvs[0] + 273.15
    elif unid.lower() == 'c' and change_def[1]:
        temp += 273.15
        temps_itrvs[1] = temps_itrvs[1] + 273.15
    elif unid.lower() == 'c' and change_def[2]:
        temp += 273.15
        temps_itrvs[2] = temps_itrvs[2] + 273.15
    elif unid.lower() == 'c' and change_def[3]:
        temp += 273.15
        temps_itrvs[3] = temps_itrvs[3] + 273.15
    elif unid.lower() == 'k':
        pass
    else:
        raise ValueError(
            f"ERRO!\nA unidade de medida {unid} é desconhecida."
        )
    
    # Aplica as equações segundo o tipo e o intervalo de temperatura
    if str(type(temp)) == "<class 'numpy.ndarray'>":
        any_interval = [False for l in range(len(temps_itrvs))]

        # intevalo 1: 20 <= temp <= 600
        temps_0 = temp[
            np.greater_equal(
                temp, [temps_itrvs[0][0]]
            )
        ]
        temps_0 = temps_1[
            np.less_equal(
                temps_0, [temps_itrvs[0][1]]
            )
        ]
        cv0 = 425. + 7.73e-1 * temps0 -\
            1.69e-3 * np.power(temps0,2) -\
            2.22e-6 * np.power(temps0,3)
        if cv0.size != 0: any_interval[0] = True

        # intevalo 1: 600 < temp <= 735
        temps_1 = temp[
            np.greater(
                temp, [temps_itrvs[1][0]]
            )
        ]
        temps_1 = temps_1[
            np.less_equal(
                temps_1, [temps_itrvs[1][1]]
            )
        ]
        cv1 = 425. + 13002/(738. - temps_1)
        if cv1.size != 0: any_interval[1] = True

        # intevalo 1: 735 < temp <= 900
        temps_2 = temp[
            np.greater(
                temp, [temps_itrvs[2][0]]
            )
        ]
        temps_2 = temps_2[
            np.less_equal(
                temps_2, [temps_itrvs[2][1]]
            )
        ]
        cv2 = 545. + 17820./(temps_2 - 731.)
        if cv2.size != 0: any_interval[2] = True

        # intevalo 1: 900 < temp <= 1515
        temps_3 = temp[
            np.greater(
                temp, [temps_itrvs[3][0]]
            )
        ]
        temps_3 = temps_3[
            np.less_equal(
                temps_3, [temps_itrvs[3][1]]
            )
        ]
        cv3 = 650.
        if cv3.size != 0: any_interval[3] = True

        if not True in any_interval:
            raise ValueError(
                "ERRO!\nIntervalo de temperatura indevido."
                )
        return np.concatenate((cv0, cv1, cv2, cv3))
    elif type(temp) == float:
        if temp >=  temps_itrvs[0][0] and temp <=  temps_itrvs[0][1]:
            # intevalo 1: 20 <= temp <= 600
            cv = 425. + 7.73e-1 * temps0 -\
            1.69e-3 * np.power(temps0,2) -\
            2.22e-6 * np.power(temps0,3)
        elif temp >  temps_itrvs[1][0] and temp <=  temps_itrvs[1][1]:
            # intevalo 1: 600 < temp <= 735
            cv = 425. + 13002/(738. - temps_1)
        elif temp >  temps_itrvs[1][0] and temp <=  temps_itrvs[1][1]:
            # intevalo 1: 735 < temp <= 900
            cv = 545. + 17820./(temps_2 - 731.)
        elif temp >  temps_itrvs[1][0] and temp <=  temps_itrvs[1][1]:
            # intevalo 1: 900 < temp <= 1515
            cv =650.
        else:
            raise ValueError(
                "ERRO!\n"+\
                f"A temperatura {temp} não é definida no intervalo de"+\
                "temperaturas do modelo."
            )
    return cv