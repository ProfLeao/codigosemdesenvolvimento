import numpy as np

def edo_marcelo(
    coef_w,
    coef_h,
    coef_Tinf,
    coef_j,
    coef_r,
    f_k = cfcond_termica,
    f_cv = calesp_vol,
    f_rho =densidade
):
    def modelo(v_temp, v_z):
        # Primeira parcela
        parc1 = cfcond_termica(v_temp[0]) *\
             v_temp[1] + v_temp[0]*np.gradient(cfcond_termica(v_temp[0]))

        # Segunda parcela
        parc2 = coef_w * f_rho(v_temp[0]) *\
            (f_cv(v_temp[0]) * v_temp[1] +\
                v_temp[0] * np.gradient(f_cv(v_temp[0])))

        # Terceira parcela
        parc3 = - coef_h * (v_temp[0] - coef_Tinf)

        # Quarta parcela
        parc4 = np.power(coef_j, 2) * coef_r

        return parc1 + parc2 + parc3 + parc4
    
        
    
