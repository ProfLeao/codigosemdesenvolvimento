import numpy as np
import matplotlib.pyplot as plt
from marcelo_solver.edo_marcelo import edo_marcelo
from marcelo_solver.diversos import carrega_parametros
import pandas as pd

def main():
    print("-------------------------------------------------------------------------")
    print("|\tDistribuição de temperatura em eletrodo ao longo de seu  \t|")
    print("|\tcomprimento energizado na soldagem subaquática com arame \t|")
    print("|\ttubular.                                                 \t|")
    print("-------------------------------------------------------------------------\n")

    dict_parametros = carrega_parametros()

    parm = [float(dict_parametros[k]) for k in dict_parametros.keys() if "coef" in k]
    parm.insert(1, dict_parametros["dic_h"])
    parm = tuple(parm)

    cond_cont = (
        [
            float(dict_parametros["condicoes de contorno"]\
                  ["extremidade consumida"][0]),
            float(dict_parametros["condicoes de contorno"]\
                  ["extremidade consumida"][1])
        ],
        [
            float(dict_parametros["condicoes de contorno"]\
                  ["bico de contato"][0]),
            float(dict_parametros["condicoes de contorno"]\
                  ["bico de contato"][1])
        ]
    )

    interv = (
    float(
        dict_parametros["intervalo de integração"][0]
    ),
    float(
        dict_parametros["intervalo de integração"][1]
    ),
    float(
        dict_parametros["discretização do intervalo"]
    )
    )

    sol = edo_marcelo(parm, cond_cont, interv)

    df = pd.DataFrame(
        {
            "z" : sol.t,
            "t" : sol.y[0],
            "dtdz" : sol.y[1]
        }
    )
    df.to_csv("dados.csv", sep=',')

if __name__ == "__main__":
    main()

