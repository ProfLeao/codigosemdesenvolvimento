import numpy as np
import pandas as pd
from marcelo_solver.edo_marcelo import edo_marcelo
from marcelo_solver.diversos import carrega_parametros
import json

def main():
    dicpar = carrega_parametros()

    sol0, sol1 = edo_marcelo(
        (dicpar[k] for k in dicpar.keys() if "coef" in k), 
        (
            [
                dicpar["condicoes de contorno"]["extremidade consumida"][0],
                dicpar["condicoes de contorno"]["bico de contato"][0]
            ],
            [
                dicpar["condicoes de contorno"]["extremidade consumida"][1],
                dicpar["condicoes de contorno"]["bico de contato"][1]
            ]
        ),
        (
            dicpar["intervalo de integração"][0],
            dicpar["intervalo de integração"][1],
            dicpar["discretização do intervalo"]
        )
    )

    

