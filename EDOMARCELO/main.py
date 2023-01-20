import numpy as np
import pandas as pd
from marcelo_solver.edo_marcelo import edo_marcelo
import json

def main():
    with open("parametros.json") as fparm:
        parametros = json.load(fparm)

    y0 = 