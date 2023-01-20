

import json
def carrega_parametros(path = "parametros.json"):
    with open("parametros.json") as fparm:
        dic_param = json.load(fparm)
    
    return dic_param