from requests.auth import  HTTPBasicAuth
import requests


USER = 'LDW_group'
PASS = 'RiboVision_Desire'


""" Check out the "Index" of the protein uL02:
>>> https://ribovision3.chemistry.gatech.edu/desire-api/alignments/?name=uL02
"""
""" According to that index, uL02 -- 1: 
>>> https://ribovision3.chemistry.gatech.edu/ortholog-aln-api/1/2157,2759,2
"""



def fetchAlignment():
    resp = requests.get('https://ribovision3.chemistry.gatech.edu/ortholog-aln-api/16/2157,2759,2', auth=HTTPBasicAuth(USER, PASS))
    if resp.status_code != 200:
        print(resp)
        raise NameError
    print(len( resp.json()['AA frequencies'] ))
    
    return
