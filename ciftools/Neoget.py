from neo4j import GraphDatabase, BoltStatementResult
from dotenv import load_dotenv
import os

# Neo4j driver has to be created inside the __init__ of the object to persist with its life. Otherwise is consumed on the first call.

# SET THIS
envpath="./../.env" 
load_dotenv(envpath)
uri        =  os.getenv( 'NEO4J_URI' )
authglobal = (os.getenv( 'NEO4J_USER' ),
            os.getenv( 'NEO4J_PASSWORD' ))
current_db =  os.getenv( 'NEO4J_CURRENTDB' )

if (__name__ =='__main__'):
    print("This file is not supposed to be run as a script")

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
def Neoget(CYPHER:str)->BoltStatementResult:
    driver = GraphDatabase.driver(uri,auth=authglobal)
    def make_query(tx):
        return tx.run(CYPHER)

    with driver.session(database=current_db) as session:
        result:BoltStatementResult = session.read_transaction(make_query)
        driver.close()
    return result