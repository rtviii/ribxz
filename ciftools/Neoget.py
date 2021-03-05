from neo4j import GraphDatabase, Result
from dotenv import load_dotenv
import os

if (__name__ =='__main__'):
    print("This file is not supposed to be run as a script")

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
# Neo4j driver has to be created inside the __init__ of the object to persist with its life. Otherwise is consumed on the first call.
def _neoget(CYPHER_STRING:str)->Result:
    driver = GraphDatabase.driver(
    os.getenv( 'NEO4J_URI' ),
    auth=(os.getenv( 'NEO4J_USER' ),
    os.getenv( 'NEO4J_PASSWORD' )))

    def parametrized_query(tx, **kwargs):
        result:Result = tx.run(CYPHER_STRING, **kwargs)
        return result.values()

    with driver.session() as session:
        session.close()
        return session.read_transaction(parametrized_query)
    