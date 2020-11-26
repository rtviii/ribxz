from neo4j import GraphDatabase

uri = "bolt://ribosome.xyz:7687/"
driver = GraphDatabase.driver(uri, auth=("rt", "rrr"))

def myfunc(tx ):
    return tx.run("MATCH (n:RibosomeStructure) return n limit 1;")

with driver.session(database='ribban03') as session:
    result =session.read_transaction(myfunc);
    for x in result:
        print(x)

driver.close()
