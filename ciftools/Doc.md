# Generating positional data:

<!-- Ligands -->
1. Alpha-shape of the whole structure
    - For each protein, each residues' distance to the alpha-shape

2. A ligand's neighboring residues within some radius, their parent chains

- Not every structur lists all ligands in non-polymers:

<!-- Protein -->
1. A protein's conservation across species, residue-wise

<!-- Protein Interfaces -->

<!-- Tunnel -->
1. Exit tunnel has to be generated(somehow), then:
2. For each protein, each residues' distance to the tunnel



# In the way of comments for now:
Assigining a negative value to the moletunnel of the dataframe to indicate that:
     - mole was unable to find the correct tunnel(0)
    - something is obviously blocking the tunnel (-1)
    - the ribosome is not canonically assembled(PTC is   fragmented/constriction site is far away from ptc)( -2)