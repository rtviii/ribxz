import sys
import argparse
import xml

def makeparser():
    
    def ResidueSpecification(string: str):
        try:
            chainId, seqnum = map(str, string.split(','))
            return chainId, seqnum
        except:
            print("Malformed or absent origin Residue specification")
            raise SyntaxError

    parser = argparse.ArgumentParser(
        "See MOLE2.pdf documentation for detailed options")

    #Pdbid
    parser.add_argument('-pdbid', '--pdb_id', dest='PDBID')

    # Paths
    parser.add_argument('-output', '--output_path', dest='Output',  help='output folder for mole')
    parser.add_argument('-input',  '--input_path',   dest='Input',   help='input file to operate on')
    # Params ----

    # Cavity
    parser.add_argument('-pr',  dest='ProbeRadius',         nargs='?', required=False)
    parser.add_argument('-it',  dest='InteriorThreshold',   nargs='?', required=False)
    parser.add_argument('-md',  dest='MinDepth',            nargs='?', required=False)

    # Tunnel
    parser.add_argument('-or',  dest='OriginRadius',        nargs='?', required=False)
    parser.add_argument('-scr', dest='SurfaceCoverRadius',  nargs='?', required=False)
    parser.add_argument('-br',  dest='BottleneckRadius',    nargs='?', required=False)
    parser.add_argument('-mtl', dest='MinTunnelLength',     nargs='?', required=False)
    parser.add_argument('-mpl', dest='MinPoreLength',       nargs='?', required=False)
    parser.add_argument('-mts', dest='MaxTunnelSimilarity', nargs='?', required=False)
    parser.add_argument('-bt',  dest='BottleneckTolerance', nargs='?', required=False)

    # Export
    parser.add_argument('--exports', dest='exports', choices=['t', 'tc', 'tcp'],
                        help="""
                            Specify exports to enable thus
                            [flag -> mole option]:
                            t  -> tunnels only
                            tc -> tunnesl and cavities
                            tcp -> tunnels, cavities, poresauto, poresuser, poresmerged
                            """)

    # Origins
    parser.add_argument('-o_points', dest='Points', nargs='*', action='append',
                        help='array of points to use as origin',
                        required=False)
    parser.add_argument('-o_residues', dest='Residues', nargs='*', action='append', type=ResidueSpecification,
                        help="""pass a comma-delimited string of Chain-name and Sequence number in that chain for each residue
    Ex. ...-o_residues A,130 A,145 A,160 ... results in <Residue Chain="A" SequenceNumber="130">\n
    <Residue Chain="A" SequenceNumber="145">\n 
    <Residue Chain="A" SequenceNumber="160">etc.
    """,
                        required=False)


    return parser
