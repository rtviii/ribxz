from pymol import cmd

"""
A script to split 6OLE into two subunits and animate.

"""






cmd.fetch("6ole")

cmd.select("c. SL or c. SZ or c. Sf ")
cmd.select("not sele")
cmd.create("_6ole", "sele")
cmd.delete("6ole")

# Small subunit
cmd.select("""c. SW or c. SU or c. SG or c. SJ or c. SA or c. Se \
    or c. SV or c. SQ or c. SI or c. ST or c. Sc or c. SH or\
        c. Sa or c. SX or c. SB or c. SK or\
            c. SO or c. SN or c. SD or c. Sk or c. Sb\
            or c. SS or c. SF or c. Sz or c. Sd or c. SP \
or c. SC or c. SQ or c. SY or c. Sg or c. SE or c. SM or c. S2 or c. SR""")
cmd.create("SSU",'sele')

# Large subunit
cmd.select("_6ole and not sele")
cmd.create("LSU", "sele")

#Translate
cmd.translate("[0,-40,0]", "SSU")



