# Use Matlab Engine from Python
import matlab.engine
eng = matlab.engine.start_matlab()
eng.alpha(eng.double(5), eng.double(6))



# Use Octave
# import os
# pathToExecutable = (
#     'C:/Octave/Octave-5.2.0/mingw64/bin/octave-cli.exe'
# )
# os.environ['OCTAVE_EXECUTABLE'] = pathToExecutable
# from oct2py import octave


# from oct2py import Oct2Py
# oc = Oct2Py()
# oc.alpha(5, 6)
