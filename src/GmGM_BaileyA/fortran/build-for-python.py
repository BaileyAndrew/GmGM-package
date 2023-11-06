from numpy import f2py
import os

# For every .f90 file in this directory
for file in os.listdir('.'):
    if file.endswith('.f90'):
        # Compile it to a .so file
        module = file.split('.')[0]
        with open(file, 'r') as f:
            code = f.read()
            f2py.compile(code, modulename=module, extra_args = '--quiet', extension='.f90')