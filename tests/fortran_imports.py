import sys
sys.path.append('./src/GmGM_BaileyA')
sys.path.append('./build')

import GmGM
import fortran_core

print(fortran_core.__doc__)
print(fortran_core.fortran_core.extract_d_values.__doc__)