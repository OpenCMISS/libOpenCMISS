"""
OpenCMISS Interface Generation
==============================

This python script parses opencmiss.F90 and generates
code for the bindings to other languages.

C
-

Integer constants are included in opencmiss.h

C functions are written to opencmiss_c.F90 from subroutines in opencmiss.F90
and function declarations are written to opencmiss.h.

For interfaces where one routine takes parameters as scalars and another takes
them as arrays, only the routine that takes arrays is included. This is done
by checking the routine names for "Number0", "Number1", "Number01" etc. so
relies on this naming convention to work correctly.

Extra arguments required: path to opencmiss.h, path to opencmiss_c.F90

Python
------

Generates a Python module that wraps the lower level C extension module
generated by SWIG.

SWIG
----

Generates a C header file with extra SWIG annotations to apply typemaps.

Extra arguments required: path to opencmiss.i

Testing
-------

The bindings generation scripts can be tested by running:
    python -m tests.all_tests

The tests are also useful for illustrating how the bindings generation works.

Limitations
-----------

- Doesn't support multi-dimensional arrays of CMISS types or Logicals, but
  this can be added if required.

- Doesn't account for the difference in storage order of multi-dimensional
  arrays between C and Fortran, except in CMISSC2FStrings.

"""

import os
import sys

from c import generate as c_generate
from python import generate as python_generate
from swig import generate as swig_generate


if len(sys.argv) >= 3:
    (cm_path, language) = sys.argv[1:3]
    extra_args = sys.argv[3:]
else:
    sys.stderr.write('Usage: %s opencmiss_src_path language language_specific_arguments\n'
            % sys.argv[0])
    exit(1)

languages = {'C': c_generate, 'Python': python_generate, 'SWIG': swig_generate}
if language not in languages.keys():
    sys.stderr.write('Language must be one of:\n')
    for l in languages:
        sys.stderr.write('  %s\n' % l)
    exit(1)

languages[language](cm_path, extra_args)
