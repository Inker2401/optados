"""
Parser function parse() to parse the .odi file of OptaDOS.
"""
from __future__ import print_function

import inspect
import re
from collections import defaultdict

from . import show_output

core_chem_shift = re.compile("Mizoguchi chemical shift of E_TE\ \s*\:\s*([0-9\.-]+)\s*")

def parse(fname):
    """
    Open the file, parses it and return the values
    """
    retdict = defaultdict(list)

    if show_output:
        print("[{}.{}] Parsing file '{}'".format(
            __name__, inspect.currentframe().f_code.co_name, fname))

    with open(fname) as f:
        lines = f.readlines()

    for lno, l in enumerate(lines):

        # Return the value of the core_chemical_shift 
        match = core_chem_shift.search(l)
        if match:
            retdict["chem_shift"].append(float(match.groups()[0]))
            continue

        ###############################################################


    retdict = dict(retdict)
    if show_output:
        for k in sorted(retdict):
            print("  {}: {}".format(k, retdict[k]))
        print("-"*72)
    return retdict
