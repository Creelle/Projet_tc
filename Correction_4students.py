# test function

import numpy as np
import importlib
import sys

# insert at 1, 0 is the script path (or '' in REPL)
# adress = path'/GT.py'

GROUPS = np.arange(1) + 5

for gr in GROUPS:
    path = "Gr_%d" % gr;
    print(path)
    sys.path.insert(1, path)
    import GT_arguments as GT_arg;
    import GT as GT_student

    importlib.reload(GT_arg)
    importlib.reload(GT_student)

    ##Test function GT:
    # Create the input arguments
    inputs = GT_arg.GT_input();
    inputs.Pe = 50e3;

    # tests
    try:
        GToutputs = GT_arg.GT_outputs();
        GToutputs = GT_student.GT(inputs);
#
        # Print and save all figure generated by students
        for i in range(len(GToutputs.fig)):
            GToutputs.fig[i].show()
            figname = './%s/figures/GT%d.pdf' % (path, i)
            print(figname)
            GToutputs.fig[i].savefig(figname)
            GToutputs.fig[i].clf(True)

    except:
        print('\n Basic GT failed')

    # Now I remove the path:
    sys.path.remove(path)
    del inputs
    del GToutputs
