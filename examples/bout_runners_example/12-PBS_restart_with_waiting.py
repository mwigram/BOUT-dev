#!/usr/bin/env python

"""Driver which restarts a scan, given a restart function"""

import re
from bout_runners.bout_runners import PBS_runner

scan = [("cst", "D_perp", [1.0,5.5]),\
        ("cst", "D_par",  [1.5,2.5])
       ]

#{{{restart_from_func
def restart_from_func(dmp_folder,\
                      one_of_the_restart_paths_in_scan = None,\
                      scan_parameters = None,\
                      **kwargs):
    """
    Function which returns the restart from dmp_folder and
    one_of_the_restart_paths_in_scan

    NOTE: This will not work if the values of one of the scan parameters
          contains an underscore, or if the initial hit is in the root
          folder.

    Parameters
    ----------
    dmp_folder : str
        Given by the bout_runners. Used to find the current scan
        values.
    one_of_the_restart_paths_in_scan : str
        One of the restart paths from a previously run scan.
    scan_parameters : list
        List of strings of the names of the scan paramemters.
    kwargs : dict
        Dictionary with additional keyword arguments, given by
        bout_runners.

    Returns
    -------
    restart_from : str
        String which gives path to restart from
    """

    # Make a template string of one_of_the_restart_paths_in_scan
    restart_template = one_of_the_restart_paths_in_scan
    for scan_parameter in scan_parameters:
        hits = [m.start() for m in \
                re.finditer(scan_parameter, restart_template)]
        while(len(hits) > 0):
            # Replace the values with {}
            # The value is separated from the value by 1 character
            value_start = hits[0] + len(scan_parameter) + 1
            # Here we assume that the value is not separated by an
            # underscore
            value_len = len(restart_template[value_start:].split("_")[0])
            value_end = value_start + value_len
            # Replace the values with {}
            restart_template =\
                "{}{{0[{}]}}{}".format(\
                    restart_template[:value_start],\
                    scan_parameter,\
                    restart_template[value_end:])
            # Update hits
            hits.remove(hits[0])

    # Get the values from the current dmp_folder
    values = {}
    for scan_parameter in scan_parameters:
        hits = [m.start() for m in \
                re.finditer(scan_parameter, dmp_folder)]
        # Choose the first hit to get the value from (again we assume
        # that the value does not contain a _)
        value_start = hits[0] + len(scan_parameter) + 1
        # Here we assume that the value is not separated by an
        # underscore
        values[scan_parameter] = dmp_folder[value_start:].split("_")[0]

    # Insert the values
    restart_from = restart_template.format(values)

    return restart_from
#}}}

# Initial runs
# ===========================================================================
init_run = PBS_runner(additional = scan)

dmp_folder, PBS_ids =\
        init_run.execute_runs()

one_of_the_restart_paths_in_scan = dmp_folder[0]
# ===========================================================================


# Restart the scan
# ===========================================================================
restart_run = PBS_runner(nout         = 5                ,\
                         restart      = "overwrite"      ,\
                         restart_from = restart_from_func,\
                         additional   = scan             ,\
                        )

restart_run.execute_runs(\
                         # Declare dependencies
                         job_dependencies = PBS_ids,\
                         # Below are the kwargs given to the
                         # restart_from_func
                         one_of_the_restart_paths_in_scan =\
                         one_of_the_restart_paths_in_scan,\
                         scan_parameters = ["D_perp", "D_par"],\
                        )
# ===========================================================================
