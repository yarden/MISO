#!/usr/bin/env python
##
## Check if all the necessary modules to run MISO are available
##
import time
import os

def check_module_availability(required_modules):
    unavailable_mods = 0
    for module_name in required_modules:
	print "Checking for availability of: %s" %(module_name)
	try:
	    __import__(module_name)
            # Manually check for correct matplotlib version
            # required for sashimi_plot
            if module_name == "matplotlib":
                import matplotlib.pyplot as plt
                if not hasattr(plt, "subplot2grid"):
                    print "WARNING: subplot2grid function is not available in matplotlib. " \
                          "to use sashimi_plot, you must upgrade your matplotlib " \
                          "to version 1.1.0 or later. This function is *not* required " \
                          "for MISO use."
	except ImportError:
	    print "  - Module %s not available!" %(module_name)
	    unavailable_mods += 1
    if unavailable_mods != 0:
	print "Total of %d modules were not available. " \
              "Please install these and try again." %(unavailable_mods)
    else:
	print "All modules are available!"
    return unavailable_mods


if __name__ == '__main__':
#    from optparse import OptionParser
#    parser = OptionParser()
#    parser.add_option("--add-modules", dest="add_modules", action="store_true", default=False,
#		       help="Try to add modules using the 'modules' system.")
#    (options, args) = parser.parse_args()
    required_modules = ['numpy', 'scipy', 'simplejson', 'matplotlib',
                        'pysam']
 	
    check_module_availability(required_modules)
    
