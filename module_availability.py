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
	except ImportError:
	    print "  - Module %s not available!" %(module_name)
	    unavailable_mods += 1
    if unavailable_mods != 0:
	print "Total of %d modules were not available.  Please install these and try again." %(unavailable_mods)
    else:
	print "All modules are available!"
    return unavailable_mods

def add_modules(required_modules):
    """
    Try to add modules using the 'module' system.  Only supports Bash shell for now.
    """
    time_str = time.strftime("%m-%d-%y_%H:%M:%S")
    module_script_filename = "/tmp/add_module_%s.sh" %(time_str)
    module_script = open(module_script_filename, 'w')
    module_script.write("#!/bin/bash\n")
    for module_name in required_modules:
	print "Trying to add %s to path..." %(module_name)
	module_script.write("source /etc/profile.d/modules.sh\n")
	module_script.write("module add %s\n" %(module_name))
    module_script.close()
    os.system("chmod +x %s" %(module_script_filename))
    os.system("%s" %(module_script_filename))

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--add-modules", dest="add_modules", action="store_true", default=False,
		       help="Try to add modules using the 'modules' system.")
    (options, args) = parser.parse_args()
    
    required_modules = ['numpy', 'scipy', 'simplejson', 'jsonpickle', 'pygsl', 'matplotlib',
                        'pysam']
    
    if options.add_modules:
	add_modules(required_modules)
	
    check_module_availability(required_modules)
    
