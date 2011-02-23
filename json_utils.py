import os
import simplejson

def json_serialize(obj, filename, use_jsonpickle=True):
    f = open(filename, 'w')
    if use_jsonpickle:
	import jsonpickle
	json_obj = jsonpickle.encode(obj)
	f.write(json_obj)
    else:
	simplejson.dump(obj, f, indent=1)	
    f.close()

def json_load_file(filename, use_jsonpickle=True):
    if os.path.isdir(filename):
	raise Exception, "%s is a directory -- expected a JSON filename." %(filename)
    f = open(filename)
    if use_jsonpickle:
	import jsonpickle
	import as_events	
	json_str = f.read()
	obj = jsonpickle.decode(json_str)
	return obj
    else:
	obj = simplejson.load(f)
    return obj
