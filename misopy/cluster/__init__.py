'''
misopy.cluster

Cluster factory class for different cluster types

@author: Aaron Kitzmiller
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
@contact: aaron_kitzmiller@harvard.edu
'''

def getClusterEngine(cluster_type):
    '''
    Returns the correct cluster engine
    '''
    classname = 'misopy.cluster.%s' % cluster_type
    try:
        parts = classname.split('.')
        module = ".".join(parts[:-1])
        m = __import__( module )
        for comp in parts[1:]:
            m = getattr(m, comp)            
        return m
    except ImportError as e:
        raise Exception('Unable to import %s: %s' % (classname,str(e)))
    return None
    