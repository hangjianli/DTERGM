import numpy as np

def pretty_str(p, decimal_places=2, join_str='\n'):
    '''Pretty-print a matrix or vector.'''
    if type(p) == list:
        return join_str.join([pretty_str(x) for x in p])
    if len(p.shape) == 1:
        return vector_str(p, decimal_places)
    if len(p.shape) == 2:
        return matrix_str(p, decimal_places)
    if len(p.shape) == 3:
        return array_3d_str(p, decimal_places)
    raise Exception('Invalid array with shape {0}'.format(p.shape))
    
    
def vector_str(p, decimal_places=2):
    '''Pretty-print the vector values.'''
    style = '{0:.' + str(decimal_places) + 'f}'
    return '[{0}]'.format(", ".join([style.format(a) for a in p]))



def matrix_str(p, decimal_places=2):
    '''Pretty-print the matrix.'''
    return '[{0}]'.format("\n  ".join([vector_str(a, decimal_places) for a in p]))

def array_3d_str(p, decimal_places=2):
    '''Pretty-print the 3D array.'''
    return'[{0}]'.format("\n  ".join([matrix_str(a, decimal_places) for a in p]))


def estimate_changepts(params, cutoff=4):
    '''estimate the change points whose l2 norm of the change vector is larger than cutoff'''
    theta = params['theta']
    T = len(theta) - 1
    mask = np.linalg.norm(np.diff(theta, axis=0), ord=2, axis=1) > cutoff
    return np.arange(0, T)[mask]

