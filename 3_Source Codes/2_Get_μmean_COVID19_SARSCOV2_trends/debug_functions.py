import numpy as np
#import cupy as cp

def prv( name, arg ):
    print( f'\n{name} -> {type(arg)} -> {arg}' )


def pri( name, arg ):
    print( f'\n{name} -> {len(arg)} -> {arg}' )


def npri( name, array, axis=0 ):
    print( f'\n {name}')
    print( array, array.shape, array.size, np.count_nonzero( array, axis=axis ) )


def cpri( name, array, axis=0 ):
    print( f'\n {name}')
    print( array, array.shape, array.size,
           #cp.count_nonzero( array, axis=axis )
           )
