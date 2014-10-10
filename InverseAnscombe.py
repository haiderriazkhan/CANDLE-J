'''
Created on Oct 7, 2014

@author: haiderriaz
'''
import math

# Could this be made faster?
def loadOVSTable():
    with open('Table.dat','r') as f:
        rows=[map(float,L.strip().split(',')) for L in f]  # list of lists
    (row1 , row2) = zip(*rows)    
    return (row1 , row2)


def interp1(xx , yy):
    print "To do"



def InvAnscombe(fimg):
    
    (Efz , Ez) = loadOVSTable()

    # Estimation of the asympotic inverse AT (Classical one)
    asympinvfimg = math.pow((fimg / 2.0) , 2.0) - (1.0/8.0)
    
    # Estimation of the optimal inv AT (see paper [2]).
    # interp1 to do 
    optimalinvfimg = interp1()

    
    # Detection of image areas under and above the validity range of the optimal inverse transform
    maxEfz = max(Efz)
    asympmap = [i for i, e in enumerate(fimg) if e > maxEfz]
    OptInvConstant = 2.0* math.sqrt(3.0/8.0)
    biasedmap = [i for i, e in enumerate(fimg) if e < OptInvConstant]
    
    
    # Filling of areas under and above validity range
    for idx in asympmap:
        optimalinvfimg[idx] = asympinvfimg[idx]
        
    for idx in biasedmap:
        optimalinvfimg[idx] = 0
    
    
    
    