'''
Created on Aug 18, 2014

@author: haiderriaz
'''

import array
from array import zeros
import math


def estimate(ima,x,y,z):
    
    minimum = reduce(min, ima)
    
    if (minimum < 0):
        ima = map(lambda x: x - minimum, ima)
     
    p1 =  math.pow(2 , math.ceil(math.log(x, 2)))
    p2 =  math.pow(2 , math.ceil(math.log(y, 2)))         
    p3 =  math.pow(2 , math.ceil(math.log(z, 2)))        
                
    # zeros pading 
    
    if p1 == x and p2 ==y and p3 == z:
        pad1 = array.array('f', ima)
        
    else:
        pad1 = zeros('f' , x*y*z)
        #for i in xrange(len(pad1)):
            
        
    
    