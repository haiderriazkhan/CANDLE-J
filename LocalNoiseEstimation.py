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
    