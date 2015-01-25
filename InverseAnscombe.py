'''
Created on Oct 7, 2014

@author: haiderriaz
'''
import math
from decimal import *
#from jarray import array , zeros
from array import * 
    

def loadOVSTable():    
    row1 = []
    row2 = []
    with open('Table.dat','r') as f:
        for line in f:
            r1, r2 = line.split(",")
            row1.append(float(r1.strip()))
            row2.append(float(r2.strip())) 
            
    return (row1 , row2)           


# x and y are the keypoints 
# xx are the values you want to use to interpolate
# m = y2 - y1 / x2 - x1
# m*(x2 - x1) = y2 - y1
# y2 = y1 + m*(x2 - x1) <-- x2 is the value you want to use to interpolate, y2 is the output
def interp1(x, y, xx):
  
    out = []
    for pt in xx:
        # Case #1 - if pt < x[0]
        if pt < x[0]:
            
            x1 = x[0]
            x2 = x[1]
            y1 = y[0]
            y2 = y[1]
        # Case #2 - if pt > x[end]
        elif pt > x[-1]:
            
            x1 = x[-2]
            x2 = x[-1]
            y1 = y[-2]
            y2 = y[-1]
        else:
            
            index = binary_search(x, pt) # Figures out where we need to access for the x points
            x1 = x[index]
            x2 = x[index+1]
            y1 = y[index]
            y2 = y[index+1]
    
        # Calculate the slope 
        m = (y2 - y1) / (x2 - x1)
    
        # Interpolate the point
        interp = y1 + m*(pt - x1)
    
        # Add the point to the list
        out.append(interp)
    
    return out


# Inputs - x is an array of values - values are sorted 
# search is the value we want to search for
def binary_search(x, search):
  
    # We are defining our range of searching
    # First is the left side of the range
    first = 0
    # Last is the right side of the range
    last = len(x) - 1
  
    # Calculate middle point
    
  
    # Until we find the number we are looking for
    # or if we don't find what we're looking for
    while (first <= last):
        
        middle = (first + last) / 2
        
    # Access the middle point and see if it's
    # less than what we're looking for 
        if x[middle] < search:
            first = middle + 1 # If the value in the middle is LESS than what we're looking for, move left range up
        elif x[middle] == search: # If the middle value is what we're looking for, get out of the loop
            break
        else: # If the value in the middle is GREATER than what we're looking for, move right range down
            last = middle -1
      
    # If we don't find what we're looking for, return the index of where it SHOULD fit
        if (first > last):
            middle = last
            break
  
        # Returns location of where search is found
      
    return middle







def InvAnscombe(fimg):
    
    (Efz , Ez) = loadOVSTable()

    # Estimation of the asympotic inverse AT (Classical one)
    asympinvfimg = map(lambda x: math.pow((x / 2.0) , 2.0) - (1.0/8.0), fimg) 
    
    
    # Estimation of the optimal inv AT (see paper [2]).
    optimalinvfimg = interp1(Efz , Ez , fimg)
    
    
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
        
    
    #print optimalinvfimg
    out = array('f' , optimalinvfimg )
    
    
    
    return out    
    
    
    
    