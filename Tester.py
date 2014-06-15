try:
    import cshift3D
    import pprint
except ImportError:
    print "Something is wrong"		

from jarray import array
 
data = [[[1, 2, 3 ], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]], [[1, 2, 3], [4, 5, 6]]]
x = array(data, Class.forName('[[D'))

n = 3
distance = [[[0 for k in xrange(n)] for j in xrange(2)] for i in xrange(n)]
distance[0][0][0] = 1; distance[1][1][2] = 3
distance[0][1][2] = 4; distance[2][0][2] = 8

c = 3

l = [1,2,3]

l = l[c:] + l[:c]

print distance