def cshift3D(x ,m , d):

	N1 = len(x)
	N2 = len(x[0])
	N3 = len(x[0][0])
	
	if d == 1:
		n = list(xrange(N1))
		n = (n-m) % N1
		
	elif d == 2:
		n = list(xrange(N2))
		n = (n-m) % N2

	elif d == 3:
		n = list(xrange(N3))
		n = (n-m) % N3

	return N1,N2,N3