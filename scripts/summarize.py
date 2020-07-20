import sys

accept_method='median'




def printdata(out,data,sddmm):


	out.write('Matrix,K,')
	for mtx in data:
		for k in data[mtx]:
			for kernel in sorted(data[mtx][k]):
				if ("SDDMM" in kernel ) == False:
					#print(kernel+',',end='')
					out.write(kernel+',')
			for kernel in sorted(data[mtx][k]):
				if ("SDDMM" in kernel ) == True:
					#print(kernel+',',end='')
					out.write(kernel+',')
			break
		break
	out.write('\n')
	#print()

	for mtx in data:
		
		for k in data[mtx]:
			out.write(mtx+',')
			out.write(k+',')
			for kernel in sorted(data[mtx][k]):
				if ("SDDMM" in kernel ) == False:
					#print(data[mtx][k][kernel]+',',end='')
					out.write(data[mtx][k][kernel]+',')
			for kernel in sorted(data[mtx][k]):
				if ("SDDMM" in kernel ) == True:
					out.write(data[mtx][k][kernel]+',')
			out.write('\n')
		#print()





f = open(sys.argv[1]).readlines()

data = {}
for l in f:
	t = l.split(',')
	method = t[0]
	mtx = t[1]
	k = t[2]
	fun = t[3]
	gflop = t[4].replace('\n','')
	if mtx not in data:
		data[mtx] = {}

	if k not in data[mtx]:
		data[mtx][k] = {}

	if fun == accept_method:
		data[mtx][k][method] = gflop


out_file=sys.argv[2]
f = open(out_file, 'w')

#printdata(f,data,False)

printdata(f,data,True)

#print(data)





