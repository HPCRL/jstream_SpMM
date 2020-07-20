import sys

f = open(sys.argv[1]).readlines()


matrix = 'matrix'
papi = ''
k='128'
kernel='kernel'
func = ''
func_name=['median','max','min']
func_i=0
for l in f:
	t = l.split()
	#if './bin/mkl' == l[:9]:

		#kk = t[0].split('_')[1][9:]
		#print (t[0])
		#kernel = t[0].split('.')[1][10:]
		#print (kernel)
	if 'Matrix' == l[:6]:
		matrix = t[1].split('/')[-1].split(".mtx")[0].split(".asym")[0]
	if t[0][:5] == "event":
		#print (t)
		papi += t[1]+t[3]+t[5]
		#print (papi)
	#if len(t)>3 and t[2] == "mflop/sec:" and t[0] == "Flushed":
	if t[0] == 'MKL':
		k = t[-5]
		print("MKL,"+matrix+','+k+papi + func_name[func_i]+ ',' + t[-1] )
		func_i = (func_i + 1) % 3


