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
	if './bin/j' == l[:7]:
		kerneln=t[0]

		if "spmm" in kerneln:
			kernel = " SpMM"
		else:
			kernel = " SDDMM"
			
		if "model" in kerneln:
			kernel += " model"
		else:
			kernel += ""
#		if "spmm" in kerneln:
			
		k = t[-2] #.split('_')[1][9:]
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
	if t[0] == 'J':
		k = t[9]
		ti = t[10][3:-1]#""
		tk = t[12][:-1]#""
		kernel = t[2]
		print("J_Stream_" + (kernel)  +","+matrix+','+k+','+ t[3]+ ',' + t[-1] + "," +ti+','+tk)
		func_i = (func_i + 1) % 3


