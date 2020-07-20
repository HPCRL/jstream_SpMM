import sys
import statistics as stat

f = open(sys.argv[1]).readlines()


matrix = 'matrix'
papi = ''
k='k'

res = []
for l in f:
	t = l.split()
	if 'Matrix' == l[:6]:
		#k = t[0].split('_')[1][1:]
		matrix = t[1].split('/')[-1].split(".mtx")[0].split(".asym")[0]
	if "K=" in l:
		k = t[4]

	if t[0] == 'CSB':
		res += [float(t[-1])]
	if t[0][:5] == "event":
		#print (t)
		papi += t[1]+t[3]+t[5]
		#print (papi)
	if len(t)>3 and t[2] == "mflop/sec:" and t[0] == "Flushed":
		#print(matrix+','+k+","+papi + str(float(t[3])/1000))
		print("CSB,"+matrix+','+k+"median," + str(stat.median(res)))
		print("CSB,"+matrix+','+k+"max," + str(max(res)))
		print("CSB,"+matrix+','+k+"min," + str(min(res)))
		papi = ''
		res = []
#	if "/mtx/" in l:		matrix = t[1].split('/')[-1].split(".mtx")[0].split(".asym")[0]
	if 'mflop' in l:
		papi = ''
