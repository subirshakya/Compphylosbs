from __future__ import division
import scipy
import time

#Discrete sampling
def sample(list, prob, outputsize = 1):
	"""Resamples with replacement from original discrete distribution to give newlist"""
	output = []
	for val in range(outputsize):
		cumfreq = []
		y = 0
		for i in range(len(prob)):
			cumfreq.append(y + prob[i])
			y = cumfreq[i]
		rand = scipy.random.random()
		for i in cumfreq:
			if rand <= i:
				output.append(list[cumfreq.index(i)])
				break
	return output

def Markov(list, condMat, steps, start, timed = "yes"):
	"""Processes a Markov chain from a list to a defined number of steps"""
	output = []
	while steps!=0:
		problist = condMat[list.index(start)]
		next = sample(list, problist)
		start = next[0]
		steps-=1
		if timed == "yes":
			print next[0],
			time.sleep(1)
		else:
			output.append(next[0])
	return output
	
		
list = ("A", "T", "C", "G") #Provide list of names
InitProb = (0.25, 0.25, 0.25, 0.25) # Provide probabilities of items

#Define conditional matrix:
#CondMat = [[0.25 for y in range(len(InitProb))] for y in range(len(InitProb))]
CondMat = [[0.5, 0.5,0,0],[0,0.5,0.5,0],[0,0,0.5,0.5],[0.5,0,0,0.5]]
start = scipy.random.choice(list)

for i in range(100):
	output = Markov(list, CondMat, 100, start, timed="no")
	print start, output[99]
"""
#mutator
seq = "ATTCATTAGAGTTCCATAT"
newseq = ""
for j in range(len(seq)):
	seqmut =  Markov(list, CondMat,100, seq[j], timed="no")
	newseq += seqmut[99]
print seq,"\n",newseq
"""
