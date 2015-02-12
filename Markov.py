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
	output = [start]
	while steps!=1:
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
		
list = ("A", "B") #Provide list of names
InitProb = (0, 1) # Provide probabilities of items

#Define conditional matrix:
CondMat = [[0.50, 0.50], [0.50, 0.50]] #pA2A, pA2B, pB2A, pB2B
start = scipy.random.choice(list)
Markov(list, CondMat, 25, start)
#print Markov(list, CondMat, 25, start, timed="no")
