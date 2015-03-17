from __future__ import (division, print_function)
import random
import numpy
import math
from scipy.linalg import expm

class SequenceEvo(object):
	#Initialize object matrix
	def __init__(self, 	sequence = ["A", "C", "G", "T"], 
						matrix = 	[[-1, 1/3, 1/3, 1/3], 
									[1/3, -1, 1/3, 1/3], 
									[1/3, 1/3, -1, 1/3], 
									[1/3, 1/3, 1/3, -1]],
						time = 1):
		"""Takes a sequence list and the corresponding matrix for the model"""
		self.sequence = sequence #Set values
		self.matrix = matrix #Set transition matrix
		self.time = time #Set length of interval
		
		#Calculate probability matrix for rate change
		values = {} 
		diag = {}
		for base in sequence:
			row = matrix[sequence.index(base)]
			diag[base] = -1*min(row)
			values[base] = [x/(min(row)*-1) for x in row]
		self.values = values
		self.diag = diag
		runQT = self.calcMargProb()
		self.statfreq = runQT[0]

def sample(list, prob, outputsize = 1):
	"""Resamples with replacement from original discrete distribution to give newlist"""
	output = []
	for val in range(outputsize):
		cumfreq = []
		y = 0
		for i in range(len(prob)): #Calculate cumulative frequency list
			cumfreq.append(y + prob[i])
			y = cumfreq[i]
		rand = random.random()
		for i in cumfreq: #Discrete sampling
			if rand <= i:
				output.append(list[cumfreq.index(i)])
				break
	return output
		
class ContMarkov(SequenceEvo):
	#Continuous Markov chain
	def simulation(self, seq=None, seql = 50):
		"""Runs a simulation with a continuous Markov chain using parameters given"""
		sequence = ""
		if seq is None:
			seq = sample(self.sequence, self.statfreq, outputsize = seql)
			seq = ''.join(seq)
		else:
			seq = seq[0]
		for i in range(len(seq)):
			startchoice = seq[i]
			time_passed, time_elapsed = 0, 0
			states, timevals = [], []
			while time_passed <= self.time: #Does sampling till time elapsed is greater than total time specified
				states.append(startchoice)
				timevals.append(time_elapsed)			
				time_elapsed = (-1/self.diag[startchoice])*math.log(random.random()) #Calculate length of time after next instance
				listvals = list(filter(lambda x:x != startchoice, self.sequence)) #Get list of variables it can change to
				prob = list(filter(lambda x:x != -1, self.values[startchoice])) #Get list of probabilities for these variables
				nextval = sample(listvals, prob) #Samples new distribution
				time_passed += time_elapsed #Calculate time elapsed
				startchoice = nextval[0] #Reset the starting variable
			sequence += states[-1]
		return sequence
		
	def calcMargProb(self, time=100):
		#Calculate Marginal Probabilities of each outcome
		return (expm((numpy.array(self.matrix) * time)))
