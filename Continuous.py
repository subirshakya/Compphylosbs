from __future__ import division
from __future__ import print_function
import random
from scipy.linalg import expm
import numpy
import math

class SequenceEvo(object):
	#Initialize object matrix
	def __init__(self, sequence=["A","C","G","T"], 
				matrix=[[-1, 1/3, 1/3, 1/3], [1/3, -1, 1/3, 1/3], [1/3, 1/3, -1, 1/3], [1/3, 1/3, 1/3, -1]], 
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
		
class Markov(SequenceEvo):
	#Discrete Markov chain
	def sample(self, list, prob, outputsize = 1):
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
		
	def Markov(self, steps, start, timed = "yes"):
		"""Processes a Markov chain from a list to a defined number of steps"""
		output = []
		while steps!=0:
			problist = self.matrix[self.sequence.index(start)]
			next = sample(self.sequence, problist)
			start = next[0]
			steps-=1
			if timed == "yes":
				print (next[0],)
				time.sleep(1)
			else:
				output.append(next[0])
		return output
		
class ContMarkov(SequenceEvo):
	#Continuous Markov chain
	def simulation(self, start="random"):
		"""Runs a simulation with a continuous Markov chain using parameters given"""
		if start == "random": #Figure out start parameter as specific or random choice
			start = random.choice(self.sequence)
		else:
			start = start
		time_passed, time_elapsed = 0, 0
		states, timevals = [], []
		while time_passed <= self.time: #Does sampling till time elapsed is greater than total time specified
			states.append(start)
			timevals.append(time_elapsed)			
			time_elapsed = (-1/self.diag[start])*math.log(random.random()) #Calculate length of time after next instance
			listvals = list(filter(lambda x:x != start, self.sequence)) #Get list of variables it can change to
			prob = list(filter(lambda x:x != -1, self.values[start])) #Get list of probabilities for these variables
			nextval = Markov().sample(listvals, prob) #Samples new distribution
			time_passed += time_elapsed #Calculate time elapsed
			start = nextval[0] #Reset the starting variable
		return states, timevals
		
	def drawsim(self, states, timevals):
		#Visualize the state space
		multiplier = 40
		time_scaled = [int(x*multiplier/2) for x in timevals]
		for val in range(len(states)):
			print ("-"*time_scaled[val] + states[val], end="")
		print ("-"*int((self.time-sum(timevals))*multiplier/2))
		
	def endFreqs(self, start, trials):
		#Calculate end frequencies after n trials
		endvals = []
		for x in range(trials):
			output1, output2 = model.simulation(start=start)
			endvals.append(output1[len(output1)-1]) #Get last value after simulation
		freqlist = {}
		for val in endvals: #Calculate frequencies of each value
			if val in freqlist:
				freqlist[val]+=1
			else:
				freqlist[val]=1
		frequencies = []
		for val in self.sequence:
			frequencies.append(freqlist[val]/1000)
		return frequencies
		
	def calcHistProb(self, states, timevals):
		#Calculate historical probabilities
		pass
	
	def calcMargProb(self):
		#Calculate Marginal Probabilities of each outcome
		return (expm((numpy.array(self.matrix) * self.time)))
			
#statetime, endFreqs, calcHistProb, calcMargProb functions
#Use function scipy.linalg.expm(Qv)
		
model = ContMarkov(sequence = ["A", "C", "G", "T"], matrix = [[-1.916, 0.541, 0.787, 0.588], [0.148, -1.069, 0.415, 0.506], [0.286, 0.170, -0.591, 0.135], [0.525, 0.236, 0.594, -1.355]], time = 10)
#states, times = model.simulation() 
#print (states) # Use this to print the states through time
#print (states, times) #Use this to print states and times
#model.drawsim(states, times) #Use to visualize states along timeline
#print (model.endFreqs(start = "A", trials=1000)) #Use this to calculate frequency of an allele after certain generations
#print (model.calcMargProb()) #Use this to print the final stationary probabilities after given time 
