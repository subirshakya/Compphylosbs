from __future__ import division
from __future__ import print_function
import random
from scipy.linalg import expm
from scipy.stats import expon
import numpy
import math

class SequenceEvo(object):
	#Initialize object matrix
	def __init__(self, sequence=["A","C","G","T"], 
				matrix=[[-1, 1/3, 1/3, 1/3], 
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
		
class Functions():
	#Functions used in the program
	def exppdf(self, lamb, x):
		"""Calculate pdf for x for an exponential function with rate lambda"""
		try:
			pdf = lamb * math.exp(-lamb * x)
		except:
			pdf = 0
		return pdf

	def hillclimb(self, data, prob, likeli=0, diff = 0.1, oldvalue=1, direction="RIGHT"):
		"""Calculate maximum likelihood for a dataset"""
		while diff >= 0.0001: #Terminator loop
			if likeli == oldvalue: #Terminate recursion if values equal and change difference parameter to start new recursive list
				likeli = likeli
				prob = prob
				diff = diff/10
				oldvalue = 1
				hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue, direction="RIGHT")
			else: #Recurvise function to look for maximum point
				oldvalue = likeli
				problist = [prob-diff, prob, prob+diff]
				valuelist = [bernoulli(data[0], data[1], prob-diff), bernoulli(data[0], data[1], prob), bernoulli(data[0], data[1], prob+diff)]
				if direction == "LEFT":
					problist = problist[::-1]
					valuelist = valuelist[::-1]
				likeli = max(valuelist) #Choose maximum value of likelihood
				prob = problist[valuelist.index(likeli)]
				if prob<=0: #The program always goes left if two values are equal due to the max function. THis shifts it right.
					direction = "LEFT"
				hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue, direction=direction)
		return likeli, prob
		
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
		
class ContMarkov(SequenceEvo):
	#Continuous Markov chain
	def simulation(self, start="random"):
		"""Runs a simulation with a continuous Markov chain using parameters given"""
		if start == "random": #Figure out start parameter as specific or random choice
			startchoice = Functions().sample(self.sequence, self.statfreq)
			startchoice = startchoice[0]
		else:
			startchoice = start
		time_passed, time_elapsed = 0, 0
		states, timevals = [], []
		while time_passed <= self.time: #Does sampling till time elapsed is greater than total time specified
			states.append(startchoice)
			timevals.append(time_elapsed)			
			time_elapsed = (-1/self.diag[startchoice])*math.log(random.random()) #Calculate length of time after next instance
			listvals = list(filter(lambda x:x != startchoice, self.sequence)) #Get list of variables it can change to
			prob = list(filter(lambda x:x != -1, self.values[startchoice])) #Get list of probabilities for these variables
			nextval = Functions().sample(listvals, prob) #Samples new distribution
			time_passed += time_elapsed #Calculate time elapsed
			startchoice = nextval[0] #Reset the starting variable
		return states, timevals
		
	def drawsim(self, states, timevals):
		#Visualize the state space
		multiplier = 40
		time_scaled = [int(x*multiplier/2) for x in timevals]
		for val in range(len(states)):
			print ("-"*time_scaled[val] + states[val], end="")
		print ("-"*int((self.time-sum(timevals))*multiplier/2))
		
	def endFreqs(self, trials, start="random"):
		#Calculate end frequencies after n trials
		endvals = []
		for x in range(trials):
			if start == "random":
				startchoice = Functions().sample(self.sequence, self.statfreq)
				startchoice = startchoice[0]
			else:
				startchoice = start
			output1, output2 = self.simulation(start=startchoice)
			endvals.append(output1[len(output1)-1]) #Get last value after simulation
		freqlist = {}
		for val in endvals: #Calculate frequencies of each value
			if val in freqlist:
				freqlist[val]+=1
			else:
				freqlist[val]=1
		frequencies = []
		for val in self.sequence:
			try:
				frequencies.append(freqlist[val]/trials)
			except:
				frequencies.append(0)
		return frequencies
		
	def calcHistProb(self, states, timevals):
		#Calculate historical probabilities from a given sequence 
		product = self.statfreq[self.sequence.index(states[0])] #Probability of 1st event
		for count in range(len(timevals)-1):# Probabilities of waiting times
			product *= Functions().exppdf(self.diag[states[count]], timevals[count+1])
		for count in range(len(states)-1):# Probabilities of change of states
			product *= self.values[states[count]][self.sequence.index(states[count+1])]
		product *= 1-(Functions().exppdf(self.diag[states[-1]], self.time-sum(timevals)))		#Probability of end waiting period
		return product
		
	def calcMargProb(self, time=100):
		#Calculate Marginal Probabilities of each outcome
		return (expm((numpy.array(self.matrix) * time)))
		
	def calcStateChange(self, start, end, time):
		#Calculate the probability of changing from start to end in given time 
		values = self.calcMargProb(time = time)
		return values[self.sequence.index(start)][self.sequence.index(end)]
		
#model = ContMarkov(sequence = ["A", "C", "G", "T"], matrix = [[-1.916, 0.541, 0.787, 0.588], [0.148, -1.069, 0.415, 0.506], [0.286, 0.170, -0.591, 0.135], [0.525, 0.236, 0.594, -1.355]], time = 2) #Non GTR model
model = ContMarkov(sequence = ["A","C","G","T"], matrix = [[-(1.87+4.25+2.53), 1.87, 4.25, 2.53],[1.87, -(1.87+0.62+8.7), 0.62, 8.7],[4.25, 0.62, -(4.25+0.62+1), 1], [2.53, 8.7, 1, -(2.53+8.7+1)]], time = 20) # GTR model
#states, times = model.simulation() 
#print (states) # Use this to print the states through time
#print (states, times) #Use this to print states and times
#model.drawsim(states, times) #Use to visualize states along timeline
#print (model.endFreqs(trials=1000)) #Use this to calculate frequency of an allele after certain generations
#print (model.calcMargProb()) #Use this to print the final stationary probabilities after given time 
#print (ContMarkov(time = 100).endFreqs(trials=100)) #Use the original matrix
#print (model.calcHistProb(["A","T","A","G"],[0, 0.09, 0.84, 0.15])) #Calculate probability of obtaining this series of data with the corresponding waiting times.
#print (model.calcStateChange(start="G", end="A", time=2)) #Calculate probability of getting end state from start state in given time
