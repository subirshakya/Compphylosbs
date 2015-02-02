from __future__ import division
import random

def multiplier(num1, num2=1):
	""" Multiplies all consecutive numbers between num1 and num2 """
	start = max(num1, num2) #Get larger number
	stop = min(num1, num2) #Get smaller number
	product = stop #Initialize product as 1
	while int(start) != int(stop): #Loop while larger number is not equal to one less than smaller number
		product *= start #Multiply the results
		start -= 1 #Decrease larger number by 1
	return product #Return product
		
def bincoeffb(k, n):
	""" Calculate binomial coefficient (factoring first) from n sets choosing k subsets"""
	try:
		coeff = multiplier(n, n-k+1)/multiplier(k)
	except:
		coeff = 0
	return coeff
	
def bernoulli(k, n, p):
	"""Calculate a Probability Mass Function (PMF) of a binomial distribution with k successes in n Bernoulli trials with probability p"""
	PMF = bincoeffb(k, n)*pow(p,k)*pow((1-p),(n-k))
	return PMF

def likelihood (k, n, p):
	"""Calculate likelihood"""
	likelival = pow(p,k)*pow((1-p),(n-k))
	return likelival

def hillclimb(data, prob, likeli=0, diff = 0.1, oldvalue=1, direction="RIGHT"):
	"""Calculate maximum likelihood for a dataset"""
	while diff >= 0.001: #Terminator loop
		if likeli == oldvalue: #Terminate recursion if values equal and change difference parameter to start new recursive list
			likeli = likeli
			prob = prob
			diff = diff/5
			oldvalue = 1
			hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue, direction="RIGHT")
		else: #Recurvise function to look for maximum point
			oldvalue = likeli
			problist = [prob-diff, prob, prob+diff]
			valuelist = [likelihood(data[0], data[1], prob-diff), likelihood(data[0], data[1], prob), likelihood(data[0], data[1], prob+diff)]
			if direction == "LEFT":
				problist = problist[::-1]
				valuelist = valuelist[::-1]
			likeli = max(valuelist) #Choose maximum value of likelihood
			prob = problist[valuelist.index(likeli)]
			if prob<=0: #The program always goes left if two values are equal due to the max function. THis shifts it right.
				direction = "LEFT"
			hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue, direction=direction)
	return likeli, prob
	
"""Enter data below as pass and fail separated by commas"""
data = ["pass","pass","pass","pass","pass"]
data = [y.lower() for y in data]
dataval = data.count("pass")
numTrials = len(data)
#data = (dataval, numTrials)
data =(12,20)

if data[0] == data[1]: #Some special cases
	maxlikeli,  probability = 1.0, 1.0
elif data[0] == 0: # Special case 2
	maxlikeli, probability = 1.0, 0.0
else:
	prob = random.random() #To randomly choose starting probability
	maxlikeli, probability = hillclimb(data, 0.5)
print "The maximum likelihood for given data is {0:.4e} for probability of {1}".format(maxlikeli, round(probability, 5))
