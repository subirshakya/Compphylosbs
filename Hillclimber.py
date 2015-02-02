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
		
def bincoeffb(n, k):
	""" Calculate binomial coefficient (factoring first) from n sets choosing k subsets"""
	try:
		coeff = multiplier(n, n-k+1)/multiplier(k)
	except:
		coeff = 0
	return coeff
	
def bernoulli(k, n, p):
	"""Calculate a Probability Mass Function (PMF) of a binomial distribution with k successes in n Bernoulli trials with probability p"""
	PMF = bincoeffb(n, k)*pow(p,k)*pow((1-p),(n-k))
	return PMF

def hillclimb(data, prob, likeli=0, diff = 0.1, oldvalue=1):
	"""Calculate maximum likelihood for a dataset"""
	while diff >= 0.001: #Terminator loop
		if likeli == oldvalue: #Terminate recursion if values equal and change difference parameter to start new recursive list
			likeli = likeli
			prob = prob
			diff = diff/5
			oldvalue = 1
			hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue)
		else: #Recurvise function to look for maximum point
			oldvalue = likeli
			problist = [prob-diff, prob, prob+diff]
			valuelist = [bernoulli(data[0], data[1], prob-diff), bernoulli(data[0], data[1], prob), bernoulli(data[0], data[1], prob+diff)]
			likeli = max(valuelist) #Choose maximum value of likelihood
			if valuelist[0] == 0 and valuelist[1] == 0: #The program always goes left if two values are equal due to the max function. THis shifts it right.
				likeli = valuelist[2]
			prob = problist[valuelist.index(likeli)]
			hillclimb(data, prob=prob, likeli=likeli, diff=diff, oldvalue=oldvalue)
	return likeli, prob
	
"""Enter data below as pass and fail separated by commas"""
data = ["pass","pass","pass","pass","pass"]
data = [y.lower() for y in data]
dataval = data.count("pass")
numTrials = len(data)
data = (dataval, numTrials)

if data[0] == data[1]: #Some special cases
	maxlikeli,  probability = 1.0, 1.0
elif data[0] == 0: # Special case 2
	maxlikeli, probability = 1.0, 0.0
else:
	startlist = [y/10 for y in range(0,10,1)]
	berns = [bernoulli(data[0],data[1], y) for y in startlist] 
	prob = startlist[berns.index(random.choice(berns))] #To randomly choose starting probability
	maxlikeli, probability = hillclimb(data, prob)
print "The maximum likelihood for given data is {0} for probability of {1}".format(round(maxlikeli, 5), round(probability, 5))
