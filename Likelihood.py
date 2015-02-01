from __future__ import division
import scipy
import matplotlib.pyplot as plt

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
	

"""Enter data below as pass and fail separated by commas"""
data = ["pass","pass","pass","pass","pass"]
data = [y.lower() for y in data]
dataval = data.count("pass")
numTrials = len(data)

list = [x/100 for x in range(0,101,5)] #List with probabilities in range 0 to 1 spaced by 0.05

likeli = []

for val in list: #Calculate likelihoods
	num = bernoulli(dataval, numTrials, val)
	likeli.append(num)

maxval = max(likeli)
print "\nThe maximum likelihood is {0} for data with probability {1}\n".format(round(maxval,4), list[likeli.index(maxval)])

for numbers in range(0,len(list)): #Print data
	print "{0:0=3d}% balls black \t likelihood = {1:.4f} \t likelihood ratio = {2:.4f}".format(int(list[numbers]*100), round(likeli[numbers],4), round(likeli[numbers]/maxval, 4))


input = raw_input("Graph data (Y/N)") #Graph data Y/N
if input.upper() == "Y": 
	plt.plot(list, likeli)
	plt.xlabel("Probability")
	plt.ylabel("Likelihood")
	plt.show()
	
