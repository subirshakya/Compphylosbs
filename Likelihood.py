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

"""Enter data below as pass and fail separated by commas"""
data = ["pass","pass","pass","pass","fail"]
data = [y.lower() for y in data]
dataval = data.count("pass")
numTrials = len(data)

list = [x/100 for x in range(0,101,5)] #List with probabilities in range 0 to 1 spaced by 0.05

likelist = []

for val in list: #Calculate likelihoods
	num = bernoulli(4,5,val)
	likelist.append(num)

maxval = max(likelist)
print "\nThe maximum likelihood is {0:.4e} for data with probability {1}\n".format(maxval, list[likelist.index(maxval)])

for numbers in range(0,len(list)): #Print data
	print "{0:0=3d}% balls black \t likelihood = {1:.4e} \t ratio = {2:.4f}".format(int(list[numbers]*100), likelist[numbers], round(likelist[numbers]/maxval, 4))

input = raw_input("Graph data (Y/N)") #Graph data Y/N
if input.upper() == "Y": 
	plt.plot(list, likelist)
	plt.xlabel("Probability")
	plt.ylabel("Likelihood")
	plt.show()
	
