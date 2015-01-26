"""
*** Discrete Sampling Practice ***

---> Creating useful functions <---

(1) Write a function that multiplies all consecutively decreasing numbers between a maximum and a minimum supplied as arguments. (Like a factorial, but not necessarily going all the way to 1). This calculation would look like

max * max-1 * max-2 * ... * min

(2) Using the function you wrote in (1), write a function that calculates the binomial coefficient (see Definition 1.4.12 in the probability reading). Actually, do this twice. The first time (2a) calculate all factorials fully. Now re-write the function and cancel as many terms as possible so you can avoid unnecessary multiplication (see the middle expression in Theorem 1.4.13).

(3) Try calculating different binomial coefficients using both the functions from (2a) and (2b) for different values of n and k. Try some really big values there is a noticeable difference in speed between the (2a) and (2b) function. Which one is faster? By roughly how much?

(4) Use either function (2a) or (2b) to write a function that calculates the probability of k successes in n Bernoulli trials with probability p. This is called the Binomial(n,p) distribution. See Theorem 3.3.5 for the necessary equation. [Hint: pow(x,y) returns x^y (x raised to the power of y).]

(5) Now write a function to sample from an arbitrary discrete distribution. This function should take two arguments. The first is a list of arbitrarily labeled events and the second is a list of probabilities associated with these events. Obviously, these two lists should be the same length.

---> Sampling sites from an alignment <---

Imagine that you have a multiple sequence alignment with two kinds of sites. One type of site pattern supports the monophyly of taxon A and taxon B. The second type supports the monophyly of taxon A and taxon C.

(6) For an alignment of 400 sites, with 200 sites of type 1 and 200 of type 2, sample a new alignment (a new set of site pattern counts) with replacement from the original using your function from (5). Print out the counts of the two types.

(7) Repeat (6) 100 times and store the results in a list.

(8) Of those 100 trials, summarize how often you saw particular proportions of type 1 vs. type 2. 

(9) Calculate the probabilities of the proportions you saw in (8) using the binomial probability mass function (PMF) from (4).

(10) Compare your results from (8) and (9).

(11) Repeat 7-10, but use 10,000 trials.
"""
from __future__ import division
import numpy
from scipy import stats
import matplotlib.pyplot as plt

def multiplier(num1, num2):
	""" Multiplies all consecutive numbers between num1 and num2 """
	start = max(num1, num2) #Get larger number
	stop = min(num1, num2) #Get smaller number
	product = 1 #Initialize product as 1
	while int(start) != int(stop)-1: #Loop while larger number is not equal to one less than smaller number
		product *= start #Multiply thre results
		start -= 1 #Decrease larger number by 1
	return product #Return product
	
def bincoeffa(n, k):
	""" Calculate binomial coefficient (multiplying first) from n sets choosing k subsets""" 
	try:
		coeff = multiplier(n, 1)/(multiplier(n-k, 1)*multiplier(k, 1)) 
	except:
		coeff = 0
	return coeff
	
def bincoeffb(n, k):
	""" Calculate binomial coefficient (factoring first) from n sets choosing k subsets"""
	try:
		coeff = multiplier(n, n-k+1)/multiplier(k, 1)
	except:
		coeff = 0
	return coeff

"""	
# To test speeds
print bincoeffa(100000, 200) #Roughly 14 secs
print bincoeffb(100000, 200) #Roughly <1 sec
"""

def bernoulli(k, n, p):
	"""Calculate a Probability Mass Function (PMF) of a binomial distribution with k successes in n Bernoulli trials with probability p"""
	PMF = bincoeffb(n, k)*pow(p,k)*pow((1-p),(n-k))
	return PMF

"""
#Plot a PMF distribution	
list = []	
for x in range(0,11):
	list.append(bernoulli(x, 10, 1/8))
	
plt.bar([y for y in range(0,11)], list)
plt.xlabel("k (Successes)")
plt.ylabel("PMF")
plt.xlim(0,10.0)
plt.show()	
"""

def sample(list, prob, outputsize):
	"""Resamples with replacement from original discrete distribution to give newlist"""
	new = stats.rv_discrete(name ="new", values=(list, prob))#Scipy function to sample from discrete distribution
	return new.rvs(size = outputsize)
		
def discsampl(sampsize):
	list = (1, 2) #List of events
	prob = (0.5, 0.5) #Probability of events
	val1 = [] #Empty list to hold occurrence of event 1 after sampling
	val2 = [] #Empty list to hold occurrence of event 2 after sampling
	for i in range(sampsize):
		val = sample(list, prob, 400)	#Call function
		value = []
		for all in val:
			value.append(all)
		val1.append(value.count(1)) #Count instances of event 1
		val2.append(value.count(2)) #Count instances of event 2

	"""print "Replicate\tEvent1\tEvent2"
	prop = [] #New list to store proportions of event 1 vs event 2
	for all in range(sampsize):
		#print all+1, "\t\t", val1[all], "\t", val2[all] #Prints all values
		prop.append(val1[all]/val2[all]) #Find proportions of all values
	"""	
	propfreq = {} #New dictionary to calculate frequencies
	for all in val1:
		if all in propfreq:
			propfreq[all]+=1
		else:
			propfreq[all]=1

	proportions = sorted(propfreq.keys())

	print "Proportion of events\t\tNo. of events\t\tProbability (PMF)"
	freqlist = []
	pmffreqlist = []
	for value in proportions:
		freqlist.append(propfreq[value]/sampsize)
		pmffreqlist.append(bernoulli(value, 400, 0.5))
		print "{0} of event 1:{1} of event 2".format(value, 400-value), "\t\t", propfreq[value]/sampsize,"\t\t", round(bernoulli(value, 400, 0.5),5) #Get Binomial PMF for observed occurrences in 400 sites with probability 0.5
	a = raw_input("Press enter to graph data")
	plt.plot(proportions, freqlist, 'ro', label="Observed (Sampled)")
	plt.plot(proportions, pmffreqlist, label="Expected (PMF)")
	plt.xlabel("Number of type 1 site")
	plt.ylabel("Probability")
	plt.legend()
	plt.show()
discsampl(100)
a = raw_input("Press enter to repeat with 10000 trials")
discsampl(1000)
