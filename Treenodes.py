from __future__ import (division, print_function)
import ctmc #Need to be in python folder Lib
from scipy.linalg import expm
import numpy

# ---> Defining Node and Tree classes <---
class Node:

	def __init__(self,name="",parent=None,children=None, branchlength = 0, sequence = None):
		self.name = name
		self.parent = None
		if children is None:
			self.children = []
		else:
			self.children = children
		self.brl = branchlength
		if sequence is None:
			self.sequence = []
		else:
			self.sequence = sequence
		self.likeli = [] #Saves marginal probs for the four bases

class Tree:
	"""
	Defines a class of phylogenetic tree, consisting of linked Node objects.
	"""
	
	def __init__(self, data, model = None, seqmat = None):
		"""
		The constructor really needs to be more flexible, but for now we're 
		going to define the whole tree structure by hand. This just uses
		the same statements we used above. By next Thurs (3/19), see if you can
		write a constructor that takes a parenthetical tree as its argument and 
		builds the corresponding tree in memory. 
		"""

		self.root = Node("root") #Define root
		self.newicksplicer(data, self.root) #Splice newick data
		if seqmat is None:
			if model is None:
				self.setModels(self.root)
			else:
				self.setModels(self.root, Model = model) #Default seqlength is 50, you can add any seqlength, default Model is JK, you can add model matrix
		else:
			seqdict = {} #Create dictionary of sequences if sequence provided not simulated
			for a in seqmat:
				seqdict[a[0].replace(" ", "")]=a[1]
			self.assignseqs(self.root, seqdict)
		
	def newicksplicer(self, data, base): 
		""" 
		Splices newick data to create a node based tree. Takes a base argument which is the root node
		"""
		
		data = data.replace(" ", "")[1: len(data)] 	 #Get rid of all spaces and removes first and last parenthesis
		n = 0
		if data.count(",") != 0: #While there is no more comma separated taxa
			for key in range(len(data)): #Find the corresponding comma for a given parenthesis (n will be 0 for the correct comma)
				if data[key] == "(":
					n += 1 #Increase index of n by 1 for 1 step into new node
				elif data[key] == ")":
					n -= 1 #Decrease index of n by 1 for 1 step outout node
				elif data[key] == ",": 
					if n == 0: #To check for correct comma
						vals = (data[0:key], data[key+1:len(data)-1]) #Break newick into left and right datasets
						for unit in vals: #For each entry of dataset
							if unit[-1] != ")": #For cases with branch lengths
								data = unit[0:unit.rfind(":")] #get rid of trailing branchlength if provided
								node_creater = Node(data, parent = base) #Create node entry
								node_creater.brl = float(unit[unit.rfind(":")+1:]) #Append branch length of that branch
								base.children.append(node_creater) #Create children
								self.newicksplicer(data, node_creater) #Recursive function
							else: #For case with no branch lengths
								node_creater = Node(data, parent = base)
								base.children.append(node_creater)
								self.newicksplicer(data, node_creater)
						break #Terminate loop, we don't need to look any further
	
	def printNames(self,node):
		"""
		A method of a Tree object that will print out the names of its
		terminal nodes.
		"""
		
		if node.children == []: #Identifies terminal node
			print (node.name) 
		else:
			for child in node.children:
				self.printNames(child)

	def treeLength(self,node):
		"""
		A method to calculate and return total tree length.
		"""
		
		tot_len = 0 
		if node.children == []: #Terminal branch returns branch length
			return node.brl
		else:
			tot_len += node.brl #Add length of internal branch
			for child in node.children:
				tot_len += self.treeLength(child) #Add length of terminal branch
			return tot_len
		
	def newick(self,node):
		"""
		A method of a Tree object that will print out the Tree as a 
		parenthetical string (Newick format).
		"""
		
		newick = "(" #Opening bracket
		if node.children == []: #Terminal branch returns name
			return node.name + ":" + str(node.brl)
		else:
			for child in node.children:
				if node.children[-1] == child: #Don't add commas to last entry
					newick += self.newick(child)
				else:
					newick += self.newick(child) + "," #Adds commas to non-last entries
			if node.brl != 0:
				newick += "):" + str(node.brl) #Adds closing bracket
			else:
				newick += ")"
			return newick    
	
	def setModels(self,node, Model = None, seqlength = None):
		"""
		This method of a Tree object defines a ctmc object associated with all
		nodes that have a branch length (i.e., all but the root).
		"""
		
		if Model is None and seqlength is None: #Set starting sequence for root and change model or sequence length if required
			node.sequence = [ctmc.ContMarkov().simulation()]
		elif Model is not None and seqlength is None:
			node.sequence = [ctmc.ContMarkov(matrix = Model).simulation()]
		elif Model is None and seqlength is not None:
			node.sequence = [ctmc.ContMarkov().simulation(seql=seqlength)]
		else:
			node.sequence = [ctmc.ContMarkov(matrix = Model).simulation(seql=seqlength)]
		self.simulate(node)

	def simulate(self,node):
		"""
		This method simulates evolution along the branches of a tree, taking
		the root node as its initial argument.
		"""
		
		if node.children == []:
			pass
		else:
			for child in node.children:
				child.sequence = [ctmc.ContMarkov(time = child.brl).simulation(seq=node.sequence)]
				self.simulate(child)
		
	def printSeqs(self, node):
		"""
		This method prints out the names of the tips and their associated
		sequences as an alignment (matrix).
		"""
		
		if node.children == []:
			print (node.name, "\t", node.sequence)
		for child in node.children:
			self.printSeqs(child)
		
	def assignseqs(self,node, dict):
		"""Assign given sequences to correct node of tree"""
		if node.name in dict:
			node.sequence = dict[node.name] #Save sequence
			node.likeli = self.statelistgen(node) #Save marginal probs
		else:
			for child in node.children:
				self.assignseqs(child, dict)
	
	def statelistgen(self, node):
		ambisites = 	{"A": [1,0,0,0], "C": [0,1,0,0], "G": [0,0,1,0], "T":[0,0,0,1],
						"M": [1,1,0,0], "R": [1,0,1,0], "W": [1,0,0,1], "S": [0,1,1,0], "Y": [0,1,0,1], "K": [0,0,1,1],
						"V": [1,1,1,0], "H": [1,1,0,1], "D":[1,0,1,1], "B": [0,1,1,1],
						"N": [1,1,1,1], "-":[1,1,1,1]}	
		output = []
		for i in range(len(node.sequence)):
			output.append(ambisites[node.sequence[i]])
		return numpy.array(output)#Returns a numpy array of the genetic code in terms of marginal probs
				
class Likelihood:
	"""Calculate likelihood of a tree given the data"""
	
	def __init__(self, matrix, node, likeli = 0):
		self.matrix = numpy.array(matrix)
		self.P = self.calcMargProb()
		self.generator(node) #Fills all nodes with respective conditional probs
		self.calclikeli(node) #Calculate likelihood
		
	def calcMargProb(self, time=100):
		#Calculate Marginal Probabilities of each outcome
		return (expm((numpy.array(self.matrix) * time)))	
		
	def generator(self, node):
		value = [] #To hold probabilities
		brl = [] #To hold branch lengths
		for child in node.children:
			if child.likeli == []:
				self.generator(child)
			value.append(child.likeli)
			brl.append(child.brl)
		node.likeli = self.calculator(value[0], value[1], brl[0], brl[1])
			
	def calculator(self, val1, val2, brl1, brl2):
		"""Calculates probability for node given probabilities for the descendent nodes"""
		matrix1 = self.calcMargProb(time=brl1).transpose() #Calculate matrix and transpose it to fit maths model for branch 1
		matrix2 = self.calcMargProb(time=brl2).transpose() #Calculate matrix and transpose it to fit maths model for branch 2
		get1 = [x.dot(matrix1) for x in val1] #Multiply Q matrix with all items in one base, and do it for all bases in sequence for branch 1
		get2 = [x.dot(matrix2) for x in val2] #Multiply Q matrix with all items in one base, and do it for all bases in sequence for branch 2
		masternext = [] 
		for i in range(len(get1)):
			nextcond = []
			for j in range(len(get1[i])):
				nextcond.append(get1[i][j] * get2[i][j]) #Just to set up the maths, prob for next function is the product of each row entry
			masternext.append(nextcond)
		return numpy.array(masternext) #Return new probability matrix for ancestral node
	
	def calclikeli(self, node):
		all = (self.P[0].dot(node.likeli.transpose())) #Multiply with marginal probabilities and add the terms to get total likelihood for each site
		likeli = 1
		for num in all:
			likeli *= num #Multiply likelihood of each site to give total likelihood for sequences
		print ("The likelihood of the given data is", likeli) #Print likelihood
	
data = "((Sp A:0.1, Sp B:0.13):0.13 , (Sp C: 0.15, (Sp D:0.1, Sp E:0.15):0.3):0.02)"
seqmat = 	[["Sp A", "NACA"],
			["Sp B", "NACC"],
			["Sp C", "NAGG"],
			["Sp D", "NATT"],
			["Sp E", "NATA"]]
matrix = 	[[-0.96, 0.38,0.39,0.19],
			[0.29,-0.78,0.39,0.1],
			[0.58,0.76,-1.44,0.1],
			[0.58,0.38,0.19,-1.15]]
Sim = Tree(data, model = matrix, seqmat=seqmat)
Likelihood(matrix, Sim.root)
#Sim.printSeqs(Sim.root)
#Sim.printNames(Sim.root)
#print (Sim.treeLength(Sim.root))
#print (Sim.newick(Sim.root))
