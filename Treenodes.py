from __future__ import (division, print_function)
import ctmc #Need to be in python folder Lib

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

class Tree:
	"""
	Defines a class of phylogenetic tree, consisting of linked Node objects.
	"""
	
	def __init__(self, data):
		"""
		The constructor really needs to be more flexible, but for now we're 
		going to define the whole tree structure by hand. This just uses
		the same statements we used above. By next Thurs (3/19), see if you can
		write a constructor that takes a parenthetical tree as its argument and 
		builds the corresponding tree in memory. 
		"""

		self.root = Node("root") #Define root
		self.newicksplicer(data, self.root) #Splice newick data
		self.setModels(self.root) #Default seqlength is 50, you can add any seqlength, default Model is JK, you can add model matrix
		
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
			node.sequence = [ctmc.ContMarkov(model = Model).simulation()]
		elif Model is None and seqlength is not None:
			node.sequence = [ctmc.ContMarkov().simulation(seql=seqlength)]
		else:
			node.sequence = [ctmc.ContMarkov(model = Model).simulation(seql=seqlength)]
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
		
data = "((Species A:1, Species B:1):2 , (Species C: 3, (Species D:4, Species E:4):5):6)"

Sim = Tree(data)
#Sim.printSeqs(Sim.root)
#Sim.printNames(Sim.root)
#print (Sim.treeLength(Sim.root))
#print (Sim.newick(Sim.root))
