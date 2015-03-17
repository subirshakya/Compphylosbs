from __future__ import (division, print_function)
import ctmc

# ---> Defining Node and Tree classes <---
class Node:

	def __init__(self,name="",parent=None,children=None, branchlength = 0):
		self.name = name
		self.parent = None
		if children is None:
			self.children = []
		else:
			self.children = children
		self.brl = branchlength

class Tree:
	"""
	Defines a class of phylogenetic tree, consisting of linked Node objects.
	"""
	
	def __init__(self):
		"""
		The constructor really needs to be more flexible, but for now we're 
		going to define the whole tree structure by hand. This just uses
		the same statements we used above. By next Thurs (3/19), see if you can
		write a constructor that takes a parenthetical tree as its argument and 
		builds the corresponding tree in memory. 
		"""
		
		self.root = Node("root") 
		self.spC = Node("SpeciesC",parent=self.root)
		self.root.children.append(self.spC)
		self.ancAB = Node("ancAB",parent=self.root)
		self.root.children.append(self.ancAB)
		self.spA = Node("SpeciesA",parent=self.ancAB)
		self.spB = Node("SpeciesB",parent=self.ancAB)
		self.ancAB.children.append(self.spA)
		self.ancAB.children.append(self.spB)
		# Now, let's add branch lengths to our Node objects (remember, these fields
		# can be added arbitrarily in Python). In the future, we should probably include
		# branch lengths in the Node constructor.
		self.spA.brl = 0.1
		self.spB.brl = 0.1
		self.spC.brl = 0.2
		self.ancAB.brl = 0.1
		self.root.brl = 0
		# We're also going to add lists to each node that will hold simulated sequences.
		self.spA.seq = []
		self.spB.seq = []
		self.spC.seq = []
		self.ancAB.seq = []
		self.root.seq = []
		self.setModels(self.root)
		
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


	# Now, let's write a recursive function to simulate sequence evolution along a
	# tree. This amounts to simply simulating evolution along each branch 
	# from the root towards the tips. We'll need to use our ctmc class for setting the 
	# conditions of our simulation, which is why we imported it above our tree 
	# class definition. In this case, we've stored the definition of our ctmc 
	# class in a separate file (ctmc.py) to keep our tree code compact.
	# Now, let's add a ctmc object to each internal node in our tree (except the
	# root). Again, it would be best to add the ctmcs as part of the Node
	# constructor, if we know that we'll be simulating data.
	
	# Try to get this simulator and associated functions working by next Thurs. (3/19)    
	
	def setModels(self,node):
		"""
		This method of a Tree object defines a ctmc object associated with all
		nodes that have a branch length (i.e., all but the root).
		"""
		self.root.seq = [ctmc.ContMarkov().simulation()]
		self.spC.seq = [ctmc.ContMarkov(time = self.spC.brl).simulation(seq=self.root.seq)]
		self.ancAB.seq = [ctmc.ContMarkov(time = self.ancAB.brl).simulation(seq=self.root.seq)]
		self.spA.seq = [ctmc.ContMarkov(time = self.spA.brl).simulation(seq=self.ancAB.seq)]
		self.spB.seq = [ctmc.ContMarkov(time = self.spB.brl).simulation(seq=self.ancAB.seq)]

	def simulate(self,node):
		"""
		This method simulates evolution along the branches of a tree, taking
		the root node as its initial argument.
		"""
		pass
		
	def printSeqs(self):
		"""
		This method prints out the names of the tips and their associated
		sequences as an alignment (matrix).
		"""
		print (self.root.seq)
		print (self.spC.seq)
		print (self.ancAB.seq)
		print (self.spA.seq)
		print (self.spB.seq)

#data = "(Species C:0.2, (Species A:0.1, Species B: 0.1):0.1)"
#print (newicksplicer(data))

Sim = Tree()
Sim.printSeqs()
#Sim.printNames(Sim.root)
#print (Sim.treeLength(Sim.root))
#print (Sim.newick(Sim.root))
