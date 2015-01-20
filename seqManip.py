#Written for Python 2.7
#Windows 8

"""
With the skills you have now learned in the first part of this script, try this
exercise:
*** Sequence Manipulation Exercise ***
- Create a new Python script (text file)
- At the beginning of the script, define a DNA sequence (taken from
https://github.com/jembrown/CompPhylo_Spr2015/blob/master/CodingSeq.txt)
- Print the length of the sequence to the screen along with text explaining
the value
- Create and store the RNA equivalent of the sequence, then print to screen.
- Create and store the reverse complement of your sequence, then print to
screen.
- Extract the bases corresponding to the 13rd and 14th codons from the
sequence, then print them to the screen.
- Create a function to translate the nucleotide sequence to amino acids
using the vertebrate mitochondrial genetic code (available from
https://github.com/jembrown/CompPhylo_Spr2015/blob/master/VertMitTransTable.txt).
- Translate the sequence and print it to the screen.
- Be sure you've added comments to explain what this script is and what the
different bits of code mean.
- Save this script as "seqManip.py" and commit it to your class GitHub repo.
"""

def translate(sequence, codontable): #Function to translate sequence
	codonvals = codontable
	seq = sequence
	seqtrans = "" #Empty string
	for code in range(0,len(seq),3):
		codon = seq[code:code+3]
		if codon.upper() in codonvals: #Check to see if codon is in table
			seqtrans = seqtrans + codonvals[codon.upper()]
		else: #If codon shorter than 3 letters then continues
			continue
	return seqtrans
			

seq = "aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggctg"

print "The sequence is:\n", seq #print sequence
print

#LENGTH
print "The length of the sequence is:\n", len(seq) #print length of sequence
print

#RNA equivalent
RNAseq = seq.replace("t", "u") #Replace T's with U's
print "The RNA version of the string is:\n", RNAseq #print RNA version
print

#Reverse complement
compseq = "" #New list
for letter in range(len(seq)-1, -1, -1): 
	compseq = compseq + seq[letter] #Reverses sequence
compseq =  compseq.replace("a","G").replace("t","C").replace("c","T").replace("g","A").lower() #Complements sequence
print "The reverse complement of the sequence is:\n", compseq

#13th and 14th codon position
codonlist = [] #New list for codons
for codon in range(0,len(seq),3): #Go through sequences 3 at a time
	codonlist.append(seq[codon:codon+3]) #Add 3 points at a time
print
print "13th codon position:", codonlist[12]
print
print "14th codon position:", codonlist[13]
print

#Codon table
AAs = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG" #Copied from text
Starts = "--------------------------------MMMM---------------M------------"
Base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
Base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
Base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"

codontable = {} #Codon table list
for posn in range(len(Base1)): #Go through each position
	codontable[Base1[posn]+Base2[posn]+Base3[posn]]=AAs[posn] #Compile library with the 3 codon bases assigned to amino acid respectively

#TRANSLATE
translatedseq= translate(seq, codontable) #Call function
print "The translated sequence is:\n", translatedseq
