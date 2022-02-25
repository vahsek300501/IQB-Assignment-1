import pdb
class Sequence(object):
	"""Class to define sequence and function related to it"""
	def __init__(self, p_sequence1,p_sequence2,p_matchScore,p_missScore,p_gapScore):
		self.sequence1 = p_sequence1
		self.sequence2 = p_sequence2
		self.matchScore = p_matchScore
		self.missMatchScore = p_missScore
		self.gapScore = p_gapScore
		self.n = len(self.sequence1)
		self.m = len(self.sequence2)
		self.dpMat = []

		for i in range(0,self.m+1):
			tmpRow = []
			for j in range(0,self.n+1):
				tmpRow.append(0)
			self.dpMat.append(tmpRow)

	# Function to print DP Matrix
	def printDPMatrix(self):
		for row in self.dpMat:
			print(row)

	# Function to fill DP matrix in bottom up approach
	def fillDPMatrix(self):
		for i in range(1,self.n+1):
			self.dpMat[0][i] = self.dpMat[0][i-1] + self.gapScore
		for i in range(1,self.m+1):
			self.dpMat[i][0] = self.dpMat[i-1][0] + self.gapScore


mySeq = Sequence("abcde","abc",1,-2,-1)
mySeq.fillDPMatrix()
mySeq.printDPMatrix()