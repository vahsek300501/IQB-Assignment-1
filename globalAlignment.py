import pdb
class GlobalAlignment(object):
	"""Class to define sequence and function related to it"""
	def __init__(self, p_sequence1,p_sequence2,p_matchScore,p_missScore,p_gapScore):
		self.sequence1 = p_sequence1
		self.sequence2 = p_sequence2
		self.matchScore = p_matchScore
		self.missMatchScore = p_missScore
		self.gapScore = p_gapScore
		self.n = len(self.sequence1)
		self.m = len(self.sequence2)
		self.dpMatGlobal = []
		self.globalSequences = []

		for i in range(0,self.m+1):
			tmpRow = []
			for j in range(0,self.n+1):
				tmpRow.append(0)
			self.dpMatGlobal.append(tmpRow)

	def printMatrix(self):
		for row in self.dpMatGlobal:
			print(row)


	# Function to fill DP matrix in bottom up approach for global alignment
	def globalAlignmentMatrix(self):
		for i in range(1,self.n+1):
			self.dpMatGlobal[0][i] = self.dpMatGlobal[0][i-1] + self.gapScore
		for i in range(1,self.m+1):
			self.dpMatGlobal[i][0] = self.dpMatGlobal[i-1][0] + self.gapScore

		for i in range(1,self.m+1):
			for j in range(1,self.n+1):
				gap1 = self.dpMatGlobal[i-1][j] + self.gapScore
				gap2 = self.dpMatGlobal[i][j-1] + self.gapScore
				alignScore = -999999999
				if self.sequence1[j-1] == self.sequence2[i-1]:
					alignScore = self.dpMatGlobal[i-1][j-1] + self.matchScore
				else:
					alignScore = self.dpMatGlobal[i-1][j-1] + self.missMatchScore
				self.dpMatGlobal[i][j] = max(gap1,max(gap2,alignScore))

	# Function to get optimal Alignments with the optimal score
	def getOptimalAlignment(self, seq1, seq2, i, j, ansSeq1, ansSeq2, cntScore):
		if i <= 0 or j <= 0:
			self.globalSequences.append((ansSeq1,ansSeq2,cntScore))
			return

		if self.dpMatGlobal[i-1][j] + self.gapScore == self.dpMatGlobal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j, "_"+ansSeq1, seq2[i-1]+ansSeq2, cntScore + self.gapScore)

		if self.dpMatGlobal[i][j-1] + self.gapScore == self.dpMatGlobal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i, j-1, seq1[j-1]+ansSeq1, "_"+ansSeq2, cntScore + self.gapScore)

		if seq1[j-1] == seq2[i-1] and self.dpMatGlobal[i-1][j-1] + self.matchScore == self.dpMatGlobal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j-1, seq1[j-1]+ansSeq1, seq2[i-1]+ansSeq2,cntScore + self.matchScore)

		if seq1[j-1] != seq2[i-1] and self.dpMatGlobal[i-1][j-1] + self.missMatchScore == self.dpMatGlobal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j-1, seq1[j-1]+ansSeq1,seq2[i-1]+ansSeq2,cntScore + self.missMatchScore)

	# Function to call getOptimalAlignment method
	def optimalAlignment(self):
		self.getOptimalAlignment(self.sequence1,self.sequence2,self.m,self.n,"","",0)


	def printMat(self):
		for val in self.globalSequences:
			print(val)

mySeq = GlobalAlignment("ATCAGAGTA","TTCAGTA",2,-1,-1)
mySeq.globalAlignmentMatrix()
mySeq.printMatrix()
mySeq.optimalAlignment()
mySeq.printMat()
