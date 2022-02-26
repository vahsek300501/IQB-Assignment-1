class LocalAlignment(object):
	"""Class to define sequence and function related to it"""
	def __init__(self, p_sequence1,p_sequence2,p_matchScore,p_missScore,p_gapScore):
		self.sequence1 = p_sequence1
		self.sequence2 = p_sequence2
		self.matchScore = p_matchScore
		self.missMatchScore = p_missScore
		self.gapScore = p_gapScore
		self.n = len(self.sequence1)
		self.m = len(self.sequence2)
		self.dpMatLocal = []
		self.localSequences = []

		for i in range(0,self.m+1):
			tmpRow = []
			for j in range(0,self.n+1):
				tmpRow.append(0)
			self.dpMatLocal.append(tmpRow)

	# Function to fill DP matrix in bottom up approach for local alignment
	def localAlignmentMatrix(self):
		for i in range(1,self.n+1):
			self.dpMatLocal[0][i] = 0
		for i in range(1,self.m+1):
			self.dpMatLocal[i][0] = 0

		for i in range(1,self.m+1):
			for j in range(1,self.n+1):
				gap1 = self.dpMatLocal[i-1][j] + self.gapScore
				gap2 = self.dpMatLocal[i][j-1] + self.gapScore
				alignScore = -999999999
				if self.sequence1[j-1] == self.sequence2[i-1]:
					alignScore = self.dpMatLocal[i-1][j-1] + self.matchScore
				else:
					alignScore = self.dpMatLocal[i-1][j-1] + self.missMatchScore
				self.dpMatLocal[i][j] = max(gap1,max(gap2,max(alignScore,0)))

	def getOptimalAlignment(self, seq1, seq2, i, j, ansSeq1, ansSeq2, cntScore):
		if i <= 0 or j <= 0:
			self.localSequences.append((ansSeq1,ansSeq2,cntScore))
			return
		if self.dpMatLocal[i][j] == 0:
			self.localSequences.append((ansSeq1,ansSeq2,cntScore))

		if self.dpMatLocal[i-1][j] + self.gapScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j, "_"+ansSeq1, seq2[i-1]+ansSeq2, max(cntScore + self.gapScore,0))

		if self.dpMatLocal[i][j-1] + self.gapScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i, j-1, seq1[j-1]+ansSeq1, "_"+ansSeq2, max(cntScore + self.gapScore,0))

		if seq1[j-1] == seq2[i-1] and self.dpMatLocal[i-1][j-1] + self.matchScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j-1, seq1[j-1]+ansSeq1, seq2[i-1]+ansSeq2, max(cntScore + self.matchScore,0))

		if seq1[j-1] != seq2[i-1] and self.dpMatLocal[i-1][j-1] + self.missMatchScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j-1, seq1[j-1]+ansSeq1,seq2[i-1]+ansSeq2, max(cntScore + self.missMatchScore,0))

	def optimalAlignment(self):
		self.getOptimalAlignment(self.sequence1,self.sequence2,self.m,self.n,"","",0)

	def printMat(self):
		for row in self.dpMatLocal:
			print(row)


mySeq = LocalAlignment("ATCAGAGTA","TTCAGTA",2,-1,-1)
mySeq.localAlignmentMatrix()
mySeq.printMat()
mySeq.optimalAlignment()
for row in mySeq.localSequences:
	print(row)