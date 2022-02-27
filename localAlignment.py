import numpy as np
import pandas as pd

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

	# Function to get optimal local alignment
	def getOptimalAlignment(self, seq1, seq2, i, j, ansSeq1, ansSeq2, cntScore):
		if self.dpMatLocal[i][j] == 0:
			self.localSequences.append((ansSeq1,ansSeq2,cntScore))
			return

		if self.dpMatLocal[i-1][j] + self.gapScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j, "_"+ansSeq1, seq2[i-1]+ansSeq2, max(cntScore + self.gapScore,0))

		if self.dpMatLocal[i][j-1] + self.gapScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i, j-1, seq1[j-1]+ansSeq1, "_"+ansSeq2, max(cntScore + self.gapScore,0))

		if seq1[j-1] == seq2[i-1] and self.dpMatLocal[i-1][j-1] + self.matchScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j-1, seq1[j-1]+ansSeq1, seq2[i-1]+ansSeq2, max(cntScore + self.matchScore,0))

		if seq1[j-1] != seq2[i-1] and self.dpMatLocal[i-1][j-1] + self.missMatchScore == self.dpMatLocal[i][j]:
			self.getOptimalAlignment(seq1, seq2, i-1, j-1, seq1[j-1]+ansSeq1,seq2[i-1]+ansSeq2, max(cntScore + self.missMatchScore,0))

	# function to call getOptimalAlignment method
	def optimalAlignment(self):
		for i in range(0,self.m+1):
			for j in range(0,self.n+1):
				if self.dpMatLocal[i][j] == self.dpMatLocal[self.m][self.n]:
					self.getOptimalAlignment(self.sequence1,self.sequence2,i,j,"","",0)

	def printMatrix(self):
		colIndex = [ch for ch in self.sequence1]
		rowIndex = [ch for ch in self.sequence2]
		rowIndex.insert(0," ")
		colIndex.insert(0," ")
		df = pd.DataFrame(np.array(self.dpMatLocal), columns = colIndex, index = rowIndex)
		print()
		print("Matrix Alignment")
		print(df)
		print()
		print()

	def printSequences(self):
		for seq in self.localSequences:
			for i in range(0,len(seq[0])):
				print(seq[0][i],end = " ")
			print()
			for i in range(0,len(seq[1])):
				print(seq[1][i],end = " ")
			print()
			print("Score Value: "+str(seq[2]))
			print()
			print()


def main():

	print("LOCAL SEQUENCE ALIGNMENT")
	print()
	print()

	print("Scoring Values-1")
	print("Match Score: 2   Miss Score: -1   Gap Score: -1")
	print()

	mySeq1 = LocalAlignment("ATCAGAGTA","TTCAGTA",2,-1,-1)
	mySeq1.localAlignmentMatrix()
	mySeq1.printMatrix()
	mySeq1.optimalAlignment()
	mySeq1.printSequences()
	print()

	print("Scoring Values-2")
	print("Match Score: 2   Miss Score: -1   Gap Score: -2")
	print()

	mySeq2 = LocalAlignment("ATCAGAGTA","TTCAGTA",2,-1,-2)
	mySeq2.localAlignmentMatrix()
	mySeq2.printMatrix()
	mySeq2.optimalAlignment()
	mySeq2.printSequences()

main()