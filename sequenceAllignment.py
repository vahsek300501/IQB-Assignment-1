import pdb
from prettytable import PrettyTable

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
		self.dpMatGlobal = []
		self.dpMatLocal = []
		self.globalSequences = []
		self.localSequences = []

		for i in range(0,self.m+1):
			tmpRow = []
			for j in range(0,self.n+1):
				tmpRow.append(0)
			self.dpMatGlobal.append(tmpRow)

		for i in range(0,self.m+1):
			tmpRow = []
			for j in range(0,self.n+1):
				tmpRow.append(0)
			self.dpMatLocal.append(tmpRow)

	# Function to print DP Matrix
	def printdpMatGlobal(self):
		tmpMat = self.dpMatGlobal.copy()
		t = PrettyTable()
		for row in tmpMat:
			t.add_row(row)
		print(t)

	def printdpMatLocal(self):
		for row in self.dpMatLocal:
			for val in row:
				print(val,end = " ")
			print()

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

	# Recursive Function to generate all the optimal global alignment sequences
	def generateAllPossibleGlobalAlignmentUtil(self,seq1,seq2,ansSeq1,ansSeq2,cntScore):
		if len(seq1) == 0 and len(seq2) == 0:
			self.globalSequences.append((ansSeq1,ansSeq2,cntScore))
			return

		if len(seq1) == 0:
			self.generateAllPossibleGlobalAlignmentUtil(seq1, seq2[1:], ansSeq1+"_", ansSeq2+seq2[0], cntScore+self.gapScore)
			return

		if len(seq2) == 0:
			self.generateAllPossibleGlobalAlignmentUtil(seq1[1:], seq2, ansSeq1 + seq1[0], ansSeq2+"_", cntScore + self.gapScore)
			return

		if seq1[0] == seq2[0]:
			self.generateAllPossibleGlobalAlignmentUtil(seq1[1:], seq2[1:], ansSeq1 + seq1[0], ansSeq2 + seq2[0], cntScore + self.matchScore)
		else:
			self.generateAllPossibleGlobalAlignmentUtil(seq1[1:], seq2[1:], ansSeq1 + seq1[0], ansSeq2 + seq2[0], cntScore + self.missMatchScore)
		
		self.generateAllPossibleGlobalAlignmentUtil(seq1[1:], seq2, ansSeq1 + seq1[0], ansSeq2+"_", cntScore + self.gapScore)

		self.generateAllPossibleGlobalAlignmentUtil(seq1, seq2[1:], ansSeq1+"_", ansSeq2+seq2[0], cntScore+self.gapScore)

	def generateAllGlobalAlignments(self):
		self.generateAllPossibleGlobalAlignmentUtil(self.sequence1,self.sequence2,"","",0)
		ans = []
		for val in self.globalSequences:
			if val[2] == self.dpMatGlobal[self.m][self.n]:
				ans.append(val)

		for val in ans:
			for i in range(0,len(val[0])):
				print(val[0][i]+" ",end = "")
			print()
			for i in range(0,len(val[1])):
				print(val[1][i]+" ",end = "")
			print()
			print()
def main():
	mySeq = Sequence("ATCAGAGTA","TTCAGTA",2,-1,-1)
	print()
	print("GLOBAL SEQUENCE ALIGNMENT")
	print()
	print("**************2-D Array for optimal Global Sequence Alignment************")	
	print()
	mySeq.globalAlignmentMatrix()
	mySeq.printdpMatGlobal()
	print()
	print()
	print("**************Optimal Sequence Alignments*********************************")
	print()
	mySeq.generateAllGlobalAlignments()

main()
