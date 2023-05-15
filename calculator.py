"""
This is a code to calculate the similarities between protein sequences. It allows the user to input an MSA file in
phylip (PHY) or fasta (FA) format and choose the scoring matrix they want. The sequences are then compared using
the choosen scoring methodology and go through a neighbor joining and tree building process. The results are then
confirmed via boostrapping.

Several packages may need to be downloaded on computer based on your current computer packages.
One of these is the Bio package.
Download directions at: https://biopython.org/wiki/Download

Authors- Troy Hofstrand, Nicholas Sullivan, and Nikolaus Ryczek
Emails- troy.hofstrand@slu.edu, nicholas.sullivan@slu.edu, nikolaus.ryczek@slu.edu

Last Date Updated- 4/24/2023
"""

import os
import warnings
import pandas as pd
import itertools
import numpy as np
from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename
from Bio import SeqIO, Phylo, AlignIO
from Bio.Phylo.Consensus import majority_consensus
import re


warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", RuntimeWarning)

class Calculations:
	# A function to write a tree to a file as a Newick-format string 
	def write_tree(newick_path, tree, taxon_labels):
		newick_string = Calculations.make_newick_string(len(tree) - 1, tree, taxon_labels) + ";\n"

		newick_file = open(newick_path, "a")
		newick_file.write(newick_string)
		newick_file.close()
		Calculations.summaryFile()

	# Recursively build a Newick-format string from an adjacency list
	def make_newick_string(node_i, tree, taxon_labels):
		if len(tree[node_i]) == 0: # no outgoing edges, so must be a leaf node
			return taxon_labels[node_i]
		else: # an internal node
			newick_substrings = []
			for child_i in sorted(tree[node_i]):
				branch_length = tree[node_i][child_i]
				substring = Calculations.make_newick_string(child_i, tree, taxon_labels)
				newick_substrings.append("%s:%f" % (substring, branch_length))

			return "(" + ",".join(newick_substrings) + ")"

	# Read in a PHYLIP-format multiple sequence alignment
	def read_phylip(phylip_path):
		records = list(SeqIO.parse(phylip_path,'phylip'))
		alignment = AlignIO.read(phylip_path, "phylip")

		n_taxa = len(records)
		sequence_length = alignment.get_alignment_length()
		sequences = []
		sequence_labels = []
		msa = np.zeros((n_taxa, sequence_length), dtype = "uint8")

		for i in range(n_taxa):
			label = records[i].id
			sequence = str(records[i].seq)
			sequence_labels.append(label)
			sequences.append(sequence)
			msa[i] = np.fromstring(sequence, dtype = "uint8")

		return(sequences, sequence_length, sequence_labels, msa)

	# Read in a FASTA-format multiple sequence alignment
	def read_fasta(fasta_path):
		records = list(SeqIO.parse(fasta_path,'fasta'))
		alignment = AlignIO.read(fasta_path, "fasta")

		n_taxa = len(records)
		sequence_length = alignment.get_alignment_length()
		sequences = []
		sequence_labels = []
		msa = np.zeros((n_taxa, sequence_length), dtype = "uint8")

		for i in range(n_taxa):
			label = records[i].id
			sequence = str(records[i].seq)
			sequence_labels.append(label)
			sequences.append(sequence)
			msa[i] = np.fromstring(sequence, dtype = "uint8")

		return(sequences, sequence_length, sequence_labels, msa)

		'''
		KEY:
		65 - A, 82 - R, 78 - N, 68 - D, 67 - C, 81 - Q, 69 - E, 71 - G, 72 - H, 73 - I, 76 - L, 75 - K, 77 - M,
		70 - F, 80 - P, 83 - S, 84 - T, 87 - W, 89 - Y, 86 - V, 45 - -
		'''

	# Find the coordinates of whatever off-diagonal element of a matrix has the lowest value
	def get_lowest_off_diagonal_value_coordinate(matrix):
		lowest_value = None
		lowest_value_coordinate = None

		for i, row in enumerate(matrix):
			for j, value in enumerate(row):
				if i != j:
					if lowest_value == None or value < lowest_value:
						lowest_value = value
						lowest_value_coordinate = (i, j)

		return lowest_value_coordinate

	# Compute the neighbor joining Q-matrix from a distance matrix
	def compute_q_matrix(distance_matrix, matrix_length):
		q_matrix = (matrix_length - 2.0) * distance_matrix

		for i in range(matrix_length):
			for j in range(matrix_length):
				if i != j:
					q_matrix[i][j] -= np.sum(distance_matrix[i]) + np.sum(distance_matrix[j])

		return q_matrix

	def randComboCalc(msa_1,msa_2,combs=False):
		#REMOVE DOUBLE GAPS
		gap_ind = [i for i, (g, s) in enumerate(zip(msa_1, msa_2)) if g==s==45]
		seq1 = list(np.delete(msa_1,gap_ind))
		seq2 = list(np.delete(msa_2,gap_ind))

		#Get all combos of unique residues
		if combs==True:
			comb_list = list(itertools.product(np.unique(seq1), np.unique(seq2)))
			return comb_list,seq1,seq2
		else:
			return seq1,seq2
	
	#calculates the average identity score (percentage) for the entire msa 
	def msaIdentity(msa):
		seq = 0
		counter = 1
		pair_scores=[]
		while seq<len(msa)-1:
			score=0
			seq1 = []
			seq2 = []
			if counter == len(msa):
				seq += 1
				counter = seq+1
				if seq==len(msa)-1:
					break

			msa_1 = msa[seq]
			msa_2 = msa[counter]

			#double gaps removed and not equated into indentity score
			gap_ind = [i for i, (g, s) in enumerate(zip(msa_1, msa_2)) if g==s==45]
			seq1 = list(np.delete(msa_1,gap_ind))
			seq2 = list(np.delete(msa_2,gap_ind))

			length=len(seq1)
			for n in range(length):
				if seq1[n]==seq2[n]: #finds matching residuals and adds one to the score
					score+=1
			sim_score=score/length
			pair_scores.append(sim_score)
			counter += 1
		return sum(pair_scores)/len(pair_scores)
	
	#calls all scoring methods used in feng+doolittle equation to create distance matrix; used for all scoring methods other than pairwise
	def dist_mat(df,n_taxa,msa):
		real = Calculations.realscore(df,n_taxa,msa)
		print('real',real)
		rand = Calculations.randscore(df,n_taxa,msa)
		print('rand',rand)
		norm_scores = np.subtract(real , rand)
		identity = Calculations.identityscore(n_taxa,real)
		print('id',id)
		upper_norm = np.subtract(identity , rand)
		raw_dist = -np.log(np.divide(norm_scores, upper_norm))*100
		#c =  
		#distance_matrix = c * raw_dist
		print("Distance Matrix:\n", raw_dist)
		return raw_dist
	
	#function for pairwise scoring using feng+doolittle; returns distance matrix
	def pairwise(msa,n_taxa,gap):
		distance_matrix = np.zeros((n_taxa, n_taxa))
		seq = 0
		counter = 0

		#scoring matrix selection within scoring calculation since every pair uses own matrix; doesn't "matrix_selection()" like other scoring methods
		#scoring matrix selection is done in same manner as "msaIdentity()" but "msaIdentity()" processes entire MSA all at once
		while seq < n_taxa:
			score=0
			seq1 = []
			seq2 = []
			if counter == n_taxa:
				seq += 1
				counter = 0
				if seq == n_taxa:
					break

			msa_1 = msa[seq]
			msa_2 = msa[counter]

			gap_ind = [i for i, (g, s) in enumerate(zip(msa_1, msa_2)) if g==s==45]
			seq1 = list(np.delete(msa_1,gap_ind))
			seq2 = list(np.delete(msa_2,gap_ind))
			length=len(seq1)
			for n in range(length):
				if seq1[n]==seq2[n]:
					score+=1
			sim_score=score/length
			sim_score=round(sim_score,1)*100
			sim_score=round(sim_score) #round score to match scoring matrix names
			#adjust score in needed to match specific matrices we have available
			if sim_score<30:
				sim_score=30
			if sim_score==60:
				sim_score=62
			if sim_score>90:
				sim_score=90 
			print("sim:",sim_score)
			file = 'assets/matrices/BLOSUM%s.csv'%sim_score
			scoreMat = pd.read_csv(file, header=None).values
			# add penalty row and column of -2
			with_pen = np.empty([21,21])
			for i in range(len(scoreMat)):	
				with_pen[i] = np.append(scoreMat[i], [gap])
			with_pen[20] = (np.repeat(gap, 21))

			# make scoring dataframe
			labs = 'ARNDCQEGHILKMFPSTWYV-' #missing letters #BZX
			labels = np.fromstring(labs, dtype = "uint8") 
			df = pd.DataFrame(with_pen, columns = labels, index = labels)
			print(df)
			

			#distance calculation equation:feng+doolittle
			real=Calculations.realscoreP(df,seq,counter,msa)
			print('real:',real)
			rand=Calculations.randscoreP(df,seq,counter,msa)
			print('rand:',rand)
			rand=round(rand,1)
			identity=Calculations.identityscoreP(df,seq,counter,msa)
			print('id:',identity)
			norm_scores = real-rand
			upper_norm = identity-rand
			raw_dist = -np.log(norm_scores/upper_norm)*100
			raw_dist=round(raw_dist,1)
			print('raw',raw_dist)
			#add to complete distance matrix; adds one at a time unlike other scoring which does enitre matrix
			distance_matrix[seq][counter]=raw_dist
			counter += 1
		
		print(distance_matrix)
		return distance_matrix
	
	#real score in feng+doolittle for pairwise scoring
	def realscoreP(df,n,m,msa):
		real=0
		seq1,seq2 = Calculations.randComboCalc(msa[n],msa[m],combs=False)
		for j in range(len(seq1)):
			res1 = seq1[j]
			res2 = seq2[j]
			comp = df[res1][res2] #find in score mat
			real += comp 
		return real

	#random score in feng+doolittle for pairwise scoring
	def randscoreP(df,n,m,msa):
		rand_score=0
		seq1 = []
		seq2 = []

		msa_1 = msa[n]
		msa_2 = msa[m]

		comb_list,seq1,seq2 = Calculations.randComboCalc(msa_1,msa_2,combs=True)

		#Get scores for each combo
		gap_count1 = seq1.count(45) #get num of gaps seq 1
		gap_count2 = seq2.count(45)#get num of gaps seq 2
		for val in comb_list:
			comb_score = df[val[0]][val[1]] #find in score mat
			num1_occ = seq1.count(val[0]) #count num occ of first val (in first seq)
			num2_occ = seq2.count(val[1]) #count num occ of sec val (in sec seq)
			rand_score += (comb_score*num1_occ*num2_occ)  
		rand= (rand_score/len(seq1))- ((gap_count1+gap_count2)*2) #that score is spot in new_mat
		return rand

	#indentity score in feng+doolittle for pairwise scoring
	def identityscoreP(df,n,m,msa):
		seq1,seq2 = Calculations.randComboCalc(msa[n],msa[m],combs=False)
		score=0
		#compares each residue in each sequence to itself
		for i in seq1:
			score+=df[i][i]
		for i in seq2:
			score+=df[i][i]
		identity=score/2 #averages total score of two sequences
		return identity

	#real score in feng+doolittle
	def realscore(df,n_taxa,msa):
		#make distance matrix of zeros
		matrix_length = n_taxa
		distance_matrix = np.zeros((matrix_length, matrix_length))
		seq = 0
		counter = 0

		while True:
			if counter == n_taxa:
				print('seq',seq)
				seq += 1
				if seq == n_taxa:
					return distance_matrix
				counter = 0
			
			seq1,seq2 = Calculations.randComboCalc(msa[seq],msa[counter],combs=False)
			for j in range(len(seq1)):
				res1 = seq1[j]
				res2 = seq2[j]
				comp = df[res1][res2] #find in score mat
				distance_matrix[seq][counter] += comp #that score is spot in new_mat
			counter += 1


	# calculate random score between each sequence and add to random matrix
	def randscore(df,n_taxa,msa):
		distance_matrix = np.zeros((n_taxa, n_taxa))
		seq = 0
		counter = 0
		while True:
			rand_score=0
			seq1 = []
			seq2 = []
			if counter == n_taxa:
				print('seq',seq)
				seq += 1
				if seq == n_taxa:
					return distance_matrix
				counter = 0

			msa_1 = msa[seq]
			msa_2 = msa[counter]

			comb_list,seq1,seq2 = Calculations.randComboCalc(msa_1,msa_2,combs=True)

			#Get scores for each combo
			gap_count1 = seq1.count(45) #get num of gaps seq 1
			gap_count2 = seq2.count(45)#get num of gaps seq 2
			for val in comb_list:
				comb_score = df[val[0]][val[1]] #find in score mat
				num1_occ = seq1.count(val[0]) #count num occ of first val (in first seq)
				num2_occ = seq2.count(val[1]) #count num occ of sec val (in sec seq)
				rand_score += (comb_score*num1_occ*num2_occ)  
			distance_matrix[seq][counter] = (rand_score/len(seq1))- ((gap_count1+gap_count2)*2) #that score is spot in new_mat
			counter += 1

	# calculate identity score of sequence pairs by taking average of each identity score with itself
	def identityscore(n_taxa, reals):
		scores = reals
		matrix = np.zeros((n_taxa, n_taxa))
		i = 0
		while i < (n_taxa):
			for j in range(n_taxa):
				matrix[i, j] = (scores[i,i] + scores[j,j]) / 2 #gets score from real matrix of each sequence compared to itself
			i += 1
		return matrix

		

	# Read in the multiple sequence alignment
	def file_select():
		global fileName
		fileName=""
		print('Hello!\n'
		'Welcome to DNA similarity calculator!\n'
		'Select your desired phylip or fasta file:')
		Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
		filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
		print(filename + ' selected\n')
		Calculations.fileshort, file_type = os.path.splitext(filename)
		if file_type == '.txt':
			file_type = os.path.splitext(Calculations.fileshort)[1]
		if file_type == '.fasta' or file_type == '.fa' or file_type == '.fas':
			Calculations.sequences, Calculations.sequence_length, Calculations.taxon_labels, Calculations.msa = Calculations.read_fasta(filename)
		elif file_type == ".phy" or file_type == ".phylip":
			Calculations.sequences, Calculations.sequence_length, Calculations.taxon_labels, Calculations.msa = Calculations.read_phylip(filename)
		else:
			raise TypeError('Incorrect file type selected. Please choose a Fasta(.fa) or Phylip(.phy) file.')
			#tell app of error
		fileName=fileName+Calculations.fileshort
		print(fileName)

	#picks the matrix automatically or based off the user's selection
	def matrix_selection(value, gap,n_bootstrap = 1):
		global fileName
		global mat
		global fileNameF #required if application does multiple calculations without closing app
		Calculations.score = Calculations.msaIdentity(Calculations.msa)
		rounded = round(Calculations.score*10)*10 #rounded so every score is in the tens (10,20,30,etc.), easier to compare to matrices
		print('Rounded similarity score:', rounded)
		pair=False #required to determine which scoring technique is later chosen in "calculate_consensus_tree()"

		#auto assign scoring methods are unique cause program must do calculation to determine proper scoring matrix
		if value == 'Auto-assign BLOSUM based on average identity score' or value=='Auto-assign BLOSUM based on pairwise identity score':
			print("1\n")
			if value == 'Auto-assign BLOSUM based on pairwise identity score':
				fileNameF=fileName+'Pairwise' #adds unique 'Pairwise' to file name
				print(fileNameF)
				pair=True #set equal to true so proper scoring method is used
				print("1\n")

			#not elif because pairwise scoring still needs a selected matrix which is used in bootstrapping
			if rounded <= 30:
				arr = pd.read_csv('assets/matrices/BLOSUM30.csv', header=None).values
				mat="BLOSUM30"
				

			elif rounded == 60:
				arr = pd.read_csv('assets/matrices/BLOSUM62.csv', header=None).values
				mat="BLOSUM62"

			elif rounded >=90:
				arr = pd.read_csv('assets/matrices/BLOSUM90.csv', header=None).values
				mat="BLOSUM90"

			else:
				arr = pd.read_csv('assets/matrices/BLOSUM%s.csv'%rounded, header=None).values
				mat="BLOSUM%s"%rounded
				print("2\n")

			#adds matrix used to filename
			if pair==True:
				fileNameF=fileNameF+mat
			else:
				fileNameF=fileName+mat
		
		#if not auto assigned matrix is chosen prior by user
		else:
			file = 'assets/matrices/%s.csv'%value
			arr = pd.read_csv(file, header=None).values
			mat=value
			fileNameF=fileName+mat
			print(fileNameF)
			print("3\n")
		Calculations.calculate_consensus_tree(arr,gap, n_bootstrap,pair)
		
	
	def calculate_consensus_tree(score_mat,gap, n_bootstrap,pair):
		global fileNameF
		# add penalty row based on user selection
		gap=int(gap)
		with_pen = np.empty([21,21])
		for i in range(len(score_mat)):	
			with_pen[i] = np.append(score_mat[i], [gap])
		with_pen[20] = (np.repeat(gap, 21))

		# make scoring dataframe
		labs = 'ARNDCQEGHILKMFPSTWYV-' #missing letters #BZX
		labels = np.fromstring(labs, dtype = "uint8") 
		df = pd.DataFrame(with_pen, columns = labels, index = labels)

		# bootstrapping
		n_taxa = len(Calculations.taxon_labels)
		seq_length = Calculations.sequence_length
		cols = Calculations.msa.T
		newick_file = open(fileNameF+"Bootstraps.tre", "w")
		newick_file.close()

		#chose appropriate scoring method
		if pair==False:
			dist_mat = Calculations.dist_mat(df, n_taxa, Calculations.msa)
		else:
			dist_mat = Calculations.pairwise(Calculations.msa, n_taxa,gap )
		#first copy with original msa 
		Calculations.draw_tree(n_taxa, dist_mat)

		
		for i in range(n_bootstrap-1):
			#shuffle msa
			idx = np.random.randint(seq_length, size = seq_length)
			msa_new = cols[idx,:].T
			#calculate new distance matrix
			dist_mat = Calculations.dist_mat(df, n_taxa, msa_new)
			Calculations.draw_tree(n_taxa, dist_mat)
		
		trees = list(Phylo.parse(fileNameF+"Bootstraps.tre", "newick"))
		majority_tree = majority_consensus(trees)
		Phylo.write(majority_tree, fileNameF+"Consensus.tre", "newick")

		#get total distance of consensus tree and write it to summary file
		newick_string = open(fileNameF+"Consensus.tre", "r").read()
		distances = re.findall(r"(?<=:)[0-9]+(?:\.[0-9]+)?(?=[,);])", newick_string)
		total_dist = sum([float(dist) for dist in distances])
		file = open(Calculations.fileshort+"Summary.txt","a")
		file.writelines(["Total consensus tree length: ", str(total_dist), "\n"])
		file.close()

	# create the tree from the distance matrix and msa
	def draw_tree(n_taxa, distance_matrix):
		global fileNameF
		matrix_length = n_taxa
		n_nodes = n_taxa + n_taxa - 2
		# map matrix rows and columns to node indices
		matrix_map = [n for n in range(n_taxa)]

		tree = []
		for i in range(n_nodes):
			tree.append({})

		for u in range(n_taxa, n_nodes): # we call internal nodes "u"
			if u == n_nodes - 1:
				f, g = 0, 1 # when this is the seed node, don't have to find the next nodes to branch off
			else:
				q_matrix = Calculations.compute_q_matrix(distance_matrix, matrix_length)
				f, g = Calculations.get_lowest_off_diagonal_value_coordinate(q_matrix) # these are the next nodes to branch off

			fg_distance = distance_matrix[f][g]
			f_length = 0.5 * fg_distance + (np.sum(distance_matrix[f]) - np.sum(distance_matrix[g])) / (2.0 * (matrix_length - 2))
			g_length = fg_distance - f_length

			# add the edges and branch lengths
			tree[u][matrix_map[f]] = f_length
			tree[u][matrix_map[g]] = g_length

			# if this is the seed node, fill in the last root branch length and stop calculating
			if u == n_nodes - 1:
				tree[u][matrix_map[2]] = distance_matrix[0][2] - f_length
				break

			new_distance_matrix = np.zeros((matrix_length - 1, matrix_length - 1))

			# a and b are the old indices, i and j are the new indices
			i = 0
			new_matrix_map = [u]
			for a in range(matrix_length):
				if (a != f) and (a != g): # skip the rows to be merged
					i += 1
					j = 0

					new_matrix_map.append(matrix_map[a])

					ua_distance = 0.5 * (distance_matrix[f][a] + distance_matrix[g][a] - fg_distance)
					new_distance_matrix[0][i] = ua_distance
					new_distance_matrix[i][0] = ua_distance

					for b in range(matrix_length): # skip the columns to be merged
						if (b != f) and (b != g):
							j += 1
							new_distance_matrix[i][j] = distance_matrix[a][b]

			distance_matrix = new_distance_matrix
			matrix_map = new_matrix_map
			matrix_length = matrix_length - 1

		# save the result
		output_path = fileNameF+"Bootstraps.tre"
		Calculations.write_tree(output_path, tree, Calculations.taxon_labels)
		print(output_path + " file created!")

	#writes a new file with all the important info about the tree that was created
	def summaryFile():
		global mat
		file = open(Calculations.fileshort+"Summary.txt","w")
		L = ["MSA file name: ", Calculations.fileshort, "\n", "Score Matrix used: ", mat, "\n",
       "MSA average identity percentage: ", str(round(Calculations.score,2)*100), "\n"]
		file.writelines(L)
		file.close()
		Calculations.complete()

	#allows for final UI popup when entire calculation is complete
	def complete():
		return True