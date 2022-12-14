"""
This is a code to calculate the similarities between DNA sequences. It allows the user to input their phylip file and choose the 
scoring matrix they want. The sequences are then compared using the choosen scoring matrix, go through a neighbor joining and tree 
building process, and, finally, results are confirmed via boostrapping.

Authors- Troy Hofstrand, Nicholas Sullivan, and Nikolaus Ryczek
Emails- thofstrand@slu.edu, nsullivan@slu.edu, nryczek@slu.edu

Last Date Updated- 11/29/2022
"""
import numpy
import os
import warnings
from collections import Counter
import pandas as pd
import math
import numpy as np
from tqdm import tqdm
from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter.filedialog import askopenfilename

warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", RuntimeWarning)

class Calculations:
	# A function to print out a matrix (e.g. a distance matrix) for human viewing
	def print_matrix(title, matrix_map, matrix):
		print(title + ":")
		print("  " + "".join(["%7d" % (column) for column in matrix_map]))
		for row_i, row in enumerate(matrix):
			print("%2d" % (matrix_map[row_i]) + "".join(["%7.2f" % (element) for element in row]))
		print("")

	# A function to write a tree to a file as a Newick-format string 
	def write_tree(newick_path, tree, taxon_labels):
		newick_string = Calculations.make_newick_string(len(tree) - 1, tree, taxon_labels) + ";"

		newick_file = open(newick_path, "w")
		newick_file.write(newick_string)
		newick_file.close()

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
		phylip_file = open(phylip_path)
		phylip_header = phylip_file.readline().strip().split()

		n_taxa = int(phylip_header[0])
		sequence_length = int(phylip_header[1])
		sequences = []
		sequence_labels = []
		msa = numpy.zeros((n_taxa, sequence_length), dtype = "uint8")

		for i in range(n_taxa):
			line = phylip_file.readline()
			label, sequence = line.strip().split()

			sequence_labels.append(label)
			msa[i] = numpy.fromstring(sequence, dtype = "uint8")
			sequences.append(sequence)
		return (sequences, sequence_length, sequence_labels, msa)

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
					q_matrix[i][j] -= numpy.sum(distance_matrix[i]) + numpy.sum(distance_matrix[j])

		return q_matrix




	def realscore(df,n_taxa,msa,sequence_length): #calculating sigma(s1,s2) for every row and replacing with score
		#why don't we need sequence_length input here?
		#make distance matrix of zeros
		matrix_length = n_taxa
		distance_matrix = numpy.zeros((matrix_length, matrix_length))
		seq = 0
		counter = 0

		while True:
			if counter == n_taxa:
				print('seq',seq)
				seq += 1
				if seq == n_taxa:
					return distance_matrix
				counter = 0

			for j in range(sequence_length):
				msa_1 = msa[seq][j]#first residue
				msa_2 = msa[counter][j]#second reside
				if (msa_1 == 45 and msa_2 == 45): #two gaps calculate nothing
					continue
				else:
					comp = df[msa_1][msa_2] #find in score mat
					distance_matrix[seq][counter] += comp #that score is spot in new_mat
			print('counter',counter)
			counter += 1


	def randscore(df,n_taxa,msa,sequence_length):
		matrix_length = n_taxa
		distance_matrix = numpy.zeros((matrix_length, matrix_length))
		seq = 0
		counter = 0
		while True:
			rand_score=0
			comb_list = []
			seq1 = []
			seq2 = []
			if counter == n_taxa:
				print('seq',seq)
				seq += 1
				if seq == n_taxa:
					return distance_matrix
				counter = 0
			for j in range(sequence_length):
				msa_1 = msa[seq][j]#first residue
				msa_2 = msa[counter][j]#second reside
				if (msa_1 == 45 and msa_2 == 45): #two gaps calculate nothing
					continue
				
				else:
					if (msa_1,msa_2) not in comb_list: #combo already found, dont add it
						comb_list.append((msa_1,msa_2))
				seq1.append(msa_1)
				seq2.append(msa_2)
			
			for val in comb_list:
				comb_score = df[val[0]][val[1]] #find in score mat
				num1_occ = seq1.count(val[0]) #count num occ of first val (in first seq)
				num2_occ = seq2.count(val[1]) #count num occ of sec val (in sec seq)
				gap_count1 = seq1.count(45) #get num of gaps seq 1
				gap_count2 = seq2.count(45)#get num of gaps seq 2
				rand_score += (comb_score*num1_occ*num2_occ)  
			distance_matrix[seq][counter] = (rand_score/len(seq1))- ((gap_count1+gap_count2)*2) #that score is spot in new_mat
			#print('counter',counter)
			counter += 1


	def identityscore(n_taxa, reals):
		scores = reals
		matrix = numpy.zeros((n_taxa, n_taxa))
		i = 0
		while i < (n_taxa):
			for j in range(n_taxa):
				matrix[i, j] = (scores[i,i] + scores[j,j]) / 2
			i += 1
		return matrix

	def dist_mat(df,n_taxa,msa,sequence_length):
		real = Calculations.realscore(df,n_taxa,msa,sequence_length)
		rand = Calculations.randscore(df,n_taxa,msa,sequence_length)
		norm_scores = np.subtract(real , rand)
		print("Real score of first row:", real)
		print("Rand score of first row:", rand[0])
		print("Norm scores of first row:", norm_scores[0])
		identity = Calculations.identityscore(n_taxa,real)
		upper_norm = np.subtract(identity , rand)
		print("Identity score first row:", identity[0])
		print("Norm Upper Limit scores of first row:", upper_norm[0])
		raw_dist = -np.log(np.divide(norm_scores, upper_norm))*100
		print("Raw Distance from first sequence to all other sequences:", raw_dist[0])
		#c =  
		#distance_matrix = c * raw_dist
		return raw_dist


	#TODO:
		#Fix randdist
		#add in fast file readin


	# Read in the multiple sequence alignment
	#sequences, sequence_length, taxon_labels, msa = read_phylip("real_test2.phy")
	def file_select():
		print('Hello!\n'
		'Welcome to DNA similarity calculator!\n'
		'Select your desired phylip file:')
		Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
		filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
		print(filename+ ' selected\n')
		Calculations.sequences, Calculations.sequence_length, Calculations.taxon_labels, Calculations.msa = Calculations.read_phylip(filename)

		"""print('Scoring Matrix Options:\n'
		'A- BLOSUM30\n'
		'B- BLOSUM40\n'
		'C- BLOSUM50\n'
		'D- BLOSUM62\n'
		'E- PAM300\n'
		'F- PAM400\n'
		'G- PAM500\n'
		)"""
	def matrix_selection(value):
		if value == 'A':
			arr = pd.read_csv('Matrices/BLOSUM30.csv', header=None).values #importing BLOSUM 30 (need to clean this)
		if value == 'B':
			arr = pd.read_csv('Matrices/BLOSUM40.csv', header=None).values
		if value == 'C':
			arr = pd.read_csv('Matrices/BLOSUM50.csv', header=None).values
		if value == 'D':
			arr = pd.read_csv('Matrices/BLOSUM62.csv', header=None).values
		if value == 'E':
			arr = pd.read_csv('Matrices/PAM300.csv', header=None).values
		if value == 'F':
			arr = pd.read_csv('Matrices/PAM400.csv', header=None).values
		if value == 'G':
			arr = pd.read_csv('Matrices/PAM500.csv', header=None).values

		print(arr)
		with_pen = np.empty([21,21])
		for i in range(len(arr)):	
			with_pen[i] = np.append(arr[i], [-2])
		with_pen[20] = (np.repeat(-2.0, 21))

		labs = 'ARNDCQEGHILKMFPSTWYV-'
		labels = numpy.fromstring(labs, dtype = "uint8") 
		df = pd.DataFrame(with_pen,columns=labels,index= labels) #make scoring dataframe
		print(df)


		n_taxa = len(Calculations.taxon_labels)
		n_nodes = n_taxa + n_taxa - 2
		matrix_length = n_taxa
		distance_matrix = Calculations.dist_mat(df,n_taxa, Calculations.msa, Calculations.sequence_length)
		print("Raw Distance Matrix:", distance_matrix)


# p_distance_matrix = rand_dist_matrix 
# #to be replaced
# for i in range(n_taxa): 
# 	msa_i = msa[i]
# 	for j in range(n_taxa):
# 		msa_j = msa[i]
# 		identity = float(numpy.sum(msa_i == msa_j))
# 		p_distance_matrix[i][j] = 1.0 - identity / sequence_length
# #print(p_distance_matrix)
		
		matrix_map = [n for n in range(n_taxa)] # mapping matrix rows and columns to node indices
		# distance_matrix = -0.75 * numpy.log(1.0 - 1.3333333333 * p_distance_matrix) # using the Jukes-Cantor 1969 (JC69) model

		# #print_matrix("P-distance matrix", matrix_map, p_distance_matrix)
		# #print_matrix("Distance matrix", matrix_map, distance_matrix)

		tree = []
		for i in range(n_nodes):
			tree.append({})

		for u in range(n_taxa, n_nodes): # we call internal nodes "u"
			if u == n_nodes - 1:
				f, g = 0, 1 # when this is the seed node, don't have to find the next nodes to branch off
			else:
				q_matrix = Calculations.compute_q_matrix(distance_matrix, matrix_length)
				#print('q',q_matrix)
				f, g = Calculations.get_lowest_off_diagonal_value_coordinate(q_matrix) # these are the next nodes to branch off

			fg_distance = distance_matrix[f][g]
			f_length = 0.5 * fg_distance + (numpy.sum(distance_matrix[f]) - numpy.sum(distance_matrix[g])) / (2.0 * (matrix_length - 2))
			g_length = fg_distance - f_length

			# add the edges and branch lengths
			tree[u][matrix_map[f]] = f_length
			tree[u][matrix_map[g]] = g_length

			# if this is the seed node, fill in the last root branch length and stop calculating
			if u == n_nodes - 1:
				tree[u][matrix_map[2]] = distance_matrix[0][2] - f_length
				break

			new_distance_matrix = numpy.zeros((matrix_length - 1, matrix_length - 1))

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

			#print_matrix("Distance matrix", new_matrix_map, new_distance_matrix)

			distance_matrix = new_distance_matrix
			matrix_map = new_matrix_map
			matrix_length = matrix_length - 1

		# save the result
		output_path = "nj.tree" 
		def get_file():
			return output_path
		Calculations.write_tree(output_path, tree, Calculations.taxon_labels)
		
