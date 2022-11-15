import numpy
import os
import warnings
from collections import Counter
import pandas as pd
import math
import numpy as np

warnings.simplefilter("ignore", DeprecationWarning)
warnings.simplefilter("ignore", RuntimeWarning)


# A function to print out a matrix (e.g. a distance matrix) for human viewing
def print_matrix(title, matrix_map, matrix):
	print(title + ":")
	print("  " + "".join(["%7d" % (column) for column in matrix_map]))
	for row_i, row in enumerate(matrix):
		print("%2d" % (matrix_map[row_i]) + "".join(["%7.2f" % (element) for element in row]))
	print("")

# A function to write a tree to a file as a Newick-format string 
def write_tree(newick_path, tree, taxon_labels):
	newick_string = make_newick_string(len(tree) - 1, tree, taxon_labels) + ";"

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
			substring = make_newick_string(child_i, tree, taxon_labels)
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
def compute_q_matrix(distance_matrix):
	q_matrix = (matrix_length - 2.0) * distance_matrix

	for i in range(matrix_length):
		for j in range(matrix_length):
			if i != j:
				q_matrix[i][j] -= numpy.sum(distance_matrix[i]) + numpy.sum(distance_matrix[j])

	return q_matrix


#compute naiive prob value of residue for whole data set
def expected_val(sequences,sequence_length,taxon_labels):
	alphabet = []

	for i in range(len(sequences)):
		if i == 0:
			totals = Counter(sequences[i])
		else:
			totals.update(Counter(sequences[i]))
		
	for item in totals:
		alphabet.append(item)
	#alphabet=sorted(alphabet)
	#alphabet=alphabet[1::]
	probs = {
		letter: totals[letter] / (sequence_length*len(taxon_labels)) for letter in alphabet
	}
	
	return(probs,alphabet)

def expected_val_ND(sequences,sequence_length,taxon_labels):
	alphabet = []

	for i in range(len(sequences)):
		if i == 0:
			totals = Counter(sequences[i])
		else:
			totals.update(Counter(sequences[i]))
		
	for item in totals:
		alphabet.append(item)
	alphabet=sorted(alphabet)
	alphabet=alphabet[1::]
	probs = {
		letter: totals[letter] / (sequence_length*len(taxon_labels)) for letter in alphabet
	}
	
	return(probs,alphabet)

#get columns
def get_cols(sequences,sequence_length):
	cols_list = []
	for j in range(sequence_length):
		cols = []
		for i in sequences:
			cols.append(i[j])
			s = ''.join(cols)
		cols_list.append(s)
	return(cols_list)

#calculates naiive prob of each residue by column
def col_probs(cols_list,taxon_labels):
	x=0
	prob_list = []
	while x < len(cols_list):
		for i in range(len(cols_list)):
			if i == x:
				totals = Counter(cols_list[i])

		alphabet = []
		for item in totals:
			alphabet.append(item)

		probs = {
			letter: totals[letter] / len(taxon_labels) for letter in sorted(alphabet)
		}
		x+=1
		prob_list.append(probs)
	df = pd.DataFrame.from_dict(prob_list).fillna(0)
	return(df)

def realdist(n_taxa,msa): #calculating sigma(s1,s2) for every row and replacing with score
	#why don't we need sequence_length input here?
	#make distance matrix of zeros
	matrix_length = n_taxa
	distance_matrix = numpy.zeros((matrix_length, matrix_length))
	seq = 0
	counter = 0
	while seq < (n_taxa-1):
		if counter == n_taxa:
			seq += 1
			print(seq)
			counter = 0

		for j in range(sequence_length):
			msa_1 = msa[seq][j]#first residue
			msa_2 = msa[counter][j]#second reside
			if (msa_1 == '-' and msa_2 == '-'): #two gaps calculate nothing
				continue
			else:
				comp = df[msa_1][msa_2] #find in score mat
				distance_matrix[seq][counter] += comp #that score is spot in new_mat
			
		counter += 1
	return distance_matrix

def randdist(n_taxa,msa,sequence_length):
	matrix_length = n_taxa
	distance_matrix = numpy.zeros((matrix_length, matrix_length))
	seq = 0
	counter = 0
	while seq < (n_taxa-1):
		comb_list = []
		if counter == n_taxa:
			seq += 1
			print(seq)
			counter = 0

		for j in range(sequence_length):
			msa_1 = msa[seq][j]#first residue
			msa_2 = msa[counter][j]#second reside
			if (msa_1 == '-' and msa_2 == '-'): #two gaps calculate nothing
				continue
			else:
				if (msa_1,msa_2) in comb_list: #combo already found, dont add it
					continue
				else:
					comb_list.append((msa_1,msa_2)) #add to combo list found in two seq
		
		for val in comb_list:
			comb_score = df[val[0]][val[1]] #find in score mat
			num1_occ = np.count_nonzero(msa[seq]==val[0]) #count num occ of first val (in first seq)
			num2_occ = np.count_nonzero(msa[counter]==val[1]) #count num occ of sec val (in sec seq)
			gap_count1 = np.count_nonzero(msa[seq]==45) #get num of gaps seq 1
			#gap_count2 = np.count_nonzero(msa[counter]==45) #get num of gaps seq 2
			rand_score = ((comb_score*num1_occ*num2_occ))/sequence_length - (gap_count1*-2)
			distance_matrix[seq][counter] += rand_score #that score is spot in new_mat
			
		counter += 1
	return distance_matrix

def upper_limits(n_taxa, msa, sequence_length):
	scores = realdist(n_taxa, msa)
	matrix = numpy.zeros((n_taxa, n_taxa))
	i = 0
	while i < (n_taxa-1):
		j = 0
		while j < (sequence_length-1):
			matrix[i, j] = (scores[i,i] + scores[j,j]) / 2
			j = j + 1
		i = i + 1
	return matrix

def dist_mat(n_taxa,msa,sequence_length):
	similarities = realdist(n_taxa,msa)
	randoms = randdist(n_taxa,msa,sequence_length)
	norm_scores = similarities - randoms
	upper = upper_limits(n_taxa, msa, sequence_length)
	upper_norm = upper - randoms
	raw_dist = -math.log(norm_scores / upper_norm) * 100
	#c = 
	#distance_matrix = c * raw_dist
	return raw_dist


# Read in the multiple sequence alignment
sequences, sequence_length, taxon_labels, msa = read_phylip("real_test2.phy")


#TESTING EXPECTED VALUE FOR WHOLE SET
exp_probs,alphabet = expected_val(sequences,sequence_length,taxon_labels)

exp_probs_ND,alphabet2 = expected_val_ND(sequences,sequence_length,taxon_labels)



cols_list = get_cols(sequences,sequence_length)

df = col_probs(cols_list,taxon_labels)
#print(df)
#print(df.iloc[0,1])


#TODO:
#REPLACE DISTANCE FUNCTION WITH THE ONE FROM THE PAPER
#find each pair combo, find that score in Score matrix (BLOSUM, whatev), then put that score in a list, sum the scores for those two, then do it for each sequence combo
#resulting in a matrix of scores between each sequence

input = input('Score Matrix?: ')
print(input)
if input == '1':
	arr = pd.read_csv('BLOSUM30.csv', header=None).values #importing BLOSUM 30 (need to clean this)
if input == '2':
	arr = pd.read_csv('PAM500.csv', header=None).values
if input == '3':
	arr = pd.read_csv('BLOSUM62.csv', header=None).values

print(arr)
with_pen = np.empty([21,21])
for i in range(len(arr)):	
	with_pen[i] = np.append(arr[i], [-2])
with_pen[20] = np.append(np.repeat(-2.0, 20), 0.0)
print(with_pen)

labs = 'ARNDCQEGHILKMFPSTWYV-'
labels = numpy.fromstring(labs, dtype = "uint8") 
df = pd.DataFrame(with_pen,columns=labels,index= labels) #make scoring dataframe
print(df)


n_taxa = len(taxon_labels)
n_nodes = n_taxa + n_taxa - 2
#matrix_length = n_taxa
distance_matrix = dist_mat(n_taxa, msa, sequence_length)
print(distance_matrix)

# real_distance_matrix = realdist(n_taxa,msa)

# print(real_distance_matrix)
# print(real_distance_matrix[0])
# print(real_distance_matrix[1])

# #random calculation equation
# #(1/length) sum((each possible residue type)*number of the that residue in i * number of that residue in j))
# #  - number of gaps *penalty #random ij or sigma r

# rand_dist_matrix = randdist(n_taxa,msa,sequence_length)
# print(rand_dist_matrix)
# print(rand_dist_matrix[0])
# print(rand_dist_matrix[1])

# p_distance_matrix=rand_dist_matrix #temporary to prevent further errors of stuff we havent changed
# #to be replaced
# for i in range(n_taxa): 
# 	msa_i = msa[i]
# 	for j in range(n_taxa):
# 		msa_j = msa[i]
# 		identity = float(numpy.sum(msa_i == msa_j))
# 		p_distance_matrix[i][j] = 1.0 - identity / sequence_length
# #print(p_distance_matrix)

# matrix_map = [n for n in range(n_taxa)] # mapping matrix rows and columns to node indices
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
		q_matrix = compute_q_matrix(distance_matrix)
		f, g = get_lowest_off_diagonal_value_coordinate(q_matrix) # these are the next nodes to branch off

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
write_tree(output_path, tree, taxon_labels)