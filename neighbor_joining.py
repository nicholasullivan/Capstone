import numpy
import os
import warnings
from collections import Counter
import pandas as pd

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
			letter: totals[letter] / len(taxon_labels) for letter in alphabet
		}
		x+=1
		prob_list.append(probs)
	#print(prob_list)
	df = pd.DataFrame.from_dict(prob_list).fillna(0)
	return(df)


# Read in the multiple sequence alignment
sequences, sequence_length, taxon_labels, msa = read_phylip("C:/Users/keerp/Downloads/real_test2.phy")


#TESTING EXPECTED VALUE FOR WHOLE SET
exp_probs,alphabet = expected_val(sequences,sequence_length,taxon_labels)
#print(exp_probs)


cols_list = get_cols(sequences,sequence_length)

df = col_probs(cols_list,taxon_labels)
#print(df)

#getting transition matrix
sequences1 = sequences[0]
prob_matrix = {}
for i in alphabet:
    prob_matrix[i] = {}
    for j in alphabet:
        prob_matrix[i][j] = 0.0
		
for i, j in zip(sequences1[:-1], sequences1[1:]):
	prob_matrix[i][j] += 1

print(prob_matrix)


#Needs replacement
n_taxa = len(taxon_labels)
n_nodes = n_taxa + n_taxa - 2

matrix_length = n_taxa
p_distance_matrix = numpy.zeros((matrix_length, matrix_length))




#TODO:
#REPLACE DISTANCE FUNCTION WITH THE ONE FROM THE PAPER

for i in range(n_taxa): 
	msa_i = msa[i]
	for j in range(n_taxa):
		msa_j = msa[j]
		identity = float(numpy.sum(msa_i == msa_j))
		p_distance_matrix[i][j] = 1.0 - identity / sequence_length

matrix_map = [n for n in range(n_taxa)] # mapping matrix rows and columns to node indices
distance_matrix = -0.75 * numpy.log(1.0 - 1.3333333333 * p_distance_matrix) # using the Jukes-Cantor 1969 (JC69) model

#print_matrix("P-distance matrix", matrix_map, p_distance_matrix)
#print_matrix("Distance matrix", matrix_map, distance_matrix)

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
output_path = "C:/Users/keerp/Downloads/nj.tree" 
write_tree(output_path, tree, taxon_labels)




#TO DO:
#Figure out code (HOW TO ADD OUR OWN SCORING MATRICES???)
#FIND BOOTSTRAPPING!

