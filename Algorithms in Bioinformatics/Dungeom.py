#!/usr/bin/env python3

'''
                               Dungeom

Description:
Dungeom is a de Bruijn graph-based assembler that gets a DNA sequence and a
kmer size as inputs and outputs the reassembled sequence and "identical"
or "not identical" if the reassembled sequence is identical or not with the
input sequence.

Procedure:

1.The DNA sequence input is checked as to its length and the existance of
non-standard nucleotides. The kmer size input is also checked as to whether it
is a positive integer smaller than the length of the input sequence.
2.The sequence is broken into small fragments (kmers) and the unique kmers are
saved into a set to secure their uniqueness (nodes). Also a list of tuples with
the edges of the nodes is created that is later used in the construction of the
de Bruijn graph (edges).
3.All the edges are shuffled in order to secure the randomness of the assembly
except the first one, in order to keep the starting point. The output will be
the list of the shuffled edge tuples (graph).
4.The process of the assembly is based on finding the Eulerian path from the graph,
that is the trail in the de Bruijn graph which visits every edge exactly once.
That requires to keep track of the number of the edges for each node in a dictionary
(freq_dict) and finish the assembly when all the edges are used (values in the
dictionary are 0). The outcome of the procedure is a list (tour_list) with the
nodes that were visited in the order that they were visited.
5.Finally, the nodes in the tour list are merged into the assembled sequence.
Depending on whether the assembled sequence is identical to the input one,
"Assembled sequence: identical" or "not identical" will also be printed.

Usage: python Dungeom.py sequence kmer_size

Example: python Dungeom.py ATGGGTCA 2
'''


import sys
import argparse
from random import shuffle




''''
Input check:
The DNA sequence input is checked as to its length (max.1000 bp) and the
existance of non-standard nucleotides. The kmer size input is also checked as to
whether it is a positive integer smaller than the length of the input sequence.
'''

def check_positive_kmer_size(value):
    if not value.isdigit():
        raise argparse.ArgumentTypeError("{} is an invalid positive integer value".format(value))
    ivalue = int(value)
    if ivalue <= 0:
         raise argparse.ArgumentTypeError("{} is an invalid positive integer value".format(value))
    return ivalue



def check_dna_sequence(input_seq):
    for nucleotide in input_seq:
        if nucleotide not in "GCTAgcta":
            raise argparse.ArgumentTypeError("Non-standard nucleotides in the sequence")


parser = argparse.ArgumentParser(description=' Randomly produce kmers of a DNA sequence string\
 and reassemble it using a de Bruijn graph method.')

parser.add_argument('sequence string', metavar='seq', type=check_dna_sequence,
                    help='a sequence string of no more than 1000 bp')
parser.add_argument('kmer size', metavar='kmer', type=check_positive_kmer_size,
                    help='kmer size')
args = parser.parse_args()

seq = sys.argv[1].upper()
kmer = int(sys.argv[2])


sys.setrecursionlimit(len(seq)) #corrects the maximum recursion error created by
                                #the find_eulerian_tour and find_tour functions



if kmer >= len(seq):
    print("Dungeom.py: error: argument kmer: kmer size should be smaller than the sequence length")
    sys.exit(1)



'''
Functions
'''



def de_bruijn_convert(str, k):
    #The sequence is broken into small fragments (kmers) and the unique kmers are
    #saved into a set to secure their uniqueness (nodes). Also a list of tuples with
    #the edges of the nodes is created that is later used in the construction of the
    #de Bruijn graph (edges).
    edges = []
    nodes = set()
    for i in range(len(str) - k):
        edges.append((str[i:i+k], str[i+1:i+k+1]))
        nodes.add(str[i:i+k])
        nodes.add(str[i+1:i+k+1])
    return edges, nodes



#def visualize_de_bruijn_graph(str, k):
    #The de Bruijn graph is constructed and is saved into a pdf file using the
    #graphviz package.
    #from graphviz import Digraph
    #edges, nodes = de_bruijn_convert(seq, kmer)
    #f = Digraph('de_Bruijn_graph')
    #f.attr(rankdir='LR', size='15')
    #for node in nodes:
    #    f.node(node)
    #for lf, rg in edges:
    #    f.edge(lf, rg)
    #f.render('de_bruijn_out')



def shuffled_edges(graph):
    #All the edge tuples (graph) are shuffled in order to secure the randomness of the
    #assembly except the first one, in order to keep the starting point. The
    #output will be the list of the shuffled edge tuples (graph).
    first_element = True
    for edge in graph:
        if first_element:
            first_edge = edge
            graph.remove(edge)
            first_element = False
    shuffle(graph)
    graph.insert(0, first_edge)
    return graph




def dictionary_of_node_edge_frequencies(edges, nodes):
    #The process of the assembly is based on finding the Eulerian path from the graph,
    #that is the trail in the de Bruijn graph which visits every edge exactly once.
    #That requires to keep track of the number of the edges for each node in a dictionary
    #(freq_dict).
    nodes = list(nodes)
    result = [0 for i in range(len(nodes))]
    for node in nodes:
        for (a,b) in edges:
            if node == a or node == b:
                result[nodes.index(node)] += 1
    edge_freq = dict(zip(nodes, result))
    return edge_freq



def find_eulerian_tour(graph):
    #The assembly will be finished when all the edges are used (values in the
    #dictionary are 0). The outcome of the procedure is a list (tour_list) with the
    #nodes that were visited in the order that they were visited.
    edges, nodes = de_bruijn_convert(seq, kmer)
    freq_dict = dictionary_of_node_edge_frequencies(graph, nodes)
    while sum(freq_dict.values()) != 0:
        edges, nodes = de_bruijn_convert(seq, kmer)
        graph = shuffled_edges(edges)
        freq_dict = dictionary_of_node_edge_frequencies(graph, nodes)
        tour = []
        find_tour(graph[0][0], graph, tour, freq_dict)
    return tour
def find_tour(current,Edges,tour, kmer_dict):
    for (a,b) in Edges:
        if a==current:
            Edges.remove((a,b))
            kmer_dict[a] -= 1
            kmer_dict[b] -= 1
            find_tour(b,Edges,tour, kmer_dict)
    tour.insert(0,current)




def eulerian_tour_to_superstring(graph):
    #Finally, the nodes in the tour list are merged into the assembled sequence.
    #Depending on whether the assembled sequence is identical to the input one,
    #"Assembled sequence: identical" or "not identical" will also be printed.
    superstring = ""
    first_word = True
    for word in graph:
        if first_word:
            superstring += word
            first_word = False
        else:
            superstring += word[-1]
    return superstring


'''
Main
'''


edges, nodes = de_bruijn_convert(seq, kmer)
#de_Bruijn_graph = visualize_de_bruijn_graph(seq, kmer) #make sure the graphviz package is installed
graph = shuffled_edges(edges)
freq_dict = dictionary_of_node_edge_frequencies(graph, nodes)
tour_list = find_eulerian_tour(graph)
assembled_seq = eulerian_tour_to_superstring(tour_list)
print(assembled_seq)
if assembled_seq == seq:
    print("Assembled sequence: identical")
else:
    print("Assembled sequence: not identical")
