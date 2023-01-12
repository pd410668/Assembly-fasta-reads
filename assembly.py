#!/usr/bin/python3
# coding: utf-8

import sys
from Bio import SeqIO

def read_fasta(input_file):
    '''Create list of reads from fasta file'''
    reads = []
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        reads.append(sequence)
    return reads


def neighbors1mm(kmer, alpha):
    ''' Generate all neighbors at Hamming distance 1 from kmer '''
    neighbors = []
    for j in range(len(kmer) - 1, -1, -1):
        oldc = kmer[j]
        for c in alpha:
            if c == oldc: continue
            neighbors.append(kmer[:j] + c + kmer[j + 1:])
    return neighbors


def kmerHist(reads, k):
    ''' Return k-mer histogram and average # k-mer occurrences '''
    kmerhist = {}
    for read in reads:
        for kmer in [read[i:i + k] for i in range(len(read) - (k - 1))]:
            kmerhist[kmer] = kmerhist.get(kmer, 0) + 1
    return kmerhist


def correct1mm(read, k, kmerhist, alpha, thresh):
    ''' Return an error-corrected version of read.  k = k-mer length.
        kmerhist is kmer count map.  alpha is alphabet.  thresh is
        count threshold above which k-mer is considered correct. '''
    # Iterate over k-mers in read
    for i in range(len(read) - (k - 1)):
        kmer = read[i:i + k]
        # If k-mer is infrequent...
        if kmerhist.get(kmer, 0) <= thresh:
            # Look for a frequent neighbor
            for newkmer in neighbors1mm(kmer, alpha):
                if kmerhist.get(newkmer, 0) > thresh:
                    # Found a frequent neighbor; replace old kmer
                    # with neighbor
                    read = read[:i] + newkmer + read[i + k:]
                    break
    # Return possibly-corrected read
    return read


class DeBruijnGraph:
    ''' De Bruijn directed multigraph built from a collection of
        strings. User supplies strings and k-mer length k.  Nodes
        are k-1-mers.  An Edge corresponds to the k-mer that joins
        a left k-1-mer to a right k-1-mer. '''

    @staticmethod
    def chop(st, k):
        ''' Chop string into k-mers of given length '''
        for i in range(len(st) - (k - 1)):
            yield (st[i:i + k], st[i:i + k - 1], st[i + 1:i + k])

    class Node:
        ''' Node representing a k-1 mer.  Keep track of # of
            incoming/outgoing edges so it's easy to check for
            balanced, semi-balanced. '''

        def __init__(self, km1mer):
            self.km1mer = km1mer
            self.nin = 0
            self.nout = 0

        def isSemiBalanced(self):
            return abs(self.nin - self.nout) == 1

        def isBalanced(self):
            return self.nin == self.nout

        def __hash__(self):
            return hash(self.km1mer)

        def __str__(self):
            return self.km1mer

    def __init__(self, strIter, k, circularize=False):
        ''' Build de Bruijn multigraph given list of reads and k-mer
            length k '''
        self.k = k
        self.G = {}  # multimap from nodes to neighbors
        self.nodes = {}  # maps k-1-mers to Node objects
        for st in strIter:
            if circularize:
                st += st[:k - 1]
            for kmer, km1L, km1R in self.chop(st, k):
                nodeL, nodeR = None, None
                if km1L in self.nodes:
                    nodeL = self.nodes[km1L]
                else:
                    nodeL = self.nodes[km1L] = self.Node(km1L)
                if km1R in self.nodes:
                    nodeR = self.nodes[km1R]
                else:
                    nodeR = self.nodes[km1R] = self.Node(km1R)
                self.G.setdefault(nodeL, {})
                nodeR.nin += 1
                nodeL.nout += 1
                if nodeR in self.G[nodeL]:
                    self.G[nodeL][nodeR] += 1
                else:
                    self.G[nodeL][nodeR] = 1
        # Iterate over nodes; tally # balanced, semi-balanced, neither
        self.nsemi, self.nbal, self.nneither = 0, 0, 0
        # Keep track of head and tail nodes in the case of a graph with
        # Eularian walk (not cycle)
        self.head, self.tail = None, None
        for node in iter(self.nodes.values()):
            if node.isBalanced():
                self.nbal += 1
            elif node.isSemiBalanced():
                if node.nin == node.nout + 1:
                    self.tail = node
                if node.nin == node.nout - 1:
                    self.head = node
                self.nsemi += 1
            else:
                self.nneither += 1

    def nnodes(self):
        ''' Return # nodes '''
        return self.nodes

    def nedges(self):
        ''' Return # edges '''
        return len(self.G)

    def hasEulerianWalk(self):
        ''' Return true iff graph has Eulerian walk. '''
        return self.nneither == 0 and self.nsemi == 2

    def hasEulerianCycle(self):
        ''' Return true iff graph has Eulerian cycle. '''
        return self.nneither == 0 and self.nsemi == 0

    def isEulerian(self):
        ''' Return true iff graph has Eulerian walk or cycle '''
        # technically, if it has an Eulerian walk
        return self.hasEulerianWalk() or self.hasEulerianCycle()

    def eulerianWalkOrCycle(self):
        ''' Find and return sequence of nodes (represented by
            their k-1-mer labels) corresponding to Eulerian walk
            or cycle '''
        assert self.isEulerian()
        g = self.G
        if self.hasEulerianWalk():
            g = g.copy()
            g.setdefault(self.tail, []).append(self.head)
        # graph g has an Eulerian cycle
        tour = []
        src = next(iter(g.keys()))  # pick arbitrary starting node

        def __visit(n):
            while len(g[n]) > 0:
                dst = g[n].pop()
                __visit(dst)
            tour.append(n)

        __visit(src)
        tour = tour[::-1][:-1]  # reverse and then take all but last node

        if self.hasEulerianWalk():
            # Adjust node list so that it starts at head and ends at tail
            sti = tour.index(self.head)
            tour = tour[sti:] + tour[:sti]

        # Return node list
        return list(map(str, tour))

    def simplification(self):
        node_list = list(self.G.keys())
        flg = True
        visited = {}
        while flg == True:
            flg = False
            for node in node_list:
                visited[node] = True
                if len(self.G[node]) != 1:
                    continue
                else:
                    dist = list(self.G[node].keys())[0]
                    for i in self.G:
                        if dist in self.G[i] and i != node:
                            break
                    else:
                        flg = True
                        new = node.km1mer + dist.km1mer[self.k - 2:]
                        new_node = self.Node(new)
                        new_node.nin = node.nin
                        new_node.nout = dist.nout
                        self.nodes[new] = new_node
                        if dist in self.G:
                            self.G[new_node] = self.G[dist]
                            del self.G[dist]
                        del self.nodes[dist.km1mer]
                        for i in list(self.G.keys()):
                            if node in self.G[i]:
                                self.G[i][new_node] = self.G[i][node]
                                del self.G[i][node]
                        del self.G[node]
                        del self.nodes[node.km1mer]
                        node_list = list(filter(lambda x: x not in visited, self.G.keys()))
                        break

    def tips(self):
        flg = True
        while flg == True:
            flg = False
            for node in list(self.nodes.values()):
                if node.nin == 0:
                    if node.nout < 2 and len(node.km1mer) < 2 * self.k:
                        flg = True
                        for i in self.G:
                            if i in self.G[node]:
                                i.nin -= self.G[node][i]
                        del self.G[node]
                        del self.nodes[node.km1mer]

                if node.nout == 0:
                    if node.nin < 2 and len(node.km1mer) < 2 * self.k:
                        flg = True
                        for i in self.G:
                            if node in self.G[i]:
                                i.nout -= self.G[i][node]
                                del self.G[i][node]
                        del self.nodes[node.km1mer]
                        if node in self.G:
                            del self.G[node]

    def remove_buble(self):
        visited = []
        for x in list(self.G.keys()):
            if x in visited: continue
            src = {}
            queue = [x]

            def remove(path):
                for node in path:
                    if node in self.G:
                        for i in self.G[node]:
                            i.nin -= self.G[node][i]
                            if i in src and src[i] == node:
                                del src[i]
                        del self.G[node]
                    for i in self.G:
                        if node in self.G[i]:
                            i.nout -= self.G[i][node]
                            del self.G[i][node]
                    del self.nodes[node.km1mer]

            while queue:
                node = queue.pop()
                if node not in self.G: continue
                for i in self.G[node]:
                    if i not in visited:
                        queue.append(i)
                        visited.append(i)
                        src[i] = node
                    else:
                        if i not in src: continue
                        path1 = [node]
                        path2 = [src[i]]
                        flg = True
                        while (set(path1) & set(path2)) == set():
                            if path1[-1] in src:
                                path1.append(src[path1[-1]])
                            if path2[-1] in src:
                                path2.append(src[path2[-1]])
                            if path1[-1] not in src and path2[-1] not in src or i in path1:
                                flg = False
                                break
                        if flg:
                            start = list(set(path1) & set(path2))[0]
                            path1 = path1[:path1.index(start)]
                            path2 = path2[:path2.index(start)]
                            w1 = sum(sum(self.G[i][j] for j in self.G[i]) for i in path1)
                            w2 = sum(sum(self.G[i][j] for j in self.G[i]) for i in path2)
                            if w1 > w2:
                                remove(path2)
                            else:
                                remove(path1)

    def refine(self):
        self.simplification()
        self.tips()
        self.simplification()
        self.remove_buble()
        self.simplification()
        self.tips()
        self.simplification()


k = 15
thresh = 3
contigs = []

L = read_fasta(sys.argv[1])

# poprawianie błędów
khist = kmerHist(L, k)
corrected_reads = []
for each in L:
    correct_read = correct1mm(each, k, khist, 'ACGT', thresh)
    corrected_reads.append(correct_read)
g = DeBruijnGraph(corrected_reads, k)
g.refine()
for i in g.nodes:
    if len(i)> 80:
        contigs.append(i)
with open(sys.argv[2],'w') as out:
    for i, item in enumerate(contigs):
        out.write(">contig_" + str(i) + '\n')
        out.write(item+'\n')

