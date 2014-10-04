#!/usr/bin/env python
''' This file contains the code for the kerf function.  Inputs are:

1) Tree file in newick format.
2) Alignment file in afa format.
3) kerf threshold.
'''
from ete2 import TreeNode
from Bio import AlignIO
import time, itertools, sys, os

def pairwise_identity_belvu(seq0, seq1):
    '''Calculate pairwaise identity according to belvu.
    This is similar to the calculation done in the function 
    pairwise_identity_for_structural_alignments in bpg.common.BPGPWID
    but attempted to be more concise.  Taken from bpg.common.BPGPWID'''
    if len(seq0) != len(seq1):
        print 'Sequence lengths do not match. Exiting ...'
    return
    match = 0
    identity = 0
    for ind, char0 in enumerate(seq0):
        char1 = seq1[ind]
        if not (char0.isupper() and char1.isupper()):
            continue
        match += 1
        if char0 == char1:
            identity += 1
    try:
        pwid = float(identity)/match
    except:
        pwid = 0
    return pwid

def percent_ID(acc1, acc2, tree_msa_link, pid_cache):
    ''' This function calculates the percent ID between two leaf accession nodes.  It tries to return a 
    cached result, if it can't, it has to calculate it and store it for later.'''
    # if its identical, no need to calculate or store, just return
    if acc1 == acc2:
        return 1.0
    # order the tuple for search
    if acc1 < acc2:
        t1 = acc1
        t2 = acc2
    else:
        t1 = acc2
        t2 = acc1
    # try to return the cached result
    try:
        return pid_cache[(t1, t2)]
    except:
        pass
    # we need to calculate this
    pid_cache[(t1,t2)] = pairwise_identity_belvu(str(tree_msa_link[t1].seq), str(tree_msa_link[t2].seq))
    return pid_cache[(t1,t2)]

def validate_tree(tree_path, msa_path):
    ''' Tries to validate that the tree contains unique ids, that those ids exist in the MSA, and finally,
    it links each leaf to its msa row.'''
    # holds the link between the leaf accession and the row in that msa
    leaf_information = {}
    # first load the tree
    tree = TreeNode(newick=tree_path, format=0)
    # also load in the MSA
    msa = AlignIO.read(msa_path,'fasta')
    # next go through all of the leaves
    for leaf in tree.get_leaves():
        found = False
        for alignment_row in msa:
            if leaf.name in alignment_row.description:
                if leaf.name not in leaf_information:
                    leaf_information[leaf.name] = alignment_row
                    found = True
                    break
                else:
                    raise Exception("%s is found in the tree/MSA twice, accessions must be unique" % (leaf.name))
        if not found:
            raise Exception("%s is in the tree but not found in the MSA" % (leaf.name))
    return (tree, msa, leaf_information)

def clade_min_pid(node, tree_msa_link, pid_cache):    
    # calculates the minimum pairwise identity of a clade
    min_pid = 1.0
    for (n1, n2) in itertools.combinations(node.get_leaves(),2):
        npid = percent_ID(n1.name, n2.name, tree_msa_link, pid_cache)
        if npid < min_pid:
            min_pid = npid
    return min_pid

def kerf_parent(node, thresh, tree_msa_link, pid_cache):
    def pairwise_iterator(iterable):
        a, b = itertools.tee(iterable)
        next(b, None)
        return itertools.izip(a,b)

    # make the nodes list
    nodes = node.get_ancestors()
    nodes.insert(0,node)
    for (c, p) in pairwise_iterator(nodes):
        if (clade_min_pid(c, tree_msa_link, pid_cache) >= thresh) and (clade_min_pid(p, tree_msa_link, pid_cache) < thresh):
            return c
        elif not p.up:
            # if p is the root and the threshold is sufficient, then return the root
            if clade_min_pid(p, tree_msa_link, pid_cache) >= thresh:
                return p

def usage():
    print "*************************************************************************"
    print "USAGE:  kerf.py <path to newick tree file> <path to aligned fasta file> <kerf threshold as integer <=100 >"
    print ""
    print "kerf cuts a newick tree into subtrees based on minimum pairwise identity."
    print "*************************************************************************"
    print "\n"
    print "*************************************************************************"
    print "OUTPUTS:" 
    print "Returns 3 files per subtree cut: "
    print "One called .tree with the actual newick tree, one called .members"
    print "with the header lines of the members, and one called .msa with the"
    print "subtree multiple sequence alignment."
    print ""
    print "Returns one kerf.summary file with a summary of the kerf execution."
    print "*************************************************************************"
    return

def main(args):
    # for logging time
    start_time = time.time()
    
    # global for link between the tree node and the alignment row representing that node
    tree_msa_link = {}
    # cache for percent identities, should speed everything up.
    pid_cache = {}

    if len(args) != 4:
        usage()
        sys.exit()
    else:
        tree_path = os.path.abspath(sys.argv[1])
        msa_path = os.path.abspath(sys.argv[2])
        kerf_threshold = int(sys.argv[3]) 
        
    try:
        (tree, msa, tree_msa_link) = validate_tree(tree_path, msa_path)
    except Exception as e:
        print "Could not vaidate tree because: %s" % e
        sys.exit()

    # initialize the groups
    remaining_leaves = tree.get_leaves()
    # initialize the subtree counter
    counter = 0
    # kerf loop
    while remaining_leaves:
        subtree_file = 'kerf_%d_subtree_%d.tree' % (kerf_threshold, counter)
        subtree_members = open('kerf_%d_subtree_%d.members' % (kerf_threshold, counter), 'w')
        subtree_alignment = open('kerf_%d_subtree_%d.msa' % (kerf_threshold, counter), 'w')
        k = kerf_parent(remaining_leaves[0], kerf_threshold*1.0/100.0, tree_msa_link, pid_cache)
        k.write(outfile=subtree_file, format=0)
        for member in k.get_leaves():
            subtree_members.write(tree_msa_link[member.name].description + '\n')
            remaining_leaves.remove(member)
            subtree_alignment.write(tree_msa_link[member.name].format('fasta'))
        subtree_alignment.close()
        subtree_members.close()
        counter += 1
        #print "finished group %d" % counter
        #print "%d leaves remain" % (len(remaining_leaves))
    summary_file = open('kerf.summary','w')
    summary_file.write('kerf divided %d leaves into %d subtrees based on %d %% minimum pairwise identity.\n' % (len(tree.get_leaves()), counter, kerf_threshold))
    summary_file.write('Execution complete in %d (s).\n' % (time.time() - start_time))
    summary_file.close()

    print "kerf execution complete in %d s" % (time.time() - start_time) 
    sys.exit()

if __name__ == "__main__":
    main(sys.argv)
