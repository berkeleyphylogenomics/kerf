kerf
====

phylogenetic tree cutting algorithm

Kerf cuts a tree into subtrees based on minimum pairwise identity.

*************************************************************************
REQUIREMENTS:  requires AlignIO from BioPython, TreeNode from ete2
*************************************************************************


*************************************************************************
USAGE:  kerf.py <path to newick tree file> <path to aligned fasta file> <kerf threshold as integer <=100 >

kerf cuts a newick tree into subtrees based on minimum pairwise identity.
*************************************************************************


*************************************************************************
OUTPUTS:
Returns 3 files per subtree cut:
One called .tree with the actual newick tree, one called .members
with the header lines of the members, and one called .msa with the
subtree multiple sequence alignment.

Returns one kerf.summary file with a summary of the kerf execution.
*************************************************************************
