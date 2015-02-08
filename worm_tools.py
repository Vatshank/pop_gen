import glob
import numpy as np
import subprocess
from Bio import SeqIO, SeqRecord, Phylo, AlignIO
from cStringIO import StringIO
from Bio.Align import MultipleSeqAlignment
from cogent.phylo.maximum_likelihood import ML
from cogent.evolve.models import F81
from cogent import LoadTree, LoadSeqs

### does the genomic position, specified by contig and coordinate, align to expectatus
def outgroup_aligned(contig, pos_start, pos_stop, aln=None, refseqs=None):
    if aln==None:
        return False
    else:
        ids = [seqrec.id for seqrec in aln]
        mdis = np.mean(np.array(aln[ids.index(refseqs[0])])!=np.array(aln[ids.index(refseqs[1])]))
        print "distance to outgroup", mdis
        if mdis>0.02:
            return True
        else:
            return False

### parse the file with clade designations and return a dictionary
def get_strain_types():
    strains = {}
    with open('/ebio/ag-neher/share/users/rneher/PpacificusPopulationGenetics/data/clades.txt', 'r') as f:
        for line in f:
            strain_type = line.split()[2]
            strain_name = line.split()[1]
            if (strain_type not in strains):
                strains[strain_type]=[strain_name]
            else:
                strains[strain_type].append(strain_name)
    strains['all']=[]
    with open('/ebio/ag-neher/share/users/rneher/PpacificusPopulationGenetics/data/all_strains.txt', 'r') as f:
        for line in f:
            strain_name = line.strip()
            strains['all'].append(strain_name)
                
    return strains

## go over all worm files and collect the sequences of contigs into an msa
def get_contig_alignment(contig):
    implant_path = '/ebio/abt4_projects/roedel/2011-08-16_104_genomes/2012-02-08-104_alignments/alignments/variants/corrected'
    worm_list = glob.glob(implant_path+'/*.fa')

    msa = MultipleSeqAlignment([])
    for worm in worm_list:
        worm_name = worm.split('/')[-1]
        if worm_name.endswith('implant.fa'): worm_name=worm_name[:-11]
        else:  worm_name=worm_name[:-3]
        print worm_name
        with open(worm, 'r') as infile:
            found=False
            for seqrecord in SeqIO.parse(infile, 'fasta'):
                if seqrecord.id==contig:
                    msa.append(SeqRecord.SeqRecord(seq=seqrecord.seq, id=worm_name, description=''))
                    found=True
                    break

            if found == False:
                print "Contig not found in ", worm
                return None
    return msa

### get contig length
def get_contig_length():
    implant_path = '/ebio/abt4_projects/roedel/2011-08-16_104_genomes/2012-02-08-104_alignments/alignments/variants/corrected'
    worm_list = glob.glob(implant_path+'/*.fa')

    contig_list=[]
    worm = worm_list[0]
    worm_name = worm.split('/')[-1]
    if worm_name.endswith('implant.fa'): worm_name=worm_name[:-11]
    else:  worm_name=worm_name[:-3]
    print worm_name
    with open(worm, 'r') as infile:
        for seqrecord in SeqIO.parse(infile, 'fasta'):
            contig_list.append([seqrecord.id, len(seqrecord.seq)])

    return contig_list



### cogent interpretes the confidence as names, resolves name degeneracy by adding .1, .2, .3 etc. 
### biopython has a separate confidence attribute
### the following produces unique labels, includes confidence into the name and sets confidence to None
def make_names_unique(T, n):
    if T.is_terminal()==False:  #this is only nescessary for internal nodes since the leaves have unique names
        if T.name==None:
            T.name = 'edge.'+str(n)
            if T.confidence!=None:
                T.name+='_'+str(T.confidence)
                T.confidence=None
        else:
            T.name = 'edge.'+str(n)+'_'+T.name
        new_n=n+1
        for child in T.clades:
            new_n = make_names_unique(child, new_n)
    else:
        new_n = n
    return new_n
    
### get path to an alignment and a tree file name and build a tree using fast tree    
def build_ML_tree(aln_fname, tree_fname, outgroup):
    # set fasttree binary
    from socket import gethostname
    hostname = gethostname()
    if hostname=='rneher-iMac':
        fasttree_bin = 'fasttree'
    else:
        fasttree_bin = '/ebio/ag-neher/share/programs/bin/fasttree'
    
    fasttree_cmd = [fasttree_bin,"-nt", aln_fname]
    fasttree_output = subprocess.check_output(fasttree_cmd, stderr=subprocess.STDOUT)        
    #all fasttree output is returned as byte string, parse and write tree to file
    with open(tree_fname, 'w') as treefile:
        biopython_tree = Phylo.read(StringIO(fasttree_output.split('\n')[-2]), 'newick')
        n=0
        for child in biopython_tree.root.clades:
            n = make_names_unique(child,n)
        biopython_tree.root.name='root'
        Phylo.write(biopython_tree, treefile, format='newick')
        #treefile.write('#'+fasttree_output.split('\n')[-3]+'\n')

### take an alignment and a tree and code the alignment as differences into the tree labels
def add_differences_to_name(subtree, aln):
    seq_names = [seqrec.id for seqrec in aln]
    
    root_seq = np.array(aln[seq_names.index(subtree.name.split('%')[-1])])
    for child in subtree.clades:
        child_seq = np.array(aln[seq_names.index(child.name.split('%')[-1])])
        if (subtree.name!=None):
            mutations = np.where((root_seq!=child_seq)*(root_seq!='N')*(child_seq!='N'))
            label_str = "_".join([root_seq[i]+str(i)+child_seq[i] for i in mutations[0]])
            child.name = label_str+'%'+str(child.name)
        if child.is_terminal()==False:
            add_differences_to_name(child, aln)


### reconstruct ancestral sequences using cogent, 
def ancestral_sequences(aln_fname, tree_fname, anc_seq_fname, diff_tree_fname, outgroup):
    # load the tree and the alignment in the cogent format
    cogent_tree = LoadTree(tree_fname)
    cogent_aln = LoadSeqs(aln_fname)

    # set up the likelihood function (TODO figure out what the options mean
    lf = F81().makeLikelihoodFunction(cogent_tree, digits=3, space=2)
    lf.setAlignment(cogent_aln)

    # maximise the likelihood
    lf.optimise(show_progress=False, local=True)
    
    # return the likely ancestral sequences. this is an n-2 alignment 
    # of the internal nodes of the unrooted tree
    ancestors = lf.likelyAncestralSeqs()
    #return ancestors
    #convert everything back to biopython, root and ladderize
    with open(tree_fname, 'r') as infile:
        biopython_tree = Phylo.read(infile, 'newick')

    biopython_aln = AlignIO.read(StringIO(cogent_aln.toFasta()), 'fasta')
    biopython_aln.extend(AlignIO.read(StringIO(ancestors.toFasta()), 'fasta'))
    biopython_tree.root_with_outgroup(outgroup)
    biopython_tree.ladderize()
    seq_name_list = [seqrec.id for seqrec in biopython_aln]

    #pull out root sequence and write to file
    R = biopython_tree.root
    if len(R.clades)==2:
        if R.clades[0].name==outgroup:
            root_of_tree = R.clades[1]
        elif R.clades[1]==outgroup:
            root_of_tree = R.clades[0]
        else:
            raise "outgroup is not a child of root"
    else:
        raise "Root node is not bifurcating"

    #return biopython_tree, biopython_aln, seq_name_list
    root_sequence = biopython_aln[seq_name_list.index(root_of_tree.name)]
    with open(anc_seq_fname, 'w') as outfile:
        SeqIO.write(root_sequence, outfile, format='fasta')

    #encode leafs as differences from the root
    add_differences_to_name(root_of_tree, biopython_aln)
    with open(diff_tree_fname, 'w') as outfile:
        Phylo.write(biopython_tree, outfile, 'newick')

    return biopython_tree, biopython_aln, root_sequence

### recursively collect substitution and number of downstream leaves
def calc_SFS_subtree(subtree, substitutions, branch_length):
    subs = subtree.name.split('%')[0].split('_')
    n_leaves = len(subtree.get_terminals())
    branch_length.append((n_leaves, subtree.branch_length))
    for sub in subs:
        if sub!='':
            if sub not in substitutions:
                substitutions[sub]=[n_leaves]
            else:
                substitutions[sub].append(n_leaves)
    if subtree.is_terminal()==False:
        for child in subtree.clades:
            calc_SFS_subtree(child, substitutions, branch_length)


### Parse the difference tree file and enumerate all substitutions
def calc_SFS(tree_fname, outgroup):
    with open(tree_fname, 'r') as infile:
        T=Phylo.read(infile, 'newick')

    substitutions = {}
    branch_length = []
    
    R = T.root
    if len(R.clades)==2:
        if R.clades[0].name==outgroup:
            root_of_tree = R.clades[1]
        elif R.clades[1]==outgroup:
            root_of_tree = R.clades[0]
        else:
            raise "outgroup is not a child of root"
    else:
        raise "Root node is not bifurcating"

    calc_SFS_subtree(root_of_tree, substitutions, branch_length)
    return substitutions, branch_length





