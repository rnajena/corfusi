#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
import pprint
from BCBio.GFF import GFFExaminer
import os
import sys, getopt

#######################################################################################
### ORF correction of hybrid genome assemblies using short-read assembly annotation ###
#######################################################################################

fasta_file, gff_file, outputdir, prefix, t = '', '', './', '', 50


try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:g:t:o:p:",["fasta=","gff=", "threshold=", "outputdir=", "prefix="])
except getopt.GetoptError:
    print('orf_correction.py -f <fasta file> -g <gff file> -t <threshold for up-/downstream region> -p <prefix of output> -o <output directory>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('orf_correction.py -f <fasta file> -g <gff file> -t <threshold for up-/downstream region> -p <prefix of output> -o <output directory>')
        sys.exit()
    elif opt in ("-f", "--fasta"):
        fasta_file = arg
    elif opt in ("-g", "--gff"):
        gff_file = arg
    elif opt in ("-o", "--outputdir"):
        outputdir = arg
    elif opt in ("-p", "--prefix"):
        prefix = arg
    elif opt in ("-t", "--threshold"):
        t = int(arg)

if not outputdir[-1] == '/': outputdir += '/'
if prefix.split('.')[-1] == 'fa' or prefix.split('.')[-1] == 'fasta':
    ''.join(prefix.split('.')[:-1])

os.system('mkdir ' + outputdir)
os.system('mkdir '+ outputdir + 'blastn')


############### load fasta ###############
hybrid_fasta = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))


############### load gff ###############
in_handle = open(gff_file)
limit_info = dict(gff_source=["Prodigal:002006"])
gff = []
for rec in GFF.parse(in_handle, limit_info=limit_info):
    if len(rec.features) > 0:
        gff.append(rec)
in_handle.close()


############### makebladtdb ###############
os.system('makeblastdb -in ' +  fasta_file + ' -parse_seqids -dbtype nucl -out ' + outputdir + 'hybrid_blastdb')


############### find candidates ###############
candidates = []

for node in gff:
    # id = node.id
    seq = node.seq
    for feature in node.features:
        start = feature.location.start
        end = feature.location.end
        # strand = feature.location.strand
        feature_seq = seq[start:end+1]
        feature_id = feature.id
        
        ### create query ###
        os.system('touch query.fasta')
        os.system('echo ">' + feature_id + '\n' + str(feature_seq) + '"  >> query.fasta')

        
        ### blastn ###
        os.system('blastn -task blastn -outfmt 6 -max_target_seqs 1 -culling_limit 1 -query query.fasta -db ' + outputdir + 'hybrid_blastdb -out ' + outputdir + 'blastn/' + feature_id + '_results.out 2>/dev/null')
        os.system('rm query.fasta')
        
        results = open(outputdir + 'blastn/' + feature_id + '_results.out', 'r').readline().split('\t')

        if len(results) > 1:
           pident = float(results[2])
        else: pident = 0
        
        if pident > 90:
           length = int(results[3])
           qcov = length / (int(results[7]) - int(results[6]) + 1)
           scov = length /(int(results[9]) - int(results[8]) + 1)

            ### filter candidates ###
           if qcov != 1.0 or int(results[4]) != 0 or int(results[5]) != 0:
               candidates.append([node, feature_id, feature.location, pident, qcov, scov, int(results[4]), int(results[5])])

os.system('rm ' + outputdir + 'blastn/*')


############### filtering with blast upstream and downstream region ###############
up_down_all = []

for elem in candidates:
    ### extract up and downstream region ###
    start = int(elem[2].start)
    end = int(elem[2].end)
    upstream = elem[0].seq[start-t:start]
    downstream = elem[0].seq[end+1:end+t+1]
    
    ### check length of up and downstream region ###
    if len(upstream) == t and len(downstream) == t:
        ### create query ###
        os.system('touch query.fasta')
        os.system('echo ">upstream\n' + str(upstream) + '\n>downstream\n' + str(downstream) + '" >> query.fasta')

        ### blastn ###
        os.system('blastn -task blastn -max_target_seqs 1 -outfmt 6 -culling_limit 1 -query query.fasta -db ' + outputdir + 'hybrid_blastdb -out ' + outputdir + 'blastn/' + elem[1] + '_results.out 2>/dev/null')
        os.system('rm query.fasta')

        results = open(outputdir + 'blastn/' + elem[1] + '_results.out', 'r').readlines()
        
        up_down_pair = [results[0].split('\t')]
        for line in results:
            if line[0:10] == 'downstream':
                up_down_pair.append(line.split('\t'))
                break
        
        ### continue if up and downstream region matched ###
        if len(up_down_pair) == 2 and up_down_pair[0][1] == up_down_pair[1][1] and int(up_down_pair[0][3]) == t and int(up_down_pair[1][3]) == t:
            k = up_down_pair[0][0:2] + list(map(float, up_down_pair[0][2:]))
            l = up_down_pair[1][0:2] + list(map(float, up_down_pair[1][2:]))
            up_down_all.append((k, l, elem))
            
        up_down_pair = []

### sort blast results by start position ###
up_down_all.sort(key=lambda x: x[0][8])


############### update assembly sequence ###############
count = 0
log = open(outputdir + prefix + '_log.tsv', 'w')
log.write('short-read assembly prokka id\tstart\tend\told sequence\tnew sequence\n')

for elem in up_down_all:
    h_start = int(elem[0][9]) + count
    h_end = int(elem[1][8]) + count

    start = elem[2][2].start
    end = elem[2][2].end

    h_len = abs(h_end - h_start -1)
    sr_gene_len = end - start +1

    sr_gene = elem[2][0].seq[start:end+1]

    ### if downstream < upstream: change variable and generate reverse complement ###
    if h_start > h_end:
        h_start, h_end = h_end, h_start
        sr_gene = sr_gene.reverse_complement()

    ### filtering by length (20% longer or shorter than gene allowed) ###
    if abs(h_len - sr_gene_len) < sr_gene_len * 0.2:
        ### insert short-read gene in hybrid sequence ###
        hybrid_fasta[elem[0][1]].seq = hybrid_fasta[elem[0][1]].seq[:h_start+1] + sr_gene + hybrid_fasta[elem[0][1]].seq[h_end:]

        ### update index count ###
        if sr_gene_len != h_len: count += sr_gene_len - h_len
        
        ### write logfile ###
        log.write(elem[2][1] + '\t' + str(h_start+1) + '\t' + str(h_end-1) + '\t' + str(hybrid_fasta[elem[0][1]].seq[h_start+1:h_end]) + '\t' + str(sr_gene) + '\n')

log.close()

os.system('rm -r ' + outputdir + 'blastn/')
os.system('rm ' + outputdir + 'hybrid_blastdb*')


############## save new assembly ###############
with open(outputdir + prefix + '.fasta', 'w') as handle:
    SeqIO.write(hybrid_fasta.values(), handle, 'fasta')
