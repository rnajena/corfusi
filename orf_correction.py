#!/usr/bin/python3

from Bio import SeqIO
from BCBio import GFF
import pprint
from BCBio.GFF import GFFExaminer
import os
import sys, getopt

fasta_file = ''
gff_file = ''
t = 100

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:g:t:",["fasta=","gff=", "threshold="])
except getopt.GetoptError:
    print('orf_correction.py -f <fastafile> -g <gfffile> -t <threshold for up-/downstream region>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print('orf_correction.py -f <fastafile> -g <gfffile> -t <threshold for up-/downstream region>')
        sys.exit()
    elif opt in ("-f", "--fasta"):
        fasta_file = arg
    elif opt in ("-g", "--gff"):
        gff_file = arg
    elif opt in ("-t", "--threshold"):
        t = int(arg)


#############################################################################
### ORF correction of bacterial genomes using short and hybrid assemblies ###
#############################################################################

############### examining gff ###############
# examiner = GFFExaminer()
# pprint.pprint(examiner.available_limits(in_handle))
# in_handle.close()


############## load fasta ############### 
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
# hybrid = sys.argv[3]
# os.system('makeblastdb -in' +  hybrid + ' -parse_seqids -dbtype nucl -out' + 'blastdb')


############## find candidates ###############
outfile = open('outfile', 'w')
outfile.write('contig\tid\tid_pos\tipdent\tqcov\tscov\tmm\tgap\n')
candidates = []
for node in gff:
    id = node.id
    seq = node.seq
    for feature in node.features:
        start = feature.location.start
        end = feature.location.end
        strand = feature.location.strand
        feature_seq = seq[start:end+1]
        feature_id = feature.id
        
        ### create query ###
        os.system('touch query.fasta')
        os.system('echo ">' + feature_id + '\n' + str(feature_seq) + '" >> query.fasta')
        
        ### blastn ###
        os.system('blastn -task blastn -outfmt 6 -max_target_seqs 1 -query query.fasta -db blastdb -out blastn/' + feature_id + '_results.out')
        os.system('rm query.fasta')
        
        results = open('blastn/' + feature_id + '_results.out', 'r').readline().split('\t')
        
        if len(results) > 1:
            pident = float(results[2])
        else: pident = 0
        
        if pident > 90:
            length = int(results[3])
            qcov = length / (int(results[7]) - int(results[6]) + 1)
            scov = length /(int(results[9]) - int(results[8]) + 1)

            #### filter candidates ###
            if qcov != 1.0 or int(results[4]) != 0 or int(results[5]) != 0:
                outfile.write(id + '\t' + feature_id + '\t' + str(feature.location) + '\t' + str(pident) + '\t' + str(qcov) + '\t' + str(scov) + '\t' + str(results[4]) + '\t' + str(results[5]) + '\n')
                candidates.append([node, feature_id, feature.location, pident, qcov, scov, int(results[4]), int(results[5])])

outfile.close()


############## update assembly sequence ###############
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
        os.system('blastn -task blastn -max_target_seqs 1 -outfmt 6 -query query.fasta -db blastdb -out blastn2/' + elem[1] + '_results.out')
        os.system('rm query.fasta')

        results = open('blastn2/' + elem[1] + '_results.out', 'r').readlines()
        
        up_down_pair = [results[0].split('\t')]
        for line in results:
            if line[0:10] == 'downstream':
                up_down_pair.append(line.split('\t'))
                break
        
        ### continue if up and downstream region matched ###
        if len(up_down_pair) == 2 and up_down_pair[0][1] == up_down_pair[1][1] and int(up_down_pair[0][3]) == t and int(up_down_pair[1][3]) == t:
            up_down_all.append((up_down_pair[0], up_down_pair[1], elem))
            
        up_down_pair = []

### sort blast results by start position ###
up_down_all.sort(key=lambda x: x[0][8])

count = 0
for elem in up_down_all:
    h_start = int(elem[0][9]) + count
    h_end = int(elem[1][8]) + count

    start = elem[3][2].start
    end = elem[3][2].end

    h_len = h_end - h_start -1
    sr_gene_len = end - start +1

    ### filtering by length (20% longer than gene allowed) ###
    if h_len <= sr_gene_len * 1.2 and h_len >= sr_gene_len * 0.8:
        sr_gene = elem[0].seq[start:end+1]

        ### insert short-read gene in hybrid sequence ###
        hybrid_fasta[elem[0][1]].seq = hybrid_fasta[elem[0][1]].seq[:h_start] + sr_gene + hybrid_fasta[elem[0][1]].seq[h_end+1:]

        ### update index count ###
        if sr_gene_len != h_len:
            count += sr_gene_len - h_len


############## save new assembly ###############
with open('11DD0261_new.fasta', 'w') as handle:
    SeqIO.write(hybrid_fasta.values(), handle, 'fasta')
