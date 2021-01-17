#!/usr/bin/python3

from Bio import SeqIO
from BCBio import GFF
import pprint
from BCBio.GFF import GFFExaminer
import os
import sys

#############################################################################
### ORF correction of bacterial genomes using short and hybrid assemblies ###
#############################################################################

############### examining gff ###############
# examiner = GFFExaminer()
# pprint.pprint(examiner.available_limits(in_handle))
# in_handle.close()

fasta_file = sys.argv[1]
gff_file = sys.argv[2]


############### load fasta ############### 
seq = []
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    seq.append(seq_record)


############### load gff ###############
in_handle = open(gff_file)
gff = []
for rec in GFF.parse(in_handle):
    if len(rec.features) > 0:
        gff.append(rec)
in_handle.close()


############### makebladtdb ###############
# hybrid = sys.argv[3]
# os.system('makeblastdb -in' +  hybrid + ' -parse_seqids -dbtype nucl -out' + 'blastdb')


############## ... ###############
outfile = open('outfile', 'w')
outfile.write('id\tipdent\tqcov\tscov')
for node in gff:
    id = node.id
    seq = node.seq
    for feature in node.features:
        start = feature.location.start
        end = feature.location.end
        strand = feature.location.strand
        feature_seq = seq[start:end]
        feature_id = feature.id
        
        os.system('touch query.fasta')
        os.system('echo ">' + feature_id + '\n' + str(feature_seq) + '" >> query.fasta')
        
        os.system('blastn -task blastn -max_target_seqs 1 -outfmt 6 -query query.fasta -db blastdb -out blast_results/' + feature_id + '_results.out')
        
        results = open('blast_results/' + feature_id + '_results.out', 'r').readline().split('\t')
        os.system('touch blast_results/' + feature_id + '_results.out')
        
        pident = float(results[2])
        if pident > 90:
            length = int(results[3])
            qcov = length / (int(results[7]) - int(results[6]) + 1)
            scov = length /(int(results[9]) - int(results[8]) + 1)
            outfile.write(feature_id + '\t' + str(pident) + '\t' + str(qcov) + '\t' + str(scov))
        
        os.system('rm query.fasta')

outfile.close()



# h = gff[0]
# f = h.features
# print(f[0].id)
# print(f[0])
# print(f[0].location.start) # start end strand
# print(h.seq[f[0].location.start:f[0].location.end])


# results = open('C:/Users/sandr/Dropbox/Master/paper/11DD0261/results.out', 'r').readline().split('\t')
# pident = float(results[2])
# length = int(results[3])
# qcov = length / (int(results[7]) - int(results[6]) + 1)
# scov = length /(int(results[9]) - int(results[8]) + 1) 
# print(results)
# print(pident, qcov, scov)