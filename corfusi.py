#!/usr/bin/python3

from Bio import SeqIO
from BCBio import GFF
import os
import sys, getopt
import time

#######################################################################################
############ cORFusi - correction of ORFs utilizing short-read information ############
#######################################################################################

fasta_file, gff_file, outputdir, prefix, t = '', '', './', '', 50

help = """cORFusi - correction of ORFs utilizing short-read information
usage: python corfusi.py -f ASSEMBLY -g ANNOTATION -t INT -p NAME -o PATH

Parameters:
-h          show this message and exit
-f assembly, --fasta assembly
            assembly file (.fasta)
-g annotation, -gff annotation
            annotation file (.gff)
-t int, --threshold int
            length threshold determing up- & downstream region
-p name, --prefix name
            prefix for output files
-o path, --outdir path
            path to output folder
"""

try:
    opts, args = getopt.getopt(sys.argv[1:],"hf:g:t:o:p:",["fasta=","gff=", "threshold=", "outdir=", "prefix="])
except getopt.GetoptError:
    print(help)
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print(help)
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

if not os.path.isdir(outputdir): os.system('mkdir ' + outputdir)
os.system('mkdir '+ outputdir + 'blastn')

# time_start = time.process_time()


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
os.system('makeblastdb -in ' +  fasta_file + ' -parse_seqids -dbtype nucl -out ' + outputdir + 'hybrid_blastdb >/dev/null')


############### find candidates ###############
candidates = []
yes, no = 0,0

for node in gff:
    seq = node.seq
    for feature in node.features:
        start = feature.location.start
        end = feature.location.end
        feature_seq = seq[start:end]
        feature_id = feature.id
        
        ### create query ##
        os.system('touch query.fasta')
        os.system('echo ">' + feature_id + '\n' + str(feature_seq) + '"  >> query.fasta')
        
        ### blastn ###
        os.system('blastn -task blastn -outfmt 6 -max_target_seqs 1 -culling_limit 1 -query query.fasta -db ' + outputdir + 'hybrid_blastdb -out ' + outputdir + 'blastn/' + feature_id + '_results.out 2>/dev/null')
        os.system('rm query.fasta')
        
        results = open(outputdir + 'blastn/' + feature_id + '_results.out', 'r').readline().split('\t')

        if len(results) > 1:
           pident = float(results[2])
           length = int(results[3])
        else: pident = 0

        if pident > 90 and pident < 100 and abs(length - len(feature_seq)) <= 0.2 * len(feature_seq):
            # qcov = length / (abs(int(results[7]) - int(results[6])) + 1)
            # scov = length /(abs(int(results[9]) - int(results[8])) + 1)
            candidates.append([node, feature_id, feature.location]) #, pident, qcov, scov, int(results[4]), int(results[5])])

        if pident == 100: yes += 1
        else: no += 1

# os.system('rm ' + outputdir + 'blastn/*')


############### filtering with blast upstream and downstream region ###############
up_down_all = []

for elem in candidates:
    ### extract up and downstream region ###
    start = int(elem[2].start)
    end = int(elem[2].end)
    upstream = elem[0].seq[start-t:start]
    downstream = elem[0].seq[end:end+t]

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
log.write('annotation_ID\tstart\tend\told sequence\tnew sequence\n')

for elem in up_down_all:
    start = elem[2][2].start
    end = elem[2][2].end

    sr_gene = elem[2][0].seq[start:end]

    up_s = int(elem[0][8])
    up_e = int(elem[0][9])
    down_s = int(elem[1][8])
    down_e = int(elem[1][9])

    if up_s < down_s: ### upstream < downstream: define variables ###
        if up_s < up_e: h_start = up_e + count
        else: h_start = up_s + count
        if down_s < down_e: h_end = down_s + count
        else: h_end = down_e + count
    else: ### downstream < upstream: define variables and generate reverse complement ###
        if down_s < down_e: h_start = down_e + count
        else: h_start = down_s + count
        if up_s < up_e: h_end = up_s + count
        else: h_end = up_e + count
        sr_gene = sr_gene.reverse_complement()

    old = hybrid_fasta[elem[0][1]].seq[h_start:h_end]

    h_len = len(old)
    sr_gene_len = len(sr_gene)

    ### filtering by length (20% longer or shorter than gene allowed) ###
    if abs(h_len - sr_gene_len) < sr_gene_len * 0.2:
        ### insert short-read gene in hybrid sequence ###
        hybrid_fasta[elem[0][1]].seq = hybrid_fasta[elem[0][1]].seq[:h_start] + sr_gene + hybrid_fasta[elem[0][1]].seq[h_end:]

        ### update index count ###
        if sr_gene_len != h_len: count += sr_gene_len - h_len

        ### write logfile ###
        log.write(elem[2][1] + '\t' + str(h_start+1) + '\t' + str(h_end+(sr_gene_len - h_len)) + '\t' + str(old) + '\t' + str(sr_gene) + '\n')

log.close()

os.system('rm -r ' + outputdir + 'blastn/')
os.system('rm ' + outputdir + 'hybrid_blastdb*')


############## save new assembly ###############
with open(outputdir + prefix + '.fasta', 'w') as handle:
    SeqIO.write(hybrid_fasta.values(), handle, 'fasta')

print("ORF correction process complete!\n" + "Assembly " + prefix + ".fasta saved in " + outputdir)
print("Number 100% matches: ", yes, "\nNumber of incorrect matches: ", no)
# print("Time: ", time.process_time() - time_start)
print("Thank you for using cORFusi!")