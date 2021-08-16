# takes as input two paired Illumina NGS FASTQ files
# joins two related paired-end sequencing files (.fastq) into a single .txt file 
# in the following format: R1+\t+R2+\n
# reads are quality filtered
# output .txt file can be entered as an input for spc9_analyzer.py


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

# enter .fastq file names for the input files:
fastq_1 = open("Fig3F_spc9_mutations_dinB_ecoli_R1.fastq", "r")
fastq_2 = open("Fig3F_spc9_mutations_dinB_ecoli_R2.fastq", "r")

# enter .txt file name for the output file:
output = open("Spc9MutSeq_S1_L001_R1_R2_combined_quality.txt", "w")

# parse .fastq files
reads_1 = SeqIO.parse(fastq_1, 'fastq')
reads_2 = SeqIO.parse(fastq_2, 'fastq')


# isQuality takes as input a FASTQ record. Returns true if the sequences is
# high quality or false otherwise. Generaly, a "phred" quality score greater than 10
# is considered good (90% accuracy)
def isQuality(record):
    # list of qualities for each DNA letter in the record
    # JM - edited to look skip first and last 5 bp of each read
    qualities = record.letter_annotations["phred_quality"]
    qualities_short = qualities[5:-5]

    # if there is at least one quality smaller than 10
    # this record will be labeled as poor quality
    for quality in qualities_short:
        if quality < 10:
            return False

    # if all qualities are greater than 10, this is a high-quality record
    return True


timer = 0
good = 0

for record in reads_1:
    r1 = record
    r2 = reads_2.next()

    if isQuality(r1) and isQuality(r2):
        output.write(str(r1.seq) + '\t' + str(r2.seq) + '\n')
        good += 1

    timer += 1
    if (timer + 0.0) % 100000 == 0:
        print timer
        print "Q%: " + str(float(good) / timer)

# clean up
fastq_1.close()
fastq_2.close()
output.close()



