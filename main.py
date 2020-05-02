import argparse, alignment as al

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', required=True, help='fasta file')
parser.add_argument('-o', '--output', required=True, help='output file prefix')
parser.add_argument('-al', '--alignment', required=True, help='type of alignment: write "fasta" if yout want fasta alignmen, or "blast" and "n-w" in their cases')

args= parser.parse_args()

FASTA_input = open(args.input, 'r')
Sequences = FASTA_input.read()
FASTA_input.close()

Seq = []
Seq_inner = []
Seq = Sequences.split('>')
for i in Seq:
    Seq_inner.append(i.split('\n'))
Seq_inner.remove(Seq_inner[0])
sequences = []
for i in Seq_inner:
    s = ''
    for j in range(1,len(i)):
        s += i[j]
    sequences.append(s)
if len(sequences[0]) > len(sequences[1]):
    seq = sequences[1]
    seqi = sequences[0]
else:
    seq = sequences[0]
    seqi = sequences[1]

#Нидлман-Вунш
if args.alignment == 'n-w':
    output_name = args.output + '.Needleman-Wunch'
    Out = open(output_name, 'w')
    STRAX = al.ALIGN_NW(seq,seqi)
    if(len(STRAX[0]) > 60):
        for i in range(len(STRAX)):
            for k in range(int(len(STRAX[i]) / 60)):
                for j in range(60):
                    if (k * 60 + j < len(str(STRAX[i]))):
                        Out.write(str(STRAX[i][k*60 + j]))
                Out.write('\n')
            Out.write('\n')
    else:
        for i in range(len(STRAX)):
            Out.write(str(STRAX[i]))
            Out.write('\n')
    Out.close()

#Фаста
if args.alignment == 'fasta':
    output_name = args.output + '.Fasta'
    Out = open(output_name, 'w')
    STRAX = al.FASTA(seqi,seq)
    if(len(STRAX[0]) > 60):
        for i in range(len(STRAX)):
            for k in range(int(len(STRAX[i]) / 60)):
                for j in range(60):
                    if (k * 60 + j < len(str(STRAX[i]))):
                        Out.write(str(STRAX[i][k*60 + j]))
                Out.write('\n')
            Out.write('\n')
    else:
        for i in range(len(STRAX)):
            Out.write(str(STRAX[i]))
            Out.write('\n')
    Out.close()

#БЛАСТ
if args.alignment == 'blast':
    output_name = args.output + '.BLAST'
    Out = open(output_name, 'w')
    STRAX = al.BLAST(seqi, seq, 5)
    if(len(STRAX[0]) > 60):
        for i in range(len(STRAX)):
            for k in range(int(len(STRAX[i]) / 60)):
                for j in range(60):
                    if (k * 60 + j < len(str(STRAX[i]))):
                        Out.write(str(STRAX[i][k*60 + j]))
                Out.write('\n')
            Out.write('\n')
    else:
        for i in range(len(STRAX)):
            Out.write(str(STRAX[i]))
            Out.write('\n')
    Out.close()

#Левенштейн
output_name = args.output + '.Levin'
Out = open(output_name, 'w')
STRAX = al.distance(seqi, seq)
Out.write(str(STRAX))
Out.close()

