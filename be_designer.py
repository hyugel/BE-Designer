import os, sys, argparse, re
import mich

from subprocess import Popen, PIPE

placeholder = """>homo sapiens FANCM, exon 2
GGTCTACACAAGCTTCCACCAGGAAGGAAATA
TGGTGCAGTAAGAGAGTGCTTTTTCTTACACC
TCAGGTCATGGTAAATGACCTTTCTAGAGGAG
CTTGTCCCGCTGCTGAAATAAAGTGTTTAGTT
ATTGATGAAGCTCATAAAGCTCTCGGAAACTA
TGCTTATTGCCAG"""

amino_sym = {"TGG": "W", "GGG": "G", "GTT": "V", "CGG": "R", "ACT": "T", "TGA": "X", "CTA": "L", "TCC": "S","GAA": "E", "CCA": "P", "GAC": "D", "ACC": "T", "TTT": "F", "CTC": "L", "GCT": "A", "CCC": "P", "TCG": "S", "CAT": "H", "GTC": "V", "CGA": "R", "CAG": "Q", "ATA": "I", "AAG": "K", "CCG": "P", "GGA": "G", "AGC": "S", "TAT": "Y", "CTG": "L", "ACG": "T", "GAG": "E", "GCT": "A", "TGC": "C", "TGT": "C", "AGG": "R", "ATG": "M", "TTA": "L", "GCA": "A", "AAT": "N", "GTA": "V", "GGT": "G", "AGA": "R", "CGC": "R", "ATC": "I", "TAC": "Y", "TTG": "L", "ACA": "T", "GCG": "A", "CTT": "L", "ATT": "I", "CGT": "R", "CAC": "H", "TCA": "S", "CCT": "P", "TAA": "X", "GAT": "D", "GTG": "V", "AAA": "K", "AAC": "N", "GGC": "G", "TTC": "F", "CAA": "Q", "AGT": "S", "TAG": "X", "TCT": "S", "GCC": "A"}

def reverse_complement(s):
    return str(s).translate(s.maketrans('ATGC','TACG'))[::-1]

class BEDesigner:

    def __init__(self, args):

        self.query_seq = args.query_seq.upper()
        self.seed_len = args.seed_len
        self.pam_seq = args.pam_seq.upper()
        self.window_start = args.window_st
        self.window_ed = args.window_ed
        self.ref_base = args.ref_base.upper()
        self.alt_base = args.modified_base.upper()
        self.isreversed = False
        self.bulge_dna = 0
        self.bulge_rna = 0
        self.mismatch = args.mismatch
        self.ref_dir = args.ref_dir
        self.output = args.output

    def find_targets(self, targetonly = False):

        seqs = []
        if self.query_seq[0] == '>':
            lines = self.query_seq.strip().replace('\n','').replace('\t','').upper()
            for line in lines:
                if line.strip() == '':
                    continue
                if line[0] == '>':
                    seqs.append([line[1:], ''])
                else:
                    seqs[-1][1] += line.upper()
        else:
            seqs = [ ['Untitled', self.query_seq.upper()] ]
        targets = []
        ref_amino_list = ['','','']
        pamlen = len(self.pam_seq)

        for title, seq in seqs:
            for i in range(3):
                for x in range((len(seq)-i)//3):
                    ref_amino_list[i] += amino_sym[seq[i:][3*x:3*x+3]]
            if not targetonly:
                targets.append( [ title+';316;'+seq+';316;'+';316;'.join(ref_amino_list), [] ] )
            for target in mich.target_yield(seq, 60, self.seed_len, self.pam_seq, self.isreversed, []):
                mutated_pos = [[-1], [-1]]
                mutated_aa = [[0, 0, ''], [0, 0, ''], [0, 0, '']]
                pos = int(target[3])
                if target[1] == '+':
                    window_start_abs = int(target[3]) + self.seed_len - self.window_ed
                    window_end_abs = int(target[3]) + self.seed_len - self.window_start
                    if self.isreversed:
                        window_start_abs += pamlen
                        window_end_abs += pamlen
                    window_seq = seq[window_start_abs: window_end_abs + 1]
                    translation_st = window_start_abs - 2
                    translation_ed = window_end_abs + 2
                    if self.isreversed:
                        mutated_pos[0] += [pos+pamlen, pos+len(target[0]), window_start_abs, window_end_abs]
                        mutated_pos[1] += [pos, pos+pamlen]
                    else:
                        mutated_pos[0] += [pos, pos+len(target[0])-pamlen, window_start_abs, window_end_abs]
                        mutated_pos[1] += [pos+len(target[0])-pamlen, pos+len(target[0])]
                    if translation_st < 0:
                        translateion_seq = 'X'*abs(translation_st) + seq[0: translation_ed + 1]
                    else:
                        translation_seq = seq[translation_st: translation_ed + 1]
                else:
                    window_start_abs = int(target[3]) + pamlen -1 + self.window_start
                    window_end_abs = int(target[3]) + pamlen - 1 + self.window_ed
                    if self.isreversed:
                        window_start_abs -= pamlen
                        window_end_abs -= pamlen
                    window_seq = reverse_complement(seq[window_start_abs: window_end_abs +1])
                    translation_st = window_start_abs - 2
                    translation_ed = window_end_abs + 2
                    if self.isreversed:
                        mutated_pos[0] += [pos, pos+ len(target[0]) -pamlen, window_start_abs, window_end_abs]
                        mutated_pos[1] += [pos + self.seed_len, pos+ self.seed_len - pamlen]
                    else:
                        mutated_pos[0] += [pos+pamlen, pos+len(target[0]), window_start_abs, window_end_abs]
                        mutated_pos[1] += [pos, pos+pamlen]
                    if translation_ed >= len(seq):
                        translation_seq = seq[translation_st:] + 'X'*(translation_ed - len(seq) + 1)
                    else:
                        translation_seq = seq[translation_st: translation_ed + 1]
                mutated_translation_seq = translation_seq[:2]
                for i in range(2, len(translation_seq) - 2):
                    if translation_seq[i] == self.ref_base and target[1] == '+':
                        mutated_translation_seq += self.alt_base
                    elif translation_seq[i] == reverse_complement(self.ref_base) and target[1] == '-':
                        mutated_translation_seq += reverse_complement(self.alt_base)
                    else:
                        mutated_translation_seq += translation_seq[i]
                mutated_translation_seq += translation_seq[-2:]
                amino_seq = [['', '', '',0, ''], ['', '', '', 0, ''], ['', '', '', 0, '']]
                for i in range(3):
                    if self.isreversed:
                        if target[1] == '+':
                            codon = (target[3] + pamlen + self.window_start - 1 - i) % 3
                        else:
                            codon = (target[3] + self.seed_len - self.window_ed - i) % 3
                    else:
                        if target[1] == '+':
                            codon = (target[3] + self.seed_len - self.window_ed - i) % 3
                        elif target[1] == '-':
                            codon = (target[3] + pamlen + self.window_start - 1 - i ) % 3
                    if codon == 2: codon_add = 0
                    elif codon == 1: codon_add = 1
                    elif codon == 0: codon_add = 2
                    cut_mutated_translation_seq = mutated_translation_seq[codon_add:]
                    amino_seq[i][2] = cut_mutated_translation_seq[:len(cut_mutated_translation_seq) - len(cut_mutated_translation_seq) % 3]
                    amino_seq[i][3] = codon
                    cut_translation_seq = translation_seq[codon_add:]
                    amino_seq[i][4] = cut_translation_seq[:len(cut_translation_seq) - len(cut_translation_seq) % 3]
                    mutated_aa[i][0] = window_start_abs - codon
                    c = 0
                    while c < len(cut_mutated_translation_seq) + codon_add:
                        ref_amino = cut_translation_seq[c:c+3]
                        mut_amino = cut_mutated_translation_seq[c:c+3]
                        if len(ref_amino) != 3: break
                        elif c == 0 and ref_amino.find('X') != -1:
                            amino_seq[i][0] += 'St'
                            amino_seq[i][1] += 'St'
                        elif c != 0 and ref_amino.find('X') != -1:
                            amino_seq[i][0] += 'Ed'
                            amino_seq[i][1] += 'Ed'
                        else:
                            amino_seq[i][0] += amino_sym[ref_amino]
                            amino_seq[i][1] += amino_sym[mut_amino]
                            if amino_sym[mut_amino] == 'X':
                                mutated_aa[i][2] += 'Ter'
                            else:
                                mutated_aa[i][2] += '_' + amino_sym[mut_amino] + '_'
                        c += 3
                        mutated_aa[i][1] = window_start_abs - codon + len(mutated_aa[i][2])
                if not self.ref_base in window_seq:
                    continue
                if targetonly:
                    if self.isreversed:
                        found_target = target[2][len(self.pam_seq):]
                    else:
                        found_target = target[2][:-len(self.pam_seq)]
                    if not found_target in targets:
                        targets.append(found_target)
                else:
                    if self.isreversed:
                        seedseq = target[2][len(self.pam_seq):]
                    else:
                        seedseq = target[2][:-len(self.pam_seq)]
                    gc = ((seedseq.count('G') +seedseq.count('C'))* 100.0) / len(seedseq)
                    targets[-1][1].append([target[2],  # Found sequence
                                           window_seq,
                                           target[3] + 1,  # Position
                                           target[1],  # Direction of target
                                           gc,
                                           amino_seq,
                                           mutated_pos,
                                           mutated_aa])

        self.targets = targets

    def run_cas_offinder(self):

        s = [self.ref_dir]
        self.targets_dict = {}
        if self.isreversed:
            s.append(self.pam_seq + (self.seed_len * 'N') + ' '+ str(self.bulge_dna) + ' ' + str(self.bulge_rna))
        else:
            s.append((self.seed_len * 'N') + self.pam_seq + ' ' + str(self.bulge_dna) + ' ' + str(self.bulge_rna))
        for line in self.targets[0][1]:
            for c in line[0]:
                if not c in 'ACGTRYSWKMBDHYN':
                    raise Exception
            if not 15 <= len(line[0]) <= 30:
                raise Exception
            if self.isreversed:
                s.append((len(self.pam_seq) * 'N') + line[0][len(self.pam_seq):] + ' ' + str(self.mismatch))
                self.targets_dict[line[0][len(self.pam_seq):]] = {'mis': {i: [] for i in range(self.mismatch + 1)}}

            else:
                s.append(line[0][:-len(self.pam_seq)] + (len(self.pam_seq) * 'N') + ' ' + str(self.mismatch))
                self.targets_dict[line[0][:-len(self.pam_seq)]] = {'mis': {i: [] for i in range(self.mismatch + 1)}}

        with open('cas_offinder_input.txt', 'w') as f:
            f.write('\n'.join(s))
        p = Popen(('./cas-offinder-bulge', 'cas_offinder_input.txt', 'G', 'cas_offinder_out.txt'), stdout = PIPE)
        p.communicate()
        ret = p.wait()

        with open('cas_offinder_out.txt') as f:
            for line in f:
                if line[0] == '#': 
                    continue
                line_s = line.strip().split('\t')
                if len(line_s) < 8:
                    continue
                if self.isreversed:
                    seq = line_s[1][len(self.pam_seq):]
                else:
                    seq = line_s[1][:-len(self.pam_seq)]
                self.targets_dict[seq]['mis'][int(line_s[6])].append('\t'.join(line_s[1:7]))

        with open(self.output, 'w') as fw:
            fw.write('\t'.join(["CRISPR Target (5' to 3')", "Position", "Direction", "GC contents (%, w/d PAM)"]+[str(i) for i in range(self.mismatch + 1)]+['codon 0 pattern','','','','codon 1 pattern','','','','codon 2 pattern','','',''])+'\n')
            for line in self.targets[0][1]:
                if self.isreversed:
                    seq = line[0][len(self.pam_seq):]
                else:
                    seq = line[0][:-len(self.pam_seq)]
                mut_info = []
                for i in line[5]:
                    mut_info += [i[4], i[2], i[0], i[1]]
                fw.write('\t'.join([line[0], str(line[2]), line[3], str(line[4])] + [str(len(self.targets_dict[seq]['mis'][i])) for i in range(self.mismatch+1)] + mut_info)+ '\n')


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("query_seq", type = str, help = "Insert any sequence(s) where you want to search for CRISPR base editing targets")
    parser.add_argument("-s", "--seed_len", type = int, help = "length of target without PAM (default = 20)", default = 20)
    parser.add_argument("--pam_seq", type = str, help = "PAM sequence", default = "NGG")
    parser.add_argument("--window_st", type = int, help = "start position of base editing window (default = 13)", default = 13)
    parser.add_argument("--window_ed", type = int, help = "end position of base editing window (default = 17)", default = 17)
    parser.add_argument("-r", "--ref_base", type = str, help = "reference base (default = C)", default = 'C')
    parser.add_argument("-m", "--modified_base", type = str, help = "changed base (default = T)", default = 'T')
    parser.add_argument("--mismatch", type = int, help = "Mismatche", default = 2)
    parser.add_argument("ref_dir", type = str, help = "directory of reference genome")
    parser.add_argument("--output", type = str, help = "output file name", default = "result.txt")

    return parser.parse_args()

def main():

    args = parse_args()

    b = BEDesigner(parse_args())

    b.find_targets()
    b.run_cas_offinder()

if __name__ == '__main__':
    main()
#python3 be_designer.py GGTCTACACAAGCTTCCACCAGGAAGGAAATATGGTGCAGTAAGAGAGTGCTTTTTCTTACACCTCAGGTCATGGTAAATGACCTTTCTAGAGGAGCTTGTCCCGCTGCTGAAATAAAGTGTTTAGTTATTGATGAAGCTCATAAAGCTCTCGGAAACTATGCTTATTGCCAG NGG D:/Genome/GRCh38
