#!/usr/bin/env python3

class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

class SequenceIO:
    # Stolen from https://github.com/lh3/readfq/blob/master/readfq.py
    def each(self, fp): # this is a generator function
        
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs)  # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq  # yield a fasta record instead
                    break

    def each_sequence(self, fp):
        '''Like each except iterate over Sequence objects'''
        for name, seq, _ in self.each(fp):
            yield Sequence(name, seq)

    def read_fasta_file(self, path_to_fasta_file):
        seqs = []
        for name, seq, _ in self.each(open(path_to_fasta_file)):
            seqs.append(Sequence(name, seq))
        return seqs
    
    def write_fasta_file(self, sequence_objects, path_to_fasta_file):
        with open(path_to_fasta_file,'w') as f:
            self.write_fasta(sequence_objects, f)

    def write_fasta(self, sequence_objects, io):
        for s in sequence_objects:
            io.write(">")
            io.write(s.name)
            io.write("\n")
            io.write(s.seq)
            io.write("\n")
