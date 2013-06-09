#!/usr/bin/python
"""
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import string
import os
import optparse
import numpy
import math
import random
import re
import time
from subprocess import call, Popen, PIPE
from itertools import product
from libkmersvm import *
import difflib

def get_overlaps(x,y):
    x_overlaps = []
    y_overlaps = []
    #check left of x right of y
    for i in xrange(1, len(x)):
        if x[:i] == y[-i:]:
            x_overlaps.append(x[:i])
    #check left of x right of y
    for i in xrange(1, len(y)):
        if y[:i] == x[-i:]:
            y_overlaps.append(y[:i])

    return x_overlaps, y_overlaps

def get_longest_overlap(x,y,k):
    max_overlap = ''
    x_overlaps, y_overlaps = get_overlaps(x,y)
    for i in x_overlaps:
        if len(i) > len(max_overlap) and len(i) < k:
            max_overlap = i
    for i in y_overlaps:
        if len(i) > len(max_overlap) and len(i) < k:
            max_overlap = i

    return max_overlap
    

def get_pair_frequencies(seqs, k, w):

    pair_freq = {}
    total = 0.0
    for seq in seqs:
        total += (len(seq)-w+1)*(w-k)

    test = 0.0

    progress_count = 0
    t = len(seqs)
    for seq in seqs:
        if progress_count % math.ceil(t/10000.0) == 0:
            p = (float(progress_count)/t)*100
            sys.stderr.write("\r%.2f%% progress. " % p)
            sys.stderr.flush()
        progress_count += 1

        for i in xrange(0, len(seq)-w+1):
            a = seq[i:i+k]
            window = seq[i+1:i+w]
            for j in xrange(0, len(window)-k+1):
                b = window[j:j+k]
                test += 1.0
                if (a,b) in pair_freq:
                    pair_freq[(a,b)] += 1.0/total
                elif (b,a) in pair_freq:
                    pair_freq[(b,a)] += 1.0/total
                else:
                    pair_freq[(a,b)] = 1.0/total

    sys.stderr.write("\n")

    return pair_freq,total


def get_pair_t_statistic(pair_freqs,freqs,k,w,total):

    pair_t_statistic = {}

    for (p,v) in pair_freqs.items():
        (a,b) = p
        f_a = freqs[len(a)-1][a]
        f_b = freqs[len(b)-1][b]

        #calculate mu of null
        mu = 0.0

        #approximation: consider biggest overlap
        overlap = get_longest_overlap(a,b,k)
        l = len(overlap)
        if l > 0:
            if a[:l] == overlap:
                mu += f_b * freqs[len(a[l:])-1][a[l:]]
            elif b[:l] == overlap:
                mu += f_a * freqs[len(b[l:])-1][b[l:]]

        mu += f_a*f_b

        denom = math.sqrt( v/total )
        t_statistic = (v - mu)/denom

        pair_t_statistic[(a,b)] = t_statistic

    return pair_t_statistic


def save_results(pair_t_statistic, pair_freq, freq, output):
    ptstat = []
    for (p,v) in pair_t_statistic.items():
        (a,b) = p
        f_a = freq[len(a)-1][a]
        f_b = freq[len(b)-1][b]
        f_ab = 0.0
        if (a,b) in pair_freq:
            f_ab = pair_freq[(a,b)]
        else:
            f_ab = pair_freq[(b,a)]

        ptstat.append( (a,b,f_a,f_b,f_ab,v) )

    out = sorted(ptstat, key=lambda p: p[5], reverse=True)
    f = open(output, 'w')
    f.write('\t'.join(["seq1", "seq2","f1","f2","f_12","t-statistic"]) + "\n")
    for (a,b,f_a,f_b,f_ab,v) in out: 
        f.write('\t'.join([a, b, str(f_a), str(f_b), str(f_ab), str(v)]) + "\n")
    f.close()


def get_frequencies(seqs, k):
    freqs = {}

    total = 0.0
    for seq in seqs:
        total += len(seq)-k+1
    
    for seq in seqs:
        for i in xrange(0, len(seq)-k+1):
            kmer = seq[i:i+k]
            if kmer not in freqs:
                freqs[kmer] = 0.0
            freqs[kmer] += 1.0/total

    return freqs
       

def main(argv = sys.argv):
    usage = "Usage: %prog [options] SEQ_FILE OUTPUT_NAME"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-k", dest="kmerlen", type=int, \
            help="set kmer length", default=6)

    parser.add_option("-w", dest="window", type=int, \
            help="set window length", default=50)

    parser.add_option("-q", dest="quiet", default=False, action="store_true", \
            help="supress messages (default=false)")

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    if len(args) != 2:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit(0)

    seqf = args[0]
    outm = args[1]

    seqs, sids = read_fastafile(seqf)

    freqs = [{}] * (options.kmerlen)
    for i in xrange(0,options.kmerlen):
        freqs[i] = get_frequencies(seqs, i+1)

    pair_freq, total = get_pair_frequencies(seqs,options.kmerlen,options.window)

    pair_t_statistic = get_pair_t_statistic(pair_freq,freqs,options.kmerlen,options.window,total)
    save_results(pair_t_statistic,pair_freq,freqs,outm)

if __name__ =='__main__': main()


