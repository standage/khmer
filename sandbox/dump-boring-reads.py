#! /usr/bin/env python
# This file is part of khmer, https://github.com/dib-lab/khmer/, and is
# Copyright (C) 2016, The Regents of the University of California.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of the Michigan State University nor the names
#       of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written
#       permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Contact: khmer-project@idyll.org
"""
% python sandbox/id-variants.py [ -C <cutoff> ] <data1> <data2> ...

Use '-h' for parameter help.

TODO: add to sandbox README
"""
from __future__ import print_function

import argparse
import sys

import screed
import khmer
from khmer.khmer_args import build_nodegraph_args, create_nodegraph


def main():
    parser = build_nodegraph_args()
    parser.add_argument('-m', '--allowed-mismatches', type=int, default=0)
    parser.add_argument('reference')
    parser.add_argument('readfiles', nargs='+')
    args = parser.parse_args()

    if not args.quiet:
        print('\nPARAMETERS:', file=sys.stderr)
        print(' - kmer size =    %d \t\t(-k)' % args.ksize, file=sys.stderr)
        print(' - n hashes =     %d \t\t(-N)' % args.n_tables, file=sys.stderr)
        print(' - min hashsize = %-5.2g \t(-x)' % \
            args.max_tablesize, file=sys.stderr)
        print('', file=sys.stderr)
        print('Estimated memory usage is %.2g bytes ' \
            '(n_hashes x min_hashsize)' % \
            (args.n_tables * args.max_tablesize), file=sys.stderr)
        print('-' * 8, file=sys.stderr)

    refrgraph = create_nodegraph(args)
    nreads, nconsumed = refrgraph.consume_fasta(args.reference)
    if not args.quiet:
        message = 'consumed {:d} reads and {:d} bp from {:s}'.format(
            nreads,
            nconsumed,
            args.reference
        )
        print(message, file=sys.stderr)

    readskept, readsdiscarded = 0, 0
    for readfile in args.readfiles:
        for n, record in enumerate(screed.open(readfile)):
            kmers_present = 0
            for kmer in refrgraph.get_kmers(record.sequence):
                if refrgraph.get(kmer):
                    kmers_present += 1
                    if kmers_present > args.allowed_mismatches:
                        break
            if kmers_present > args.allowed_mismatches:
                readskept += 1
                print('@', record.name, '\n', record.sequence, '\n+\n',
                      record.quality, sep='')
            else:
                readsdiscarded += 1

    if not args.quiet:
        total = readskept + readsdiscarded
        keptperc = float(readskept) / float(total) * 100
        discperc = float(readsdiscarded) / float(total) * 100
        message = (
            'kept {} reads ({:.1f}%%), discarded {} reads ({:.1f}%%)'.format(
                readskept, keptperc,
                readsdiscarded, discperc
            )
        )
        print(message, file=sys.stderr)


if __name__ == '__main__':
    main()
