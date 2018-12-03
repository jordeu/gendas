#
#   Copyright 2018 Jordi Deu-Pons
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not use this
#   file except in compliance with the License. You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software distributed under
#   the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF
#   ANY KIND, either express or implied. See the License for the specific language
#   governing permissions and limitations under the License.
#

import gzip
import struct
from os import path
from pprint import pprint

import click


class TabixIndex:
    """
    A pure python implementation to query a tabix file
    """
    MAX_BIN = 37450
    TAD_MIN_CHUNK_GAP = 32768
    TAD_LIDX_SHIFT = 14

    def __init__(self, indexfile):

        with gzip.open(indexfile, 'rb') as idx:

            # Magic number
            magic = idx.read(4)
            if magic != b'TBI\x01':
                raise RuntimeError("[tabix index] wrong magic number")

            # idx->n
            self.mSeq = struct.unpack('i', idx.read(4))[0]

            # typedef struct
            # {
            #     int32_t preset;
            #     int32_t sc, bc, ec; // seq col., beg col. and end col.
            #     int32_t meta_char, line_skip;
            # } ti_conf_t;
            self.conf = {
                'preset': struct.unpack('i', idx.read(4))[0],
                'sc': struct.unpack('i', idx.read(4))[0],
                'bc': struct.unpack('i', idx.read(4))[0],
                'ec': struct.unpack('i', idx.read(4))[0],
                'meta_char': struct.unpack('i', idx.read(4))[0],
                'line_skip': struct.unpack('i', idx.read(4))[0]
            }

            # Read target names
            l = struct.unpack('i', idx.read(4))[0]
            buf = idx.read(l)

            self.names = []
            buf_str = ""
            for i in range(l):
                if buf[i] == 0:
                    self.names.append(buf_str)
                    buf_str = ""
                else:
                    buf_str += chr(buf[i])

            self.binning = {}
            self.linear = {}
            for i in range(self.mSeq):

                # Load binning index
                n_bin = struct.unpack('I', idx.read(4))[0]

                bins = {}
                for j in range(n_bin):
                    bin = struct.unpack('I', idx.read(4))[0]
                    bin_l = struct.unpack('I', idx.read(4))[0]

                    chunks = []
                    for k in range(bin_l):
                        u = struct.unpack('Q', idx.read(8))[0]
                        v = struct.unpack('Q', idx.read(8))[0]
                        chunks.append((u, v))
                    bins[bin] = chunks

                self.binning[self.names[i]] = bins

                # load linear index
                l = struct.unpack('I', idx.read(4))[0]
                offsets = []
                for k in range(l):
                    offsets.append(struct.unpack('Q', idx.read(8))[0])

                self.linear[self.names[i]] = {'size': l, 'offset': offsets}

            # Check other bytes
            b = idx.read(1)
            c = 1 if b else 0
            while b:
                b = idx.read(1)
                c += 1
            self.more = c

    @staticmethod
    def _reg2bins(beg, _end):

        bins = []
        end = _end

        if beg >= end:
            return bins

        if end >= 1 << 29:
            end = 1 << 29

        end -= 1
        bins.append(0)

        k = 1 + (beg >> 26)
        while k <= 1 + (end >> 26):
            bins.append(k)
            k += 1

        k = 9 + (beg >> 23)
        while k <= 9 + (end >> 23):
            bins.append(k)
            k += 1

        k = 73 + (beg >> 20)
        while k <= 73 + (end >> 20):
            bins.append(k)
            k += 1

        k = 585 + (beg >> 17)
        while k <= 585 + (end >> 17):
            bins.append(k)
            k += 1

        k = 4681 + (beg >> 14)
        while k <= 4681 + (end >> 14):
            bins.append(k)
            k += 1

        return bins

    def query(self, seq, begin: int, end: int):
        """
        TODO Implement the tabix query
        """
        if seq in self.names:

            # Bins
            idx_b = self.binning[seq]
            idx_l = self.linear[seq]

            bins = self._reg2bins(begin, end)

            # Linear offset
            l_length = idx_l['size']
            l_offsets = idx_l['offset']
            if l_length > 0:
                if begin >> TabixIndex.TAD_LIDX_SHIFT >= l_length:
                    min_off = l_offsets[l_length - 1]
                else:
                    min_off = l_offsets[begin >> TabixIndex.TAD_LIDX_SHIFT]
            else:
                min_off = 0

            n_off = 0
            for bin in bins:
                if bin in idx_b:
                    n_off += len(idx_b[bin])

            if n_off == 0:
                return

            return n_off


@click.command()
@click.argument('filename', nargs=1, type=click.Path(exists=True))
@click.option('-s', '--sequence', type=str, help='Genomic sequence')
@click.option('-b', '--begin', type=int, help='Begin position')
@click.option('-e', '--end', type=int, help='End position')
def cmdline(filename, sequence, begin, end):
    idx = TabixIndex("{}.tbi".format(path.expanduser(filename)))
    pprint(idx.query(sequence, begin, end))


if __name__ == "__main__":
    cmdline()
