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

"""
    Experimental sources implementations
"""

import os
import mmap
import bgdata

from gendas.sources import GendasSource


class HG19Source(GendasSource):
    """
    A data source that represents the Human Genome release hg19
    """
    def __init__(self):
        super().__init__(sequence="CHR", begin="BEGIN", end="END", header=["CHR", "BEGIN", "END", "SEQ"],
                         ctypes=[str, int, int, str])

        self.hg19 = bgdata.get_path('datasets', 'genomereference', 'hg19')
        self.mmap = {}

    def intersect(self, sequence, begin, end):
        yield sequence, begin, end

    def query(self, sequence, begin, end):
        yield HG19Sequence(self, sequence, begin, end)

    def __len__(self):
        pass

    def _seq_mmap(self, seq):
        if seq not in self.mmap:
            chr_file = "chr{0}.txt".format(seq)
            fd = open(os.path.join(self.hg19, chr_file), 'rb')
            self.mmap[seq] = mmap.mmap(fd.fileno(), 0, access=mmap.ACCESS_READ)
        return self.mmap[seq]

    def get_ref(self, seq, start, size=1):
        mm_file = self._seq_mmap(seq)
        mm_file.seek(start - 1)
        return mm_file.read(size).decode().upper()


class HG19Sequence:

    def __init__(self, source: HG19Source, sequence: str, begin: int, end: int):
        self.source = source
        self.sequence = sequence
        self.begin = begin
        self.end = end

    def __getitem__(self, item):
        if item == "CHR":
            return self.sequence
        if item == "BEGIN":
            return self.begin
        if item == "END":
            return self.end
        if item == "SEQ":
            return self.source.get_ref(self.sequence, self.begin + 1, self.end - self.begin)

        if type(item) == slice:
            start = self.begin + item.start
            end = self.end + item.stop
            return self.source.get_ref(self.sequence, start + 1, end - start)

        raise KeyError()
