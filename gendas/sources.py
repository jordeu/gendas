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
Module with all the available data sources
"""

import csv
import gzip
import logging
import tabix
from collections import OrderedDict
from collections import defaultdict

from intervaltree import IntervalTree
from gendas.utils import _skip_partitions, _skip_comments

logger = logging.getLogger("gendas")


class GendasSource:
    """
    Abstract class to define the source interface
    """
    def __init__(self, sequence=None, begin=None, end=None, header=None, ctypes=None):
        """
        Initialize a source

        Args:
            sequence: Sequence identifier
            begin: lower genomic position
            end: higher genomic position
            header: List of column identifiers
            ctypes: List of column data types
        """
        self.label = None
        self.sequence = sequence
        self.begin = begin
        self.end = end
        self.header = header
        self.ctypes = ctypes

    def index(self, label: str):
        """
        An iterable over all the possible values of an indexed column and a list of genomic regions that contain
        that value.

        Example of a coding genes dataset indexed by gene:
        ('TP53', [('17', 7565260, 7565332), ... , ('17', 7569527, 7569562)]),
        ('BRCA1', [('17', 41197698, 41197819), ..., ('17', 41277199, 41277202)]) ...

        Args:
            label: The label of an indexed data column
        """
        raise NotImplementedError()

    def query(self, sequence, begin, end):
        """
        Returns a generator that iterates all the rows in a 'sequence' from 'begin' to
        'end' (both included). Each row it's a dictionary following this format:
        { 'field_01_key': field_01_value, ... }

        Args:
            sequence: Sequence identifier
            begin: Start position in the sequence
            end: End position in the sequence
        """
        raise NotImplementedError()

    def intersect(self, sequence, begin, end):
        """
        Get all the available intervals that contains some data in a given region

        Args:
            sequence:
            begin:
            end:

        Returns:
            A generator of tuples like (sequence_id, begin, end)
        """
        raise NotImplementedError()

    def __iter__(self, p=None):
        """
        Iterate the whole data source

        Args:
            p:
        """
        raise NotImplementedError()

    def __len__(self):
        """
        Returns: Total number of rows in the dataset
        """
        raise NotImplementedError()


class TabixSource(GendasSource):
    """
    Source to query tabix formatted genomic files.

    Tabix files are normal sorted tabulated text files compressed using block-gzip and indexed by block.
    In a way that, given a region, you only need to decompress the blocks that contain that region to access
    them.
    """
    def __init__(self, filename, sequence=None, begin=None, end=None, header=None, ctypes=None, indices=None):
        """
        Initialize a tabix source

        Args:
            filename: Path to a tabix file
            sequence: Header label that identifies the sequence column
            begin: Header label that identifies the begin position column
            end: Header label that identifies the end position column
            header: Ordered list with all the column headers
            ctypes: Ordered list with all the column data types
            indices: List with all the columns that we want to create an index
        """
        super().__init__(sequence=sequence, begin=begin, end=end, header=header, ctypes=ctypes)

        self.sequence_idx = self._idx(sequence)
        self.begin_idx = self._idx(begin)
        self.end_idx = self._idx(end)
        self.tb = None
        self.filename = filename

        self.indices = {} if indices is None else {self._idx(i): OrderedDict() for i in indices}
        if len(self.indices) > 0:
            with gzip.open(filename, 'rt') as fd:
                reader = csv.reader(fd, delimiter='\t')
                for row in reader:
                    for i in self.indices.keys():
                        if row[i] not in self.indices[i]:
                            self.indices[i][row[i]] = []
                        self.indices[i][row[i]].append(
                            (row[self.sequence_idx], int(row[self.begin_idx]), int(row[self.end_idx])))

    def index(self, label: str):
        return self.indices[self._idx(label)].items()

    def _tabix(self):
        try:
            if self.tb is None:
                self.tb = tabix.open(self.filename)
        except tabix.TabixError:
            msg = "Error opening tabix file {}".format(self.filename)
            logger.error(msg)
            raise RuntimeError(msg)

        return self.tb

    def query(self, sequence, begin, end):
        try:
            for row in self._tabix().query(sequence, begin, end):
                yield {h: c(v) for c, v, h in zip(self.ctypes, row, self.header)}
        except tabix.TabixError:
            logger.error("Fail tabix query {}:{}-{} at {}".format(sequence, begin, end, self.filename))

    def intersect(self, sequence, begin, end):
        for row in self._tabix().query(sequence, begin, end):
            yield row[self.sequence_idx], int(row[self.begin_idx]), int(row[self.end_idx]) + 1

    def _idx(self, label):
        if type(label) == int:
            return label

        return self.header.index(label)

    def __iter__(self, p=None):
        with gzip.open(self.filename, 'rt') as fd:

            if p is None:
                it = fd
            else:
                it = _skip_partitions(fd, p)

            # Skip comments
            it = _skip_comments(it, '#')

            reader = csv.reader(it, delimiter='\t')
            for r in reader:
                yield {h: c(v) for c, v, h in zip(self.ctypes, r, self.header)}

    def __getstate__(self):
        state = dict(self.__dict__)
        state['tb'] = None
        state['indices'] = {}
        return state


class IntervalTreeSource(GendasSource):
    """
    Source that loads a typical genomic regions data file all in memory as an interval tree.

    The current implementation only supports simple tabulated text files, but more standard formats like
    VCF and MAF files will be added later.
    """
    def __init__(self, filename, sequence=None, begin=None, end=None, header=None, ctypes=None, indices=None):
        """
        Initialize a regions source dataset

        Args:
            filename: Path to a gzip compressed tabulated text file
            sequence: Header that identifies the sequence column
            begin: Header that identifies the begin column
            end: Header that identifies the end column
            header: Ordered list with all the column headers
            ctypes: Ordered list with all the column data types
            indices: List with all the columns that we want to create an index
        """
        super().__init__(sequence=sequence, begin=begin, end=end, header=header, ctypes=ctypes)

        self.sequence_idx = self._idx(sequence)
        self.begin_idx = self._idx(begin)
        self.end_idx = self._idx(end)
        self.tb = None
        self.filename = filename

        self.indices = {} if indices is None else {self._idx(i): OrderedDict() for i in indices}
        if len(self.indices) > 0:
            with gzip.open(filename, 'rt') as fd:
                reader = csv.reader(fd, delimiter='\t')
                for row in reader:
                    for i in self.indices.keys():
                        if row[i] not in self.indices[i]:
                            self.indices[i][row[i]] = []
                        self.indices[i][row[i]].append(
                            (row[self.sequence_idx], int(row[self.begin_idx]), int(row[self.end_idx])))

        self._trees = defaultdict(IntervalTree)
        with gzip.open(filename, 'rt') as fd:
            reader = csv.reader(fd, delimiter='\t')
            for r in reader:
                self._trees[r[self.sequence_idx]][int(r[self.begin_idx]):int(r[self.end_idx]) + 1] = \
                    {h: c(v) for c, v, h in zip(self.ctypes, r, self.header)}

    def index(self, label: str):
        return self.indices[self._idx(label)].items()

    def query(self, sequence, begin, end):
        for row in self._trees[sequence][begin:end]:
            yield row.data

    def intersect(self, sequence, begin, end):
        for row in self._trees[sequence][begin:end]:
            yield sequence, row.begin, row.end + 1

    def _idx(self, label):
        if type(label) == int:
            return label

        return self.header.index(label)

    def __iter__(self):
        with gzip.open(self.filename, 'rt') as fd:
            reader = csv.reader(fd, delimiter='\t')
            for r in reader:
                yield {h: c(v) for c, v, h in zip(self.ctypes, r, self.header)}


class PandasSource(GendasSource):
    """
    A data source that queries a Pandas dataframe loaded in memory
    """
    def __init__(self, df, sequence=None, begin=None, end=None):
        """
        Args:
            df: A pandas dataframe
            sequence: The sequence column label (Ex: CHROMOSOME)
            begin: The begin position column label (Ex: START)
            end: The end position column label (Ex: STOP)
        """
        super().__init__(sequence=sequence, begin=begin, end=end, header=df.columns.tolist(), ctypes=df.dtypes.tolist())
        self.df = df

    def index(self, label: str):
        return label

    def query(self, sequence, begin, end):
        raise NotImplementedError()

    def intersect(self, sequence, begin, end):
        raise NotImplementedError()

    def __iter__(self, p=None):
        for r in self.df.iterrows():
            yield r[1].to_dict()

    def __len__(self):
        return len(self.df)
