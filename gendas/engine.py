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

import logging
import os
from os.path import join, dirname

from configobj import ConfigObj, Section
from pathos.pools import ProcessPool, ParallelPool
from tqdm import tqdm

from gendas.sources import GendasSource, TabixSource, IntervalTreeSource
from gendas.statistics import count
from gendas.utils import _get_chunks, _overlap_intervals

logger = logging.getLogger("gendas")

SOURCE_TYPES = {
    'tabix': TabixSource,
    'mem': IntervalTreeSource
}


class Gendas:
    """
        Gendas main engine that represents all the loaded datasets.
        All the queries start here.
    """

    def __init__(self, configfile: 'str' = None, workers: 'int' = os.cpu_count(), servers=None, progress: 'int' = 20):
        """
        Initialize a gendas engine

        Args:
            configfile: A file that contains the datasets definitions. Defaults to None, an empty engine.
            workers: Total number of workers to parallelize the computations. Defaults to total number of cores.
            servers: A list of servers where to distribute the parallelization. Defaults to only localhost.
            progress: A smaller number means that gendas will report progress more often. Defaults to 20.

        """
        self.workers = workers
        self.servers = servers
        self.progress = progress
        self.sources = {}

        if configfile is not None:
            if not os.path.exists(configfile):
                raise FileNotFoundError("File {} not found".format(configfile))

            # Load datasets from config file
            config = ConfigObj(configfile)
            for key, section in config.items():

                # Skip general parameters
                if type(section) != Section:
                    continue

                # Load the source type
                source = SOURCE_TYPES[str(section['type']).lstrip().lower()]

                # Create a source instance from the configuration
                self[key] = source(
                    join(dirname(configfile), section['file']),
                    header=section.get('header', None),
                    ctypes=None if section.get('ctypes', None) is None else [eval(t) for t in section['ctypes']],
                    sequence=section['sequence'],
                    begin=section['begin'],
                    end=section['end'],
                    indices=section.get('indices', None)
                )

    def __setitem__(self, label: 'str', source: object) -> object:
        """
        Add a source

        Args:
            label: Source name
            source: Source object
        """
        source.label = label
        self.sources[label] = source

    def __getitem__(self, source: 'str') -> 'GendasDataset':
        """
        Returns a dataset view of the given source

        Args:
            source: The source key

        Returns: A dataset view of the source

        """
        return GendasDataset(self.sources[source], self)

    def groupby(self, field) -> 'GendasGroupBy':
        """
        Create a groupby view by the given field

        Args:
            field: The field to group

        Returns: A groupby view of the current dataset

        """
        return GendasGroupBy(field, self)

    def pool(self):
        """
        Returns: The computing pool to process run the queries

        """
        if self.servers is None:
            return ProcessPool(nodes=self.workers)
        else:
            return ParallelPool(nodes=self.workers, servers=self.servers)


class GendasSlice:
    """
    A gendas slice is a view of only some genomic regions (segments) of the whole genome.
    """

    def __init__(self, manager: Gendas, segments: list):
        """

        Args:
            manager: A gendas engine
            segments: The genomic segments of interest. A list of tuples like (chromosome, start, end)
        """
        self.manager = manager
        self.segments = segments

    def __getitem__(self, source):
        """

        Args:
            source: The source of interest

        Returns: A slice view limiting to the given source

        """
        return GendasSliceDataset(self.manager.sources[source], self)


class GendasDataset:
    """
        A gendas dataset is the interface (a dataframe like structure) to access
        and manipulate the data from a gendas source
    """

    def __init__(self, source: GendasSource, manager: Gendas):
        self.manager = manager
        self.source = source

    def __getitem__(self, field):
        """

        Args:
            field: The column of interest

        Returns: A view of only one column of the dataset

        """
        return GendasColumn(field, self)

    def merge(self, right, on=None):
        """
        Join this dataset with another dataset. By default (without any details) the datasets are always merge using
        the genomic coordinates, but you can add more restrictions using the 'on' arguments.

        Args:
            right: The other dataset to join with.
            on: A list of extra fields to add to the join (both datasets need to have the same field names)

        Returns: A gendas merge view of both datasets

        """
        return GendasMerge(self, right, on=on)

    def map(self, fn):
        """
         Apply a function to all the rows in this dataset.

        Args:
            fn: The function to apply to each row.

        Returns: A generator to the results.

        """
        return self._map_par(fn)

    def _map_seq(self, fn):
        """
        Sequential implementation of the map (for debugging purposes only)
        """
        return map(fn, self)

    def _map_par(self, fn):
        """
        Parallel implementation of the map
        """
        cores = self.manager.workers
        logger.debug("ready to create pool")
        with self.manager.pool() as executor:
            logger.debug("pool created")
            partitions = cores * self.manager.progress
            mapfn = lambda p: list(map(fn, self.__iter__(p=(p, partitions))))
            for items in executor.imap(mapfn, range(cores)):
                for item in items:
                    yield item

    def filter(self, fn):
        """
        Filter the rows of this dataset

        Args:
            fn: A function that return true/false to filter the rows

        Returns:
            A filtered view of this dataset

        """
        return GendasDatasetFilter(self, fn)

    def count(self, progress=False):
        """
        Count how many rows has this dataset

        Args:
            progress: True to show progress

        Returns:
            Total number of rows

        """
        return self._count_seq(progress=progress)

    def _count_seq(self, progress=False):
        """
        count sequential implementation (for debugging purposes)
        """
        logger.debug("Count sequencial")
        if progress:
            return count(tqdm(self))
        return count(self)

    def _count_par(self, progress=False):
        """
        count parallel implementation
        """
        logger.debug("Count parallell")
        cores = self.manager.workers
        logger.debug("ready to create pool")
        with self.manager.pool() as executor:
            logger.debug("pool created")
            partitions = cores * self.manager.progress
            mapfn = lambda p: list(map(count, self.__iter__(p=(p, partitions))))
            it = executor.uimap(mapfn, range(cores))
            if progress:
                it = tqdm(it, total=cores)
            return sum(it)

    def head(self, n=10):
        """
        Returns a generator to iterate the first 'n' rows.

        Args:
            n: Total number of rows to return
        """
        for i, r in enumerate(self, start=1):
            if i > n:
                break
            yield r
            if i == n:
                break

    def _rows(self, p=None):
        """
        Private implementation to iterate over this dataset

        Args:
            p: Partition parameter. Useful when iterating in parallel.

        """
        return self.source.__iter__(p=p)

    def __iter__(self, p=None):
        """
        Iterate this dataset rows
        """
        for r in self._rows(p=p):
            yield r

    def __len__(self):
        """
        Returns: How many rows has this dataset
        """
        return self.count()


class GendasDatasetFilter(GendasDataset):
    """
        A filtered view of a Gendas dataset
    """

    def __init__(self, dataset: 'GendasDataset', filter):
        """
        Args:
            dataset: Dataset to filter
            filter: Filtering function
        """
        super().__init__(dataset.source, dataset.manager)
        self.dataset = dataset
        self.filter = filter

    def __iter__(self, p=None):
        return filter(self.filter, self.dataset)


class GendasMergeDataset(GendasDataset):
    """
        A dataset view of only one dataset in a gendas merges
    """

    def __init__(self, merge: 'GendasMerge', source: 'GendasSource'):
        super().__init__(source, merge.left.manager)
        self.merge = merge

    def _rows(self, p=None):
        for r in self.merge.__iter__(p=None):
            yield r[self.source.label]

    def __len__(self):
        return len(self.merge)


class GendasMerge(GendasDataset):
    """
        A view that joins two datasets
    """

    def __init__(self, left: 'GendasDataset', right: 'GendasDataset', on: list):
        """

        Args:
            left: The left dataset
            right: The right dataset
            on: Extra columns to do the join (both datasets need to have the same column label)
        """
        super().__init__(None, left.manager)

        self.left = left
        self.right = right
        self.on = on
        self.sources = {
            left.source.label: left.source,
            right.source.label: right.source
        }

    def __getitem__(self, source):
        """
        Get dataset view of only one of the merged datasets.

        Args:
            source: The source label

        Returns: A gendas dataset view

        """
        return GendasMergeDataset(self, self.sources[source])

    def merge(self, right, on=None):
        """
        Merge more datasets with this merge view.

        Args:
            right: The right dataset to merge
            on:  A list with extra columns to do the join

        Returns: A gendas merge view

        """
        return GendasMultipleMerge(self, right, on=on)

    def filter(self, fn):
        """
        Filter a merge view

        Args:
            fn: The filtering function

        Returns: A gendas merge view filtered

        """
        return GendasMergeFilter(self, fn)

    def _rows(self, p=None):
        """
        Iterate the left most dataset computing the 'inner join' merge
        of other datasets

        Args:
            p: partition. Internal parameter to use when doing iterations in parallel
        """
        for l_row in self.left.__iter__(p=p):
            seq = l_row[self.left.source.sequence]
            begin = l_row[self.left.source.begin]
            end = l_row[self.left.source.end]

            if self.on is not None:
                l_key = [l_row[o] for o in self.on]
            else:
                l_key = None

            r_rows = []
            for r in self.right.source.query(seq, begin - 1, end):
                if self.on is not None:
                    r_key = [r[o] for o in self.on]
                    if l_key != r_key:
                        continue

                r_rows.append(r)

            # Inner join
            for r_row in r_rows:
                yield {
                    self.left.source.label: l_row,
                    self.right.source.label: r_row
                }

    def __iter__(self, p=None):
        """
            Return a generator to iterate this
        """
        for r in self._rows(p=p):
            yield r

    def __len__(self):
        # TODO Implement the gendas merge count (without iterating everything)
        raise NotImplementedError("GendasMerge.__len__ not implemented")


class GendasMergeFilter(GendasMerge):
    """
    Filtered view of a gendas merge
    """

    def __init__(self, merge: 'GendasMerge', filter):
        super().__init__(merge.left, merge.right, merge.on)
        self.merge = merge
        self.filter = filter

    def __iter__(self, p=None):
        return filter(self.filter, self.merge)


class GendasMultipleMerge(GendasMerge):
    """
    Merge more than two datasets.
    """

    def __init__(self, merge: 'GendasMerge', right: 'GendasDataset', on: list):
        super().__init__(merge.left, right, on)
        self.merge = merge

        # Add other sources
        for k, v in merge.sources.items():
            self.sources[k] = v

    def _rows(self, p=None):
        for m_row in self.merge.__iter__(p=p):

            l_row = m_row[self.left.source.label]
            seq = l_row[self.left.source.sequence]

            begin, end = _overlap_intervals(
                [(m_row[s.label][s.begin], m_row[s.label][s.end])
                 for s in self.merge.sources.values()]
            )

            if self.on is not None:
                m_key = []
                for o in self.on:
                    for label, source in self.merge.sources.items():
                        if o in source.header:
                            m_key.append(m_row[label][o])
            else:
                m_key = None

            r_rows = []
            for r in self.right.source.query(seq, begin - 1, end):
                if self.on is not None:
                    r_key = [r[o] for o in self.on]
                    if m_key != r_key:
                        continue

                r_rows.append(r)

            # Inner join
            for r_row in r_rows:
                res = {k: v for k, v in m_row.items()}
                res[self.right.source.label] = r_row
                yield res


class GendasSliceDataset(GendasDataset):
    """
    A dataset view of a source filtered by a gendas slice (a genomic regions definition)
    """

    def __init__(self, source: 'GendasSource', slice: 'GendasSlice'):
        super().__init__(source, slice.manager)
        self.slice = slice
        self.rows = None

    def _rows(self, p=None):
        if self.rows is None:
            self.rows = []
            for s in self.slice.segments:
                self.rows += list(self.source.query(s[0], s[1], s[2]))

        return self.rows


class GendasColumn:
    """
    A view of only one column of a dataset
    """

    def __init__(self, label: str, dataset: GendasDataset):
        self.dataset = dataset
        self.label = label

    def __iter__(self):
        for r in self.dataset:
            yield r[self.label]

    def __len__(self):
        return len(self.dataset)


class GendasGroupBy:
    """
    A grouped view of a dataset. That lets you apply an agreggating function to each group.
    """

    def __init__(self, field: 'GendasColumn', manager: 'Gendas'):
        """
        Args:
            field: The column that defines a group
            manager: A gendas manager
        """
        self.field = field
        self.manager = manager

    def aggregate(self, aggregator, **kwargs):
        """
        Returns a generator that returns the result of apply the aggregator function
        to each group.

        Args:
            aggregator: An aggregation function
            **kwargs: Extra parameters to pass to the aggregation function

        Returns: A generator

        """
        return self._aggregate_par(aggregator, **kwargs)

    def _compute(self, aggregator, args, groups) -> dict:
        label, segments = groups
        v = {self.field.label: label}
        partition = GendasSlice(self.manager, segments)
        if type(aggregator) == dict:
            for f, aggregator in aggregator.items():
                if len(args) > 0:
                    v[f] = aggregator(partition, **args)
                else:
                    v[f] = aggregator(partition)
        else:
            if len(args) > 0:
                v = aggregator(partition, v, **args)
            else:
                v = aggregator(partition, v)
        return v

    def _compute_par(self, aggregator, args, groups):
        result = [self._compute(aggregator, args, group) for group in groups]
        return result

    def _mapfn(self, r):
        return self._compute_par(self.aggregator, self.kwargs, r)

    def _aggregate_par(self, aggregator: dict, **kwargs):
        """
        Parallel implementation of the aggregate method
        """

        cores = self.manager.workers
        regions = self.field.dataset.source.index(self.field.label)
        logger.debug("Retrive valid column names")
        labels = set(self.field)
        logger.debug("Retrive regions to aggregate")
        regions = list(filter(lambda r: r[0] in labels, regions))
        regions_size = len(regions)
        partitions = cores * self.manager.progress
        chunksize = (regions_size // partitions) + 1
        regions = list(_get_chunks(regions, size=chunksize))
        logger.debug(
            "{} chunks of {} regions (total {}) to run in {} partitions at {} cores".format(len(regions), chunksize,
                                                                                            regions_size, partitions,
                                                                                            cores))
        # mapfn = lambda r: self._compute_par(aggregator, kwargs, r)

        self.aggregator = aggregator
        self.kwargs = kwargs

        with self.manager.pool() as executor:
            logger.debug("pool created")
            for items in executor.uimap(self._mapfn, regions):
                for item in items:
                    yield item

    def _aggregate_seq(self, fields, **kwargs):
        """
        Sequential implementation of the aggregate method (for testing/debugging purposes)
        """
        regions = self.field.dataset.source.index(self.field.label)
        labels = set(self.field)
        regions = list(filter(lambda r: r[0] in labels, regions))
        for label, segments in regions:
            yield self._compute(fields, kwargs, (label, segments))
