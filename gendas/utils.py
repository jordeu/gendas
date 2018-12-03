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
    Help functions used in several modules.
"""


def flatten(iterable):
    """
    A generator that iterates a collection of nested lists as a single collection

    Args:
        iterable: An iterable collection with nested lists
    """
    for i in iterable:
        if type(i) == list:
            for s in i:
                yield s
        else:
            yield i


def _get_chunks(iterable, size=1000):
    s = 0
    chunk = []
    for i in iterable:
        if s == size:
            yield chunk
            chunk = []
            s = 0
        chunk.append(i)
        s += 1
    yield chunk


def _overlap_intervals(intervals):
    if len(intervals) < 1:
        raise RuntimeError("Imposible to overlap an empty list of intervals")
    b, e = intervals[0][0], intervals[0][1]
    for i in intervals:
        b, e = max(b, i[0]), min(e, i[1])

    return b, e


def _skip_partitions(iterator, p):
    for i, v in enumerate(iterator):
        if i % p[1] != p[0]:
            continue
        yield v


def _skip_comments(iterator, char):
    for v in iterator:
        if v.startswith(char):
            continue
        yield v
