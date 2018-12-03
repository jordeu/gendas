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
    Basic statistic algorithms computed over an iterable collection without loading
    all the values into memory.
"""

import builtins
import itertools
import statistics


def peek(iterable):
    """
    Get the first value of an iterator without iterating two times an iterator

    Args:
        iterable: An iterable

    Returns:
        The first value and an iterator to iterate all the values included the first one.

    """
    try:
        first = next(iterable)
    except StopIteration:
        return None, []
    return first, itertools.chain([first], iterable)


def empty(fn, values, default=None):
    """
    Returns 'default' value if 'values' list is empty. Otherwise calls 'fn' with the 'values' list.

    Args:
        fn: The function to call if the list is not empty
        values: Values to pass to the 'fn' function
        default: The value to return if the list is empty

    Returns:
        'default' if 'values' is empty, 'fn(values)' otherwise.

    """
    if values is None:
        return default

    first, iterator = peek(iter(values))

    if first is None:
        return default

    return fn(iterator)


def mean(values):
    """
    Computes the mean value

    Args:
        values: An iterable collection

    Returns:
        The mean or None if values it's empty

    """
    return empty(statistics.mean, values)


def min(values):
    """
    Computes the minimum value

    Args:
        values: An iterable collection

    Returns:
        The minimum value or None if values it's empty

    """
    return empty(builtins.min, values)


def max(values):
    """
    Computes the maximum value

    Args:
        values: An iterable collection

    Returns:
        The maximum value or None if values it's empty

    """
    return empty(builtins.max, values)


def count(iterator):
    """
    Counts the length of an iterator iterating it

    Args:
        iterator: An iterable collection

    Returns:
        How many elements you have in the iterator

    """
    return sum(1 for i in iterator)
