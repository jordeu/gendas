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

import collections
import os
import struct
import zlib

from gendas.tabix.constants import *


class BlockReader:
    """
    Class to randomly decompress a single block in a tabix block-gzipped file
    """

    def __init__(self, filename, cache_size=1000):
        """

        Args:
            filename: Path to a tabix block-gzipped file
            cache_size: Maximum number of block to keep uncompressed in the cache
        """

        # Check that the data file exists
        if not os.path.isfile(filename):
            raise IOError("The path '{}' is not a file".format(filename))

        self.__filename = filename
        self.__data = open(self.__filename, "rb")
        self.__header = None

        self.__partial_line_ends = {}
        self.__blocks_cache = LRUCache(cache_size)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        self.__header = None
        if self.__data is not None:
            self.__data.close()

    def read(self, block):
        """
        Uncompress a block and return the contents uncompressed

        Args:
            block: Block number

        Returns:
            A list with all the text lines in this block

        """
        block_address = (block >> SHIFT_AMOUNT) & ADDRESS_MASK
        block_offset = int(block & OFFSET_MASK)

        block_length, block_content = self.__read_block_bytes(block_address)

        # Remove offset
        block_content = block_content[block_offset:]

        # Cache of partial line ends
        self.__partial_line_ends[block_address] = str(block_content[:block_offset])

        # Check if the last line is partial
        is_partial_last_line = block_content[len(block_content) - 1] != '\n'

        # Split lines
        lines = block_content.splitlines()

        # Complete last line
        if is_partial_last_line:

            next_block_address = block_address + block_length
            if next_block_address not in self.__partial_line_ends:
                # Force to read next block
                next_block_length, next_block_content = self.__read_block_bytes(next_block_address)
                next_block_offset = next_block_content.find('\n')
                self.__partial_line_ends[next_block_address] = str(next_block_content[:next_block_offset])

            complete_line = lines.pop() + self.__partial_line_ends[next_block_address]
            lines.append(complete_line)

        return lines

    def __read_block_bytes(self, block_address):

        # First check cache
        if self.__blocks_cache.has_key(block_address):
            block = self.__blocks_cache.get(block_address)
            return block[0], block[1]

        self.__data.seek(block_address)

        # Read the block header
        header = self.__data.read(BLOCK_HEADER_LENGTH)
        # Extract compressed block length
        block_compressed_length = struct.unpack_from("H", header, offset=BLOCK_LENGTH_OFFSET)[0] + 1

        # Read compressed block
        block_compressed = header + self.__data.read(block_compressed_length - BLOCK_HEADER_LENGTH)

        # Decompress
        lines = zlib.decompress(block_compressed, 15 + 32)

        # Decode string
        lines = lines.decode("utf-8")

        self.__blocks_cache.set(block_address, (block_compressed_length, lines))

        return block_compressed_length, lines

    def header(self):
        """
        Returns: The first line of the file
        """
        if self.__header is None:
            self.__header = self.read(0)[0]

        return self.__header


class LRUCache:
    """
    A simple in-memory cache implementation
    """

    def __init__(self, capacity):
        self.capacity = capacity
        self.cache = collections.OrderedDict()

    def get(self, key):
        try:
            value = self.cache.pop(key)
            self.cache[key] = value
            return value
        except KeyError:
            return -1

    def set(self, key, value):
        try:
            self.cache.pop(key)
        except KeyError:
            if len(self.cache) >= self.capacity:
                self.cache.popitem(last=False)
        self.cache[key] = value

    def has_key(self, key):
        return key in self.cache
