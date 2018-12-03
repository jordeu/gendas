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

# Block gzip constants
BLOCK_HEADER_LENGTH = 18
BLOCK_LENGTH_OFFSET = 16
BLOCK_FOOTER_LENGTH = 8

# Virtual file pointer constants
SHIFT_AMOUNT = 16
OFFSET_MASK = 0xffff
ADDRESS_MASK = 0xFFFFFFFFFFFF

# Compression level
COMPRESSION_LEVEL = 5

# Gzip overhead is the header, the footer, and the block size (encoded as a short).
GZIP_OVERHEAD = BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH + 2

# If Deflater has compression level == NO_COMPRESSION, 10 bytes of overhead (determined experimentally).
NO_COMPRESSION_OVERHEAD = 10

# Push out a gzip block when this many uncompressed bytes have been accumulated.
# This size is selected so that if data is not compressible, if Deflater is given
# compression level == NO_COMPRESSION, compressed size is guaranteed to be <= MAX_COMPRESSED_BLOCK_SIZE.
DEFAULT_UNCOMPRESSED_BLOCK_SIZE = 64 * 1024 - (GZIP_OVERHEAD + NO_COMPRESSION_OVERHEAD)
