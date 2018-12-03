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
import time

from gendas.statistics import count
from gendas.engine import Gendas

logging.basicConfig(format='[%(name)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)

t = time.time()
gd = Gendas('data/gendas.conf')
print("Engine ready {:.3f} s".format(time.time()-t))

t = time.time()

# Count mutations that have a CADD score above 20 and are in a + strand gene
positive_above20 = count(
    gd['variants'].merge(gd['cadd'], on=['REF', 'ALT']).merge(gd['genes']).filter(
        lambda r: r['cadd']['PHRED'] > 20 and r['genes']['STRAND'] == '+'
    )
)
print("Positive: {}".format(positive_above20))
print("{:.3f} s".format(time.time()-t))
t = time.time()

# The same in a faster way
fpositive_above20 = count(
    gd['genes'].merge(gd['variants']).merge(gd['cadd'], on=['REF', 'ALT']).filter(
        lambda r: r['cadd']['PHRED'] > 20 and r['genes']['STRAND'] == '+'
    )
)
print("Fast positive: {}".format(fpositive_above20))
print("{:.3f} s".format(time.time()-t))
t = time.time()

# Count mutations that have a CADD score above 20 and are in a - strand gene
negative_above20 = count(
    gd['variants'].merge(gd['cadd'], on=['REF', 'ALT']).merge(gd['genes']).filter(
        lambda r: r['cadd']['PHRED'] > 20 and r['genes']['STRAND'] == '-'
    )
)
print("Negative:{}".format(negative_above20))
print("{:.3f} s".format(time.time()-t))
t = time.time()

# The same in a faster way
fnegative_above20 = count(
    gd['genes'].merge(gd['variants']).merge(gd['cadd'], on=['REF', 'ALT']).filter(
        lambda r: r['cadd']['PHRED'] > 20 and r['genes']['STRAND'] == '-'
    )
)
print("Fast negative: {}".format(fnegative_above20))
print("{:.3f} s".format(time.time()-t))

