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

import time
import logging
import pandas as pd
from tqdm import tqdm

from gendas.engine import Gendas
from gendas.statistics import mean, max, min

logging.basicConfig(format='[%(name)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)

t = time.time()
gd = Gendas('data/gendas.conf')
print("Engine ready {:.3f} s".format(time.time()-t))

# Total genes
genes_size = len(set(gd['exons']['GENE']))

# Compute the max, mean, min of CADD score by gene.
t = time.time()

bkg_cadd_mean_by_gene = gd.groupby(gd['exons']['GENE']).aggregate({
    'MEAN': lambda gd: mean(gd['cadd']['PHRED']),
    'MAX': lambda gd: max(gd['cadd']['PHRED']),
    'MIN': lambda gd: min(gd['cadd']['PHRED'])
})

# Run the query
dataset = list(tqdm(bkg_cadd_mean_by_gene, total=genes_size))

print(pd.DataFrame.from_dict(dataset, orient='columns').set_index(['GENE']).head())
print("{:.3f} s".format(time.time()-t))


# Compute CADD score max, mean, min, total mutations, total samples mutated of the
# observed PAN mutations.
def pan_stats(gd, row, test=None):

    scores = gd['variants'].merge(gd['cadd'], on=['REF', 'ALT'])['cadd']['PHRED']
    pan_muts = gd['variants']['SAMPLE']

    row['MEAN'], row['MAX'], row['MIN'] = mean(scores), max(scores), min(scores)
    row['MUTS'], row['SMUTS'] = len(pan_muts), len(set(pan_muts))
    return row

t = time.time()
pan_cadd_mean_by_gene = gd.groupby(gd['exons']['GENE']).aggregate_seq(pan_stats, test="HOLA!")

dataset = list(tqdm(pan_cadd_mean_by_gene, total=genes_size))

print(pd.DataFrame.from_dict(dataset, orient='columns').set_index(['GENE']).head())
print("{:.3f} s".format(time.time()-t))