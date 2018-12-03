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
import numpy as np
import pandas as pd

from gendas.engine import Gendas


def oncodrive_fml(gd, sampling=100):

    cadds_gene = list(gd['cadd']['PHRED'])
    cadds_observed = list(gd['variants'].merge(gd['cadd'], on=['REF', 'ALT'])['cadd']['PHRED'])

    if len(cadds_observed) == 0:
        return None

    background = np.array([np.mean(np.random.choice(cadds_gene, size=len(cadds_observed))) for a in range(sampling)])
    obs = len(background[background >= np.mean(cadds_observed)])

    return max(1, obs) / sampling


logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', datefmt='%H:%M:%S', level=logging.INFO)
gd = Gendas('data/gendas.conf')

t = time.time()
pvalues = gd.groupby(gd['exons']['GENE']).aggregate({
    'PVALUE': oncodrive_fml
})

df = pd.DataFrame.from_dict(pvalues, orient='columns').set_index(['GENE'])
df = df.sort_values('PVALUE', ascending=True)
print(df.head(100))
print("{:.3f} s".format(time.time()-t))
