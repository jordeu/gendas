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

from collections import defaultdict
from gendas.engine import Gendas
from gendas.experimental import HG19Source

# Create a Gendas engine
gd = Gendas('data/gendas.conf')
gd['hg19'] = HG19Source()


def mut_rank(gd, baserow):

    # Get scores of all the position in this gene group by tri>alt
    context = defaultdict(list)
    for r in gd['cadd'].merge(gd['hg19']):
        key = "{}>{}".format(r['hg19'][-1:1], r['cadd']['ALT'])
        context[key].append(r['cadd']['PHRED'])

    rows = []
    for m in gd['variants'].merge(gd['cadd'], on=['REF', 'ALT']).merge(gd['hg19']):

        key = "{}>{}".format(m['hg19'][-1:1], m['variants']['ALT'])
        ctx_scores = context[key]

        # Create the output row
        row = dict(baserow)
        for k, v in m['variants'].items():
            row[k] = v

        row['KEY'] = key
        row['SCORE'] = m['cadd']['PHRED']
        row['CONTEXT'] = ctx_scores

        rows.append(row)

    return rows

for res in gd.groupby(gd['exons']['GENE']).aggregate_seq(mut_rank):
    for r in res:
        print(r)
    break