# 
# This file is part of the GSprofiler distribution (https://github.com/xxxx or http://xxx.github.io).
# Copyright (c) 2015 Liviu Ionescu.
# 
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import requests
import pandas as pd
import argparse
from numpy import log10
import matplotlib.pyplot as plt
import csv

""" Small program to download gProfiler term enrichment analysis"""


def args():
    """Argument parser"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile',
                        help='File with the list of genes to analyse')
    parser.add_argument('organism',
                        help='''Organism name according to gProfiler list.
                        Example: mmusculus, hsapiens, ggallus,
                        dmelanogaster.''')
    parser.add_argument('--output', '-o', 
    	help="""Output file""")
    parser.add_argument('--method', '-m', default='g_SCS',
    	help="""False discovery method: 'g_SCS', 'bonferroni' and 'fdr'""")
    parser.add_argument('--threshold', '-t', default=0.05,
    	help="""p-value threshold""")
    parser.add_argument('--underrepresented', '-u', default=False, action='store_true',
    	help="""If specified, returns the under represented terms (isntead of the over represented) """)
    return parser.parse_args()

def gprofiler(namelist, organism, user_threshold=0.05, method='g_SCS',
              measure_underrepresentation=False, background=None,
              simple_out=0):
    """Run gProfiler using POST api with a json query body

    Returns a pandas DataFrame with the result

    methods = 'g_SCS', 'bonferroni' and 'fdr'
    simple_out = 3 levels of output 0|1|2
    background = TO DO!!!!!!!!
    """
    if type(namelist) is not list:
        namelist = list(namelist)
    r = requests.post(
        url='https://biit.cs.ut.ee/gprofiler/api/gost/profile/',
        json={
            'organism': organism,
            'query': namelist,
            'user_threshold': user_threshold,
            'measure_underrepresentation': measure_underrepresentation,
        }
    )
    df = pd.DataFrame(r.json()['result'])
    # manually selecting and sorting output
    cols = ['source', 'native', 'name', 'p_value', 'description', 'query',
            'significant']
    extend = ['term_size', 'query_size', 'intersection_size',
              'effective_domain_size', 'intersections', 'parents']
    # Unkown columns, there is no direct explanation.
    wtf_cols = ['goshv', 'group_id', 'precision', 'recall',
                'source_order']
    assert simple_out in [0, 1, 2], 'Simple_out wrong value'
    if simple_out == 0:
        df = df[cols]
    elif simple_out == 1:
        df = df[cols + extend]
    elif simple_out == 2:
        df = df[cols + extend + wtf_cols]
    return df


args = args()

if args.output == None:
	args.output = args.infile + '.gprofiler'

query = []

with open(args.infile) as fp:
    reader = csv.reader(fp, delimiter='\t')
    # read first column with Ensembl IDs and skip header line
    query = list(zip(*reader))[0][1:]

print('Arguments')
print(args)

result = gprofiler(query, args.organism,
	user_threshold=args.threshold,
	method=args.method,
	measure_underrepresentation=args.underrepresented)

result.to_csv(args.output, sep='\t')


print('Some plots')
sources = result.source.unique()

result = result.set_index('name')

for source in sources:
    subdf = (-log10(result[result.source == source].p_value))
    # Top 10
    #subdf = subdf.iloc[:10]
    # All significant processes
    subdf = subdf.iloc[:]
    y = 0.2 * len(subdf)
    plt.figure(figsize=(10,y))
    # subdf.plot.barh(color=source_colors[source])
    subdf.plot.barh()
    plt.gca().invert_yaxis()
    plt.ylabel('process')
    plt.xlabel('-log10(p-value)')
    plt.title(source)
    plt.savefig(args.output + source + '.svg', bbox_inches = "tight")
    plt.savefig(args.output + source + '.png', bbox_inches = "tight")


print('[DONE] Be happy!!!')
