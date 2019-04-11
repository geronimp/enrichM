import sys
import re
import pickle
KEGG = '^(K\d{5})'

def ex(string):
    chunk_split = chunk.split('-')

output = sys.argv[2]
output_dict = dict()
for line in open(sys.argv[1]):
    module, description = line.strip().split('\t')
    chunks = []
    for chunk in description.split(' '):

        if '-' in chunk:
            if chunk == '--':continue
            to_keep = []

            while '-' in chunk:
                chunk_split = chunk.split('-')
                to_keep.append(chunk_split[0])
                to_assess = '-'.join(chunk_split[1:])

                if re.match(KEGG, to_assess):
                    chunk = chunk.replace('-' + to_assess[:6], "")

                else:

                    if to_assess.startswith('('):
                        chunk=chunk.replace('-' + to_assess[0:to_assess.index(')')+1], "")
        chunks.append(chunk)
    description = ' '.join(chunks) 
    if '(' in description:
        inside = False
        for i in description:
            if inside:
                count += i
            if i == '(':
                count = ''
                inside = True
            elif i == ')':

                if len(count[:-1]) == 6:
                    description= description.replace(
                        "(%s)" % count[:-1], count[:-1])
                inside = False
                count = ''

    while '  ' in description:
        description= description.replace("  ", " ")
    while description.endswith(' '):
        description = description[:-1]
    output_dict[module] = description

out_io = open(output, 'wb')
pickle.dump(output_dict, out_io)
