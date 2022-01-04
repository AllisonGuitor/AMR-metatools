#!/usr/bin/env python3
import sys

assert sys.argv[1] is not None
seqs = []
with open(sys.argv[1], 'r') as open_handle:
    item = None
    for line in open_handle.readlines():
        if line[0] == '>':
            if item is not None:
                item['len'] = len(item['seq'])
                seqs.append(item)
            item = {'name': line.strip().lstrip('>'), 'seq': ""}
        else:
            item['seq'] += line.strip()
    item['len'] = len(item['seq'])
    seqs.append(item)

lens = [i['len'] for i in seqs]

total = sum(lens)
sorted_lens = sorted(lens, reverse=True)
running_sum = 0.0
running_frac = 0.0
longest = 0
for index, length in enumerate(sorted_lens):
    if length > longest:
        longest = length
    running_sum += float(length)
    running_frac = (running_sum / total)
    if running_frac >= 0.5:
        print(
            sys.argv[1] + ': ' + str(len(lens)) + ' contigs, total: ' +
            str(total) + ' bp'
        )
        print('Largest contig: {0} bp'.format(longest))
        print(
            '(' + str(index - 1) + ') N50 - 1:' + str(sorted_lens[index - 1]) +
            ' bp'
        )
        print(
            '(' + str(index) + ') N50    :' + str(sorted_lens[index]) + ' bp'
        )
        print(
            '(' + str(index + 1) + ') N50 + 1:' + str(sorted_lens[index + 1]) +
            ' bp'
        )
        break
if len(sys.argv) == 3:
    print(sys.argv[2])
    max_len = int(sys.argv[2])
    with open('{0}.filtered'.format(sys.argv[1]), 'w') as outh:
        for s in seqs:
            if len(s['seq']) > max_len:
                outh.write('>{0}\n{1}\n'.format(s['name'], s['seq']))
