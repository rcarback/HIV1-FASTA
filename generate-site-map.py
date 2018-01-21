#!/usr/bin/python

import click
import json
import yaml

from collections import defaultdict

def parse_fasta(pf):
    """The first line in fasta files are a description that starts
    with >, the second is the peptide sequence.

    We assume lines that start with > are description lines and the
    line after is the sequence, and return a map from desc ->
    sequence"""
    seqs = dict()
    cur_desc = ''
    for line in pf:
        if line.startswith('>'):
            cur_desc = line.rstrip()[1:]
        else:
            seqs[cur_desc] = line.rstrip()
    return seqs

def parse_epitope_csv(ecsv):
    for line in ecsv:
        peptide,protein,startend,epitope = line.rstrip().split(',')
        start, end = startend.split('-')
        yield {
            'peptide': peptide.upper(),
            'protein': protein.upper(),
            'start': start,
            'end': end,
            'epitope': epitope.upper()
        }

def find_proteins(seqs, peptide):
    """Find the peptide in the sequences db and return the
    corresponding protein names"""
    proteins = []
    for seq, prots in seqs.items():
        if peptide in seq:
            proteins.extend(prots)

    return proteins

@click.command()
@click.option('-m', '--protein-map', type=click.File('rb'),
              help='Filename to proteins map')
@click.option('-s', '--sequence-csv', type=click.File('rb'),
              help='Filename to sequence, protein, site, epitope csv')
def main(protein_map, sequence_csv):
    # Read the file map to the protein names
    pmap = yaml.load(protein_map)

    # Open each fasta file, and load the entire sequence into a map
    # for the protein names, mapping full sequence to the list of proteins
    sequences = dict()
    for fastafn, proteins in pmap.items():
        with open(fastafn, 'rb') as fasta:
            for desc,seq in parse_fasta(fasta).items():
                sequences[seq] = [p.upper() for p in proteins]

    # Read the epitope map csv (peptide,protein,start-end,epitope)
    emap = parse_epitope_csv(sequence_csv)

    # For each site&peptide, mark the corresponding epitope and proteins
    sitemap = defaultdict(list)
    peptideinfo = dict()
    for e in emap:
        start = int(e['start'])
        end = int(e['end'])
        # Multiple proteins map to a given peptide
        proteins = find_proteins(sequences, e['peptide'])
        if not proteins:
            print('Could not find peptide,protein in sequences db: '
                  '{},{},{},{}-{}'.format(e['peptide'], e['protein'],
                                          e['epitope'],
                                          e['start'], e['end']))
        # Also include the protein from the epitope mapping
        proteins.append(e['protein'])
        proteins = list(set(proteins))

        key = '{},{},{}'.format(e['peptide'], e['protein'], e['epitope'])
        peptideinfo[key] = {
            'proteins': [p for p in proteins],
            'epitope': e['epitope'],
            'start': e['start'],
            'end': e['end']
            # 'peptide': e['peptide']
        }

        for site in range(start, end + 1):
            sitemap[site].append(key)

    # Save the site mappings
    with open('sitemapdb.json', 'wb') as out:
        json.dump(dict(sitemap), out)
    with open('peptideinfo.json', 'wb') as out:
        json.dump(peptideinfo, out)

if __name__ == '__main__':
    main()
