"""Restore each expanded feature's own coordinates before writing MAG inputs."""
from pathlib import Path
import re

import pandas as pd


def restore_coordinates(features, annotation_table, feature_table=None):
    """Use MP's pre-compaction feature table; never inherit template coordinates.

    MP stores start <= end and an explicit strand. Pathologic PF uses descending
    STARTBASE/ENDBASE for the reverse strand. All checks precede output creation.
    """
    if feature_table is None:
        annotation_table = Path(annotation_table)
        suffix = '.ORF_annotation_table.txt'
        if not annotation_table.name.endswith(suffix):
            raise ValueError('Cannot locate coordinates: supply --feature_table with the MP .ptinput.tsv')
        feature_table = annotation_table.with_name(
            annotation_table.name[:-len(suffix)] + '.ptinput.tsv')
    source = Path(feature_table)
    if not source.is_file():
        raise ValueError(f'Missing coordinate table {source}; supply --feature_table with the MP .ptinput.tsv')
    table = pd.read_csv(source, sep='\t', dtype=str, keep_default_na=False)
    required = {'orf_id', 'seqname', 'start', 'end', 'strand', 'contig_length'}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f'{source}: missing columns: {", ".join(sorted(missing))}')
    if table['orf_id'].eq('').any() or table['orf_id'].duplicated().any():
        raise ValueError(f'{source}: empty or duplicate feature IDs')
    table = table.set_index('orf_id')
    result = features.copy()
    starts, ends = [], []
    for row in result.itertuples():
        feature_id = row.ORF_ID[3:]
        if feature_id not in table.index:
            raise ValueError(f'{source}: missing coordinates for {feature_id}')
        feature = table.loc[feature_id]
        if not feature.seqname or feature.seqname != row.CONTIG_ID:
            raise ValueError(f'{feature_id}: coordinate-table contig does not match ORF annotation table')
        values = [feature.start, feature.end, feature.contig_length]
        if any(re.fullmatch(r'[0-9]+', value) is None for value in values):
            raise ValueError(f'{feature_id}: coordinates and contig length must be positive integers')
        start, end, length = map(int, values)
        if not 1 <= start <= end <= length:
            raise ValueError(f'{feature_id}: coordinates {start}..{end} outside contig length {length}')
        if feature.strand not in ('+', '-'):
            raise ValueError(f'{feature_id}: invalid strand {feature.strand!r}')
        if feature.strand == '-':
            start, end = end, start
        starts.append(f'STARTBASE\t{start}')
        ends.append(f'ENDBASE\t{end}')
    result['STARTBASE'] = starts
    result['ENDBASE'] = ends
    return result
