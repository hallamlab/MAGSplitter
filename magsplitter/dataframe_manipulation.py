import pandas as pd
from typing import List


# These functions are to combine each of the mapping and ePGDB dataframes into a list of dataframes split by MAGs.


# Since we now have the names of deleted ORFs in a clean format, we can add
# deleted ORF reactions back into the rxn dataframe
def undo_orf_removal_test(df_rxn: pd.DataFrame, duplicate_orf_map: List[List[str]]) -> pd.DataFrame:
    """
    Add in extra rows that were deleted from metapathways output, done previously to optimize usage of pathway tools
    Coordinates are restored separately from the authoritative feature table.
    :param df_rxn: rxn dataframe
    :param duplicate_orf_map: duplicate orf list of lists
    :return: rxn dataframe with extra rows
    """
    # Open the orf map file, and only choose the maps that have at least two ORFs in each row
    # The ones with only one ORF indicate that there are no duplicates of that rxn
    if df_rxn['ORF_ID'].duplicated().any():
        raise ValueError("Duplicate ORF IDs in Pathway Tools input")
    representatives = df_rxn.set_index('ORF_ID', drop=False)
    rows = [df_rxn]
    for cluster in duplicate_orf_map:
        if cluster[0] not in representatives.index:
            continue
        template = representatives.loc[[cluster[0]]]
        for member in cluster[1:]:
            rows.append(template.assign(ORF_ID=member, NAME="NAME" + member[2:]))
    expanded = pd.concat(rows, ignore_index=True)
    if expanded['ORF_ID'].duplicated().any():
        raise ValueError("Duplicate ORF IDs after expansion; check the duplicate ORF map")
    return expanded


# Combine the rxn dataframe with the orf:contig and contig:mag dataframes
def combine_rxn_contig_map(df_rxn: pd.DataFrame, df_contig_map: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the rxn dataframe with the contig map dataframe

    :param df_rxn: rxn dataframe
    :param df_contig_map: contig map dataframe
    :return: combined rxn and contig map dataframe
    """
    merged = df_rxn.merge(df_contig_map, on='ORF_ID', how='left', validate='one_to_one')
    if merged['CONTIG_ID'].isna().any():
        missing = merged.loc[merged['CONTIG_ID'].isna(), 'ORF_ID'].tolist()
        raise ValueError("Missing ORF-to-contig assignments: " + ", ".join(missing[:10]))
    return merged


def combine_rxn_mag_map(df_rxn: pd.DataFrame, df_mag_map: pd.DataFrame,
                        df_c2c_map: pd.DataFrame) -> pd.DataFrame:
    """
    Combine the rxn dataframe with the contig MAG map dataframe
    Replace non-binned ORFs with BIN ID of 'non_binned_metagenome'
    :param df_rxn: rxn dataframe
    :param df_mag_map: contig map dataframe
    :param df_c2c_map: contig to contig map dataframe
    :return: combined rxn and contig map dataframe
    """

    df = df_rxn.merge(df_c2c_map, on='CONTIG_ID', how='left')
    df = df.merge(df_mag_map, left_on='CONTIG_ORIG_ID', right_on='CONTIG_ID', how='left')
    df['BIN_ID'] = df["BIN_ID"].fillna(value='non_binned_metagenome')
    return df


def split_full_rxn_df_by_mag(df_rxn: pd.DataFrame) -> dict:
    """
    Split the full rxn dataframe into a dict of rxn dataframes, one for each MAG {BIN_ID: df}
    For all rxns that are not in a MAG, put them in a separate dataframe (metagenome not in MAG)
    :param df_rxn: full rxn dataframe
    :return: list of rxn dataframes, one for each MAG
    """
    mags = df_rxn['BIN_ID'].unique()
    df_rxn_list = {}
    for mag in mags:
        df_rxn_list[mag] = df_rxn[df_rxn['BIN_ID'] == mag]

    return df_rxn_list


# Misc functions
def sample_name_grabber(rxn_df: pd.DataFrame) -> pd.DataFrame:
    """
    Grab the sample name from the contig dataframe.  This is under the assumption
    that the sample name is the same for all contigs in a MAG.
    :param rxn_df: 
    :return: 
    """
    # end index is one character after the last character in the string
    full_str = rxn_df['CONTIG_ID'].iloc[0]
    end_idx = len(full_str)
    # If we only have the sample name without contig ID, turn <is_only_sample> to True
    is_only_sample = False
    while end_idx > 0 and is_only_sample is False:
        if full_str[end_idx - 1] == '_':
            is_only_sample = True
            end_idx -= 1
        else:
            end_idx -= 1
    return full_str[:end_idx]


if __name__ == '__main__':
    pass
