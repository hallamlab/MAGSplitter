from typing import List
import pandas as pd

# FUNCTIONS TO IMPORT PF FORMAT FILES


def convert_rxn_row_list_format(rxn_row: str) -> List[str]:
    """
    Convert a row in a reaction database to a list of two strings: the type of information, and the information itself
    :param rxn_row: reaction entry from rxn tools output
    :return: list of two strings: one of the type of information, and the information itself
    """
    end_idx = 1 # We add to end index until it represents the first tab in a line.
    while not rxn_row[end_idx] == '\t':
        end_idx += 1
    return [rxn_row[:end_idx], rxn_row]


def convert_pl_input_to_rxn_df(pl: str) -> pd.DataFrame:
    """
    Use the name from pl to open a pathway tools file of the same name
    Read the specific file format of pathway tools output from a Metapathways output.
    :param pl: pathway tools input file name
    :return: rxn dataframe
    """

    # Convert pl into a list with each row as an item
    with open(pl) as f:
        pl_lines = [line.rstrip() for line in f]

    # Used a two pointer style algo to determine where each rxn starts and ends,
    # then input into a dict accumulator

    test_acc = 0
    start_idx = 0  # Starting index of rows in a single rxn
    rxn_acc = []  # Accumulator that we input each single rxn into
    while start_idx < len(pl_lines):
        end_idx = start_idx  # Ending index of rows in a single rxn
        temp_rxn = {'ec': [], 'metacyc': []}  # Default values for ec and ko, # since these do not appear in every rxn
        while pl_lines[end_idx] != "//":  # Check every line to see if it's an end of a rxn entry
            entry_list_format = convert_rxn_row_list_format(pl_lines[end_idx])
            # Special case for rows signifying ec or metacyc accessions,
            # since there can be [0,inf] of these in a rxn
            if entry_list_format[0] == "EC":
                temp_rxn['ec'].append(entry_list_format[1])
            elif entry_list_format[0] == "METACYC":
                temp_rxn['metacyc'].append(entry_list_format[1])
            else:
                temp_rxn[entry_list_format[0]] = entry_list_format[1]
            end_idx += 1
        start_idx += (end_idx + 1 - start_idx)  # Update start_idx to the next rxn start
        rxn_acc.append(temp_rxn)
    df_pf = pd.DataFrame(rxn_acc)
    # Update column names to more readable and pandas usable names
    return df_pf.rename(columns={"ID": "ORF_ID",
                                 "PRODUCT-TYPE": "PRODUCT_TYPE"})


# FUNCTION TO IMPORT ORF CONTIG ANNOTATION MAP
def convert_orf_contig_map_to_df(orf_map: str) -> pd.DataFrame:
    """
    Read the ORF annotation table and change the name of ORF ID to "ID" for easier joins
    :param orf_map: ORF annotation table
    :return: ORF annotation table as a dataframe
    """

    df = pd.read_csv(orf_map, sep='\t')
    df = df.rename(columns={'# ORF_ID': 'ORF_ID'})
    df["ORF_ID"] = "ID\t" + df["ORF_ID"]
    return df[["ORF_ID", "CONTIG_ID"]]  # Return only the two columns we need


def convert_contig_contig_map_to_df(contig_map: str):
    """
    Read the contig to contig mapping file produced during the preprocessing
    step of metapathways. This will allow for MP internal IDs to be mapped
    back to the original contig IDs.
    :param contig_map: Contig mapping file
    :return: Contig mapping dataframe
    """

    df = pd.read_csv(contig_map, sep='\t', header=None)
    df.columns = ["CONTIG_ID", "CONTIG_ORIG_ID", "CONT_LEN"]
    return df[["CONTIG_ID", "CONTIG_ORIG_ID"]]

def convert_contig_mag_map_to_df(contig_mag_map: str) -> pd.DataFrame:
    """
    Read the contig mag map, and extract the CONTIG_ID and BIN_ID as a dataframe
    Input should be a TSV with no header where the CONTIG_ID is column 1 and BIN_ID
    is column 2.
    :param contig_mag_map: contig mag map
    :return: contig mag map as a dataframe
    """

    try:
        df = pd.read_csv(contig_mag_map, sep='\t', header=None)
    except: # If map is empty or wrong format create empty DF, should be improved
        df = pd.DataFrame(columns=['CONTIG_ID', 'BIN_ID'])
    df.columns = ['CONTIG_ID', 'BIN_ID']
    return df[['BIN_ID', 'CONTIG_ID']]


def convert_duplicate_orf_map_to_list(duplicate_orf_map: str) -> List[List[str]]:
    """
    Read the duplicate ORF map, and create a list of lists,
    with each item after the first in a nested list showing deleted duplicates of first item
    :param duplicate_orf_map: duplicate orf map
    :return: duplicate orf map as a dataframe
    """

    with open(duplicate_orf_map) as f:
        orf_map = [(line.rstrip().split("\t"))
                   for line in f
                   if ('\t' in line)]
    # Prepend each orf with the required rxn tools syntax
    for rxn_idx in range(len(orf_map)):
        for orf_idx in range(len(orf_map[rxn_idx])):
            orf_map[rxn_idx][orf_idx] = 'ID\t' + orf_map[rxn_idx][orf_idx]
    return orf_map
if __name__ == '__main__':
    pass
