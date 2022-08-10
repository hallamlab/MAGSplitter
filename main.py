from input_to_df import *
from dataframe_manipulation import *
from ptools_writing import *

'''
#####
Inputs
As of now, the program searches for these files within the MAGSplitter folder.  
This will be changed into a user-defined location for inputs and outputs in the future.
0.pf: pathway tools initial input from metapathways outputs
(found in /ptools/)

<sample_name>.ORF_annotation_tabel.txt: map of ORFS and contigs for an environment metagenome, outputted by metapathways 
(found in /results/annotation_table/)

orf_map.txt: map of duplicated ORFS within a metagenome, with the first ORF being the intact ORF, 
and the rest being the duplicates
(found in /results/annotation_table/)

config_info.tsv: map of contigs to MAGs created through the WGS pipeline binning process.  
(found in /binning/results/greedy)

#####
Outputs
The program will create a folder called "results" in the MAGSplitter folder.
This folder will contain folders for each MAG in the MAGSplitter folder.
Each MAG folder will contain the following files:
0.pf: pathway tools input file for pathway tools
genetic-elements.dat: pre-written file, which I should ask Ryan about
organism-params.dat: which I will also ask Ryan about
pathologic.log: pre-written file, which I should ask Ryan about
<sample_name>.dummy.txt: Dummy output file
'''
DIR_PATH = os.path.dirname(os.path.realpath(__file__))

# File name inputs, which are default set to the example folder
pf_file = 'example/0.pf'
orf_map_file = 'example/orf_map.txt'
orf_contig_map_file = "example/GAPP-5498e568-6918-4000-b27e-dbeff35eeee7.ORF_annotation_table.txt"
contig_mag_map_file = "example/contig_info.tsv"


def main():
    # Import MP/WGS pipeline outputs
    df_pf = convert_pl_input_to_rxn_df(pf_file)
    df_duplicate_orf_map = convert_duplicate_orf_map_to_list(orf_map_file)
    orf_contig_map = convert_orf_contig_map_to_df(orf_contig_map_file)
    sample_name = sample_name_grabber(orf_contig_map)
    contig_mag_map = convert_contig_mag_map_to_df(contig_mag_map_file)
    print("file imports done")
    df_pf = undo_orf_removal(df_pf, df_duplicate_orf_map)
    print("ORF de-remover done")
    df_pf_with_contig = combine_rxn_contig_map(df_pf, orf_contig_map)
    df_pf_with_mag = combine_rxn_mag_map(df_pf_with_contig, contig_mag_map)
    df_pf_split = split_full_rxn_df_by_mag(df_pf_with_mag)
    ptools_folder_creator(DIR_PATH, sample_name, df_pf_split)
    print("ptools folder creation done")

"""   
# Parse arguments
    parser = argparse.ArgumentParser(description='Convert metapathways output to ePGDB readable format')
    parser.add_argument('-i', '--input', help='metapathways output name (normally 0.pf)', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-m', '--metacyc', help='metacyc accession map', required=True)
    parser.add_argument('-c', '--contig_mag_map', help='contig mag map', required=True)
    parser.add_argument('-e', '--orf_map', help='orf map', required=True)
    args = parser.parse_args()
"""
if __name__ == '__main__':
    main()
