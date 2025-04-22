"Developed by Bernardo Almeida, version 2, April 2025"

#IMPORTS
import gzip
import shutil
import os
import pandas as pd
import csv
import re


#MERGE MULTIPLE FASTQ.GZ FILES IN A SINGLE FILE

# IMPUT PATH: Path to the folder with the fastq.gz files
folder_of_gz_files_1 = r'C:\Users\Utilizador\Desktop\Universidade\Mestrado MCB\_Dissertacao\Sequenciações\MiniON\MinION 04Abr\fastq_fail\barcode04'
# OUTPUT PATH: Path to the folder we want to save our merged file
output_file_1 = r'C:\Users\Utilizador\Desktop\Universidade\Mestrado MCB\_Dissertacao\Sequenciações\MiniON\MinION 04Abr\fastq_fail\Mergedbarcode04_fastq_fail.fastq.gz'

def fastq_gz_merge(folder_of_gz_files_1, output_file_1):
    with gzip.open(output_file_1, 'wb') as merged_file:
        for file in os.listdir(folder_of_gz_files_1):
            if file.endswith('.fastq.gz'):
                with gzip.open(os.path.join(folder_of_gz_files_1, file), 'rb') as f_in:
                    shutil.copyfileobj(f_in, merged_file)
    print(f"All files have been combined in {output_file_1}")

#fastq_gz_merge(folder_of_gz_files_1, output_file_1)



# FORMATTING OF THE .TABULAR FILE CREATED ON USEGALAXY.ORG FROM THE BAM FILE:
    # Rows: Different positions of the sequence
    # Column 1: Contig identifier
    # Column 2: Nucleotide position
    # Column 3: Reference base
    # Column 4: Read coverage
    # Column 5: Read alignment "Pileup"
    # Column 6: Base quality encoding
    # Column 7: Mapping quality

# Mapping multiple-base ambiguities to IUPAC codes
ambiguity_codes = {
    'AG': 'R',
    'CT': 'Y',
    'GT': 'K',
    'AC': 'M',
    'CG': 'S',
    'AT': 'W',
    'CGT': 'B',
    'ATG': 'D',
    'ACT': 'H',
    'ACG': 'V',
    'ACGT': 'N'
}

#file_name ends in .tabular
file_name = r'C:\Users\Utilizador\Desktop\Universidade\Mestrado MCB\_Dissertacao\Sequenciações\MiniON\MinION 5Nov\MinION 5Nov (9h30 8Nov)\PCV2\Tabular_MinION_5Nov_PCV2_GenomeDetective.tabular'
#output_file ends in .html
output_file = r'C:\Users\Utilizador\Desktop\Universidade\Mestrado MCB\_Dissertacao\Sequenciações\MiniON\MinION 5Nov\MinION 5Nov (9h30 8Nov)\PCV2\bam_GenomeDetective_MinION 5Nov barcode2_PCV2.html'

def convert_to_iupac(seq):
    if len(seq) == 1:
        return seq
    return ambiguity_codes.get(''.join(sorted(seq)), 'N')  # Default to N

def clean_pileup(pileup):
    patterns = list(re.finditer(r'[+-](\d+)', pileup))
    result = ''
    index = 0
    for pattern in patterns:
        start, finish = pattern.span()
        number = int(pattern.group(1))
        end_remove = finish + number
        result += pileup[index:start]
        index = end_remove
    result += pileup[index:]
    return result

# Auxiliary function to count characters
def count_characters(pileup):
    insertions = pileup.count('+')
    deletions = pileup.count('-')
    pileup = clean_pileup(pileup)  # clean insertions e deletions
    counts = {
        'Aa': pileup.count('A') + pileup.count('a'),
        'Cc': pileup.count('C') + pileup.count('c'),
        'Tt': pileup.count('T') + pileup.count('t'),
        'Gg': pileup.count('G') + pileup.count('g'),
        'Deletions': deletions,
        'Insertions': insertions,
        }
    # Couts characters different from: A, a, C, c, T, t, G, g, -, +
    nucleotideo_total = counts['Aa'] + counts['Cc'] + counts['Tt'] + counts['Gg']
    counts['Others'] = len(pileup) - nucleotideo_total
    return pd.Series(counts)

def determine_consensus(row):
    # Dictionary for counts and corresponding letters
    base_counts = {
        'A': row['Aa'],
        'C': row['Cc'],
        'T': row['Tt'],
        'G': row['Gg']
    }
    # Skip if all counts are zero (edge case)
    if sum(base_counts.values()) == 0:
        return 'N'
    max_count = max(base_counts.values())
    consensus_bases = []
    for base, count in base_counts.items():
        # Include the base if its count is ≥ 80% of the max count
        if count >= 0.80 * max_count and count > 0:
            consensus_bases.append(base)
    return ''.join(sorted(consensus_bases))  # Sort for consistency

# Function to style columns
def style_specific_columns(s):
    color_map = {
        'Aa': 'background-color: #3CFF14',            # Light Green
        'Cc': 'background-color: #30F5E6',            # Light Blue
        'Tt': 'background-color: #F93D15',            # Red
        'Gg': 'background-color: #CE19F6',            # Light Purple
        'Deletions': 'background-color: #BFB917',     # Dark Yellow
        'Insertions': 'background-color: #FAF68C'     # Light Yellow
    }
    return [color_map.get(s.name, '')] * len(s)

# Function to style column "Consensus Seq"
def style_consensus(value):
    color_map = {
        'A': 'background-color: #3CFF14',             # Light Green
        'C': 'background-color: #30F5E6',             # Light Blue
        'T': 'background-color: #F93D15',             # Red
        'G': 'background-color: #CE19F6',             # Light Purple
        'Deletions': 'background-color: #BFB917',     # Dark Yellow
        'Insertions': 'background-color: #FAF68C'     # Light Yellow
    }
    # If multiple letters have the same highest count, leave the cell with a white background
    if len(value) > 1:
        return 'background-color: #FFFFFF'  # White
    return color_map.get(value, '')

def read_bam(file_name, output_file):
    try:
        data = pd.read_table(file_name, header=None, sep=r'\s+', quoting=csv.QUOTE_NONE, engine='python')
        data = data.drop(columns=[0,2,5])
        data = data.rename(columns={1: 'Nt. position', 3: 'Read Coverage', 4: 'Read alignment "Pileup"'})
        # Applies the counting function to each line and adds the new columns
        character_counts = data['Read alignment "Pileup"'].apply(count_characters)
        data = pd.concat([data[['Nt. position', 'Read Coverage']], character_counts, data[['Read alignment "Pileup"']]], axis=1)
        # Adds the "Consensus Seq" column using the determine_consensus function
        data['Consensus Seq'] = data.apply(determine_consensus, axis=1)
        data = data[['Nt. position', 'Read Coverage', 'Aa', 'Cc', 'Tt', 'Gg', 'Deletions', 'Insertions', 'Others', 'Consensus Seq', 'Read alignment "Pileup"']]
        # Generate Consensus Sequence as a string
        consensus_sequence = ''.join([convert_to_iupac(seq) for seq in data['Consensus Seq']])
        # Save consensus sequence to .txt file
        txt_output_file = output_file.replace(".html", "_ConsensusSequence.txt")
        with open(txt_output_file, 'w') as f:
            f.write(consensus_sequence)
        # Applying specific style for the "Consensus Seq" column
        styled_data = data.style
        styled_data = styled_data.map(style_consensus, subset=['Consensus Seq'])
        # Define specific styles for column alignment
        styled_data.set_table_styles([
            {'selector': 'thead th', 'props': [('text-align', 'left')]}, 
            {'selector': 'td:nth-child(1)', 'props': [('text-align', 'left')]},     # Nt. position
            {'selector': 'td:nth-child(2)', 'props': [('text-align', 'left')]},     # Read Coverage
            {'selector': 'td:nth-child(10)', 'props': [('text-align', 'center')]},  # Consensus Seq
            {'selector': 'td:nth-child(11)', 'props': [('text-align', 'left')]}     # Read alignment "Pileup"
        ])

        html = styled_data.hide(axis="index").to_html(index=False, justify='left')
        custom_css = """
        <style>
            .data.col2 {
                background-color: #3CFF14;
            }
             .data.col3 {
                background-color: #30F5E6;
            }
             .data.col4 {
                background-color: #F93D15;
            }
             .data.col5 {
                background-color: #CE19F6;
            }
             .data.col6  {
                background-color: #BFB917;
            }
             .data.col7  {
                background-color: #FAF68C;
            }
        </style>
        """

        # Combine the CSS and the HTML
        final_html = custom_css + html
        # Converte o DataFrame estilizado para HTML e salva no arquivo
        with open(output_file, "w") as text_file:
            text_file.write(final_html)
        print("✅ HTML file created successfully!")
        print("✅ TXT file created successfully!")
    except Exception as e:
        print(f"❌ An error occurred: {e}")

read_bam(file_name, output_file)
