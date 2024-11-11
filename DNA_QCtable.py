# Author: Huiting Xu
import argparse

import pandas as pd


def get_info_from_fastqc(path):
    df = pd.read_csv(path, sep='\t')
    assert df[df['Sample'].str.contains('_R1_')]['Total Sequences'].values.tolist() == \
           df[df['Sample'].str.contains('_R2_')][
               'Total Sequences'].values.tolist(), "Mismatch in total sequences between R1 and R2 pairs"
    df1 = df[['Sample', 'Total Sequences']]
    df1 = df1[~df1['Sample'].str.contains("R2")]
    name = df1['Sample'].str.extract(r'^(.*)_R')[0]
    df1.set_index(name, inplace=True)
    df1 = df1.drop(columns=['Sample'], axis=1)
    return df1


def main(row_fastqc_path, trimmed_fastqc_path, bwa_path, picard_path, out_name):
    row_fastqc_df = get_info_from_fastqc(row_fastqc_path)
    row_fastqc_df.rename(columns={"Total Sequences": "input_pairs"}, inplace=True)
    trimmed_fastqc_df = get_info_from_fastqc(trimmed_fastqc_path)
    trimmed_fastqc_df['trimmed_pairs'] = row_fastqc_df['input_pairs'] - trimmed_fastqc_df['Total Sequences']

    bwa_df = pd.read_table(bwa_path, sep='\t', index_col=0)
    bwa_df['mapped_pairs'] = (bwa_df['reads_properly_paired'] / 2).astype(int)
    picard_df = pd.read_table(picard_path, sep='\t', index_col=0)
    picard_df['dup_pairs'] = picard_df['READ_PAIR_DUPLICATES'].astype(int)
    summary_df = pd.concat(
        [row_fastqc_df, trimmed_fastqc_df['trimmed_pairs'], bwa_df['mapped_pairs'], picard_df['dup_pairs']], axis=1)
    summary_df.to_csv(out_name, sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate QC table from Multiqc reports')
    parser.add_argument('-m', '--multiqc_folder', type=str, help='Path to the multiqc folder')
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    main(f'{args.multiqc_folder}/multiqc_fastqc_fastqc_raw.txt',
         f'{args.multiqc_folder}/multiqc_fastqc_fastqc_trimmed.txt',
         f'{args.multiqc_folder}/multiqc_samtools_stats.txt',
         f'{args.multiqc_folder}/multiqc_picard_dups.txt', args.output)
