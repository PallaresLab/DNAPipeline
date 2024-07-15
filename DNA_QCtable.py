import argparse
import pandas as pd


def main(general_path, bwa_path, picard_path, out_name):
    df = pd.read_table(general_path)
    df_bwa = pd.read_table(bwa_path)
    df_picard = pd.read_table(picard_path)
    df_R1 = df[df['Sample'].str.contains('R1')].reset_index(drop=True)
    df_R2 = df[df['Sample'].str.contains('R2')].reset_index(drop=True)
    col1 = 'FastQC (raw)_mqc-generalstats-fastqc_raw-total_sequences'
    assert (df_R1[col1].dropna() == df_R2[col1].dropna()).all()
    col2 = 'FastQC (trimmed)_mqc-generalstats-fastqc_trimmed-total_sequences'
    assert (df_R1[col2].dropna() == df_R2[col2].dropna()).all()
    df = df[~df['Sample'].str.contains('R1|R2')].reset_index(drop=True)

    df_table = pd.DataFrame()
    df_table['Sample'] = df['Sample']
    df_table['input_pairs'] = df_R1[col1].dropna().reset_index(drop=True).astype(int)
    df_table['aftertrim_pairs'] = df_R1[col2].dropna().reset_index(drop=True).astype(int)
    df_table['mapped_pairs'] = (df_bwa['reads_properly_paired'] / 2).astype(int)
    df_table['dup_pairs'] = df_picard['READ_PAIR_DUPLICATES'].astype(int)
    df_table.to_csv(out_name, index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate QC table from Multiqc reports')
    parser.add_argument('-i','--multiqc_folder', type=str, help='Path to the multiqc folder')
    parser.add_argument('-o','--output', type=str, help='Path to the output file')
    args = parser.parse_args()
    main(f'{args.multiqc_folder}/multiqc_general_stats.txt', f'{args.multiqc_folder}/multiqc_samtools_stats.txt',
         f'{args.multiqc_folder}/multiqc_picard_dups.txt', args.output)
