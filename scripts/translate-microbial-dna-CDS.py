import argparse
import io
import itertools

from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord


def translate(dna_cds_input_fp, protein_output_fp, untranslated_dna_output_fp, sequence_limit):
    # translate one CDS
    with open(dna_cds_input_fp, 'rt') as dna_cds_input_file:
        first_cds = [next(SeqIO.parse(dna_cds_input_file, 'fasta'))]
        print(first_cds)
        first_protein = next(make_bacterial_protein_record(first_cds, io.StringIO()))
        print(first_protein)

    with open(dna_cds_input_fp, 'rt') as dna_cds_input_file,\
            open(protein_output_fp, 'wt') as protein_output_file,\
            open(untranslated_dna_output_fp, 'wt') as untranslated_dna_file:

        if sequence_limit < 0:
            sequence_limit = None

        sequence_count = SeqIO.write(
            make_bacterial_protein_record(
                nucleotide_record_input=itertools.islice(SeqIO.parse(dna_cds_input_file, 'fasta'), sequence_limit),
                untranslated_dna_output=untranslated_dna_file
            ),
            protein_output_file,
            'fasta')
        print('wrote {} sequences to "{}"'.format(sequence_count, protein_output_fp))


def make_bacterial_protein_record(nucleotide_record_input, untranslated_dna_output):
    """Returns a new SeqRecord with the translated sequence using the Bacterial table."""
    failed_translation_count = 0

    for n, nuc_record in enumerate(nucleotide_record_input):
        if (n+1) % 10000 == 0:
            print('.')

        try:
            yield SeqRecord(
                seq=nuc_record.seq.translate(table='Bacterial', cds=True),
                id='protein_' + nuc_record.id,
                description='translation of known Bacterial CDS')
        except TranslationError as t:
            failed_translation_count += 1
            print('{}: failed to translate nucleotide record\n{}'.format(failed_translation_count, nuc_record))
            untranslated_dna_output.write(nuc_record.format('fasta'))

    print('{} un-translated sequences'.format(failed_translation_count))

def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i', '--dna-cds-input-fp')
    arg_parser.add_argument('-o', '--protein-output-fp')
    arg_parser.add_argument('-u', '--untranslated-dna-output-fp')
    arg_parser.add_argument('-l', '--sequence-limit', type=int, default=-1)

    return arg_parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    translate(**args.__dict__)
