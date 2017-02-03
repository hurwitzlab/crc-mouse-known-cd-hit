"""
  1. copy dna-CDS to work directory
  2. translate dna
  3. plot sequence length histogram for dna, untranslated dna, protein
  4. run cd-hit on protein
"""
import argparse
import os
import io
import itertools
import stat
import subprocess


def pipeline():
    args = get_args()

    qsub_params = [
        '-l', 'jobtype=serial',
        '-l', 'select=1:ncpus=1:mem=1gb',
        '-l', 'place=pack:shared',
        '-l', 'walltime=01:00:00',
        '-M', 'jklynch@email.arizona.edu',
        '-m', 'bea',
        '-q', 'standard',
        '-W', 'group_list=bhurwitz'
    ]

    work_script_dir = os.path.join(args.work_dir, 'script')
    if not os.path.exists(work_script_dir):
        os.makedirs(work_script_dir)

    _, dna_cds_known_file_name = os.path.split(args.orig_dna_cds_known_path)
    work_dna_cds_known_path = os.path.join(args.work_dir, dna_cds_known_file_name)

    ###########################################################################
    # copy Scott's original data file
    run_script(
        write_script(script_path=os.path.join(args.work_dir, 'script', 'copy.sh'), script_text="""\
            #!/bin/bash
            pwd
            mkdir -p {work_dir}
            cp {orig_dna_cds_known_path} {work_dna_cds_known_path}
            """.format(**vars(args), work_dna_cds_known_path=work_dna_cds_known_path)
        )
    )
    ###########################################################################

    ###########################################################################
    # translate dna to protein
    qsub_script(
        write_script(script_path=os.path.join(args.work_dir, 'script', 'translate.sh'), script_text= """\
            #!/bin/bash
            source activate mouse
            python {scripts_dir}/translate-microbial-dna-CDS.py \\
                -i {work_dna_cds_known_path} \\
                -o {work_dir}/crc-mouse-protein-from-known-only.fa \\
                -u {work_dir}/crc-mouse-untranslated-microbial-dna-CDS-known.fa \\
                -l {translation_limit}
            """.format(**vars(args), work_dna_cds_known_path=work_dna_cds_known_path)
        ),
        job_name='translate',
        qsub_params=qsub_params
    )
    ###########################################################################

    ###########################################################################
    # cluster proteins with CD-HIT
    qsub_script(
        write_script(script_path=os.path.join(args.work_dir, 'script', 'cluster_proteins.sh'), script_text="""\
            #!/bin/bash
            module load cdhit
            cdhit \\
                -i {work_dir}/crc-mouse-protein-from-known-only.fa \\
                -o {work_dir}/crc-mouse-cd-hit-c90-n5-protein-known.db \\
                -c 0.9 -n 5 -M 1000 -d 0 -T 1
            """.format(**vars(args))
        ),
        job_name='cdhit',
        qsub_params=qsub_params
    )
    ###########################################################################

def get_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-w', '--work-dir', default=os.path.join(os.getcwd(), 'work'))
    arg_parser.add_argument('-s', '--scripts-dir', default=os.path.join(os.getcwd(), 'scripts'))
    arg_parser.add_argument(
        '-i', '--orig-dna-cds-known-path',
        default='/rsgrps/bhurwitz/scottdaniel/extract-fasta/data/dna-of-CDS-from-known-only.fa')
    arg_parser.add_argument('-l', '--translation-limit', type=int, default=-1)

    return arg_parser.parse_args()

def write_script(script_text, script_path):
    with open(script_path, 'wt') as script_file:
        for line in io.StringIO(script_text):
            script_file.write(line.lstrip())
    os.chmod(script_path, stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR)
    return script_path

def run_script(script_path):
    print('running "{}"'.format(script_path))
    p = subprocess.run(
        [script_path],
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)
    print('stderr:\n{}'.format(p.stderr.decode('utf-8')))
    print('stdout:\n{}'.format(p.stdout.decode('utf-8')))

def qsub_script(script_path, job_name, qsub_params):
    print('qsub "{}"'.format(script_path))
    subprocess_cmd_list = ['qsub', '-N', job_name] + qsub_params + [script_path]
    print(subprocess_cmd_list)
    p = subprocess.run(
        subprocess_cmd_list,
        shell=False,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)
    print('stderr:\n{}'.format(p.stderr.decode('utf-8')))
    print('stdout:\n{}'.format(p.stdout.decode('utf-8')))


if __name__ == '__main__':
    pipeline()