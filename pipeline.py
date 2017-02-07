"""
  1. copy dna-CDS to work directory
  2. translate dna
  3. plot sequence length histogram for dna, untranslated dna, protein
  4. run cd-hit on protein
"""
import argparse
import os
import io
import stat
import subprocess


def pipeline():
    args = get_args()

    qsub_params = [
        '-l', 'place=pack:shared',
        '-M', 'jklynch@email.arizona.edu',
        '-m', 'bea',
        '-q', 'standard',
        '-W', 'group_list=bhurwitz'
    ]

    work_scripts_dir = os.path.join(args.work_dir, 'scripts')
    if not os.path.exists(work_scripts_dir):
        os.makedirs(work_scripts_dir)

    _, dna_cds_known_file_name = os.path.split(args.orig_dna_cds_known_path)
    work_dna_cds_known_path = os.path.join(args.work_dir, dna_cds_known_file_name)

    ###########################################################################
    # copy Scott's original data file
    if os.path.exists(work_dna_cds_known_path):
        print('"{}" already exists'.format(work_dna_cds_known_path))
    else:
        run_script(
            write_script(script_path=os.path.join(work_scripts_dir, 'copy.sh'), script_text="""\
                #!/bin/bash
                pwd
                mkdir -p {work_dir}
                cp {orig_dna_cds_known_path} {work_dna_cds_known_path}
                """.format(**vars(args), work_dna_cds_known_path=work_dna_cds_known_path),
                 job_name='crc-mouse-copy',
                 select=1,
                 ncpus=1,
                 mem='6gb',
                 pcmem='6gb',
                 place='free:shared',
                 walltime='00:10:00',
                 cputime='00:10:00',
                 stderr_fp='crc-mouse-copy.stderr',
                 stdout_fp='crc-mouse-copy.stdout',
                 qsub_params=qsub_params
            )
        )
    ###########################################################################

    ###########################################################################
    # translate dna to protein
    qsub_script(
        write_script(script_path=os.path.join(work_scripts_dir, 'translate.sh'), script_text= """\
            #!/bin/bash
            source activate mouse
            python {scripts_dir}/translate-microbial-dna-CDS.py \\
                -i {work_dna_cds_known_path} \\
                -o {work_dir}/crc-mouse-protein-from-known-only.fa \\
                -u {work_dir}/crc-mouse-untranslated-microbial-dna-CDS-known.fa \\
                -l {translation_limit}
            """.format(**vars(args), work_dna_cds_known_path=work_dna_cds_known_path),
            job_name='crc-mouse-translate',
            select=1,
            ncpus=1,
            mem='6gb',
            pcmem='6gb',
            place='free:shared',
            walltime='00:30:00',
            cputime='00:30:00',
            stderr_fp='mouse_translate.stderr',
            stdout_fp='mouse_translate.stdout',
            qsub_params=qsub_params
        )
    )
    ###########################################################################

    ###########################################################################
    # cluster proteins with CD-HIT
    qsub_script(
        write_script(script_path=os.path.join(work_scripts_dir, 'cluster_proteins.sh'), script_text="""\
            #!/bin/bash
            module load cd-hit
            cdhit \\
                -i {work_dir}/crc-mouse-protein-from-known-only.fa \\
                -o {work_dir}/crc-mouse-cd-hit-c90-n5-protein-known.db \\
                -c 0.9 -n 5 -M 168000 -d 0 -T 28
            """.format(**vars(args)),
            job_name='crc-mouse-cdhit',
            select=1,
            ncpus=28,
            mem='168gb',
            pcmem='6gb',
            place='pack:free',
            walltime='01:00:00',
            cputime='28:00:00',
            stderr_fp='mouse_cluster.stderr',
            stdout_fp='mouse_cluster.stdout',
            qsub_params=qsub_params
        )
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


def write_script(script_path, script_text, **kwargs):
    print(kwargs)
    with open(script_path, 'wt') as script_file:
        script_text_buffer = io.StringIO(script_text)
        script_file.write(script_text_buffer.readline().lstrip())
        script_file.write("""\
#PBS -N {job_name}
#PBS -q standard
#PBS -W group_list=bhurwitz
#PBS -l select={select}:ncpus={ncpus}:mem={mem}:pcmem={pcmem}
#PBS -l place={place}
#PBS -l cputime={cputime}
#PBS -l walltime={walltime}
#PBS -m bea
#PBS -M jklynch@email.arizona.edu
#PBS -e {stderr_fp}
#PBS -o {stdout_fp}
""".format(**kwargs))
        for line in script_text_buffer:
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


def qsub_script(script_path):
    print('qsub "{}"'.format(script_path))
    subprocess_cmd_list = ['qsub', script_path]
    print(subprocess_cmd_list)
    try:
        p = subprocess.run(
            subprocess_cmd_list,
            shell=False,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE)
        print('stderr:\n{}'.format(p.stderr.decode('utf-8')))
        print('stdout:\n{}'.format(p.stdout.decode('utf-8')))
    except FileNotFoundError as e:
        # this usually means I am testing on my laptop
        print('no qsub executable')
        run_script(script_path=script_path)


if __name__ == '__main__':
    pipeline()