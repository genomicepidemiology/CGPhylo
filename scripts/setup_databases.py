#!/usr/bin/env python3

import argparse
import sys
import os
import logging
import subprocess

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#NOTE Make sure to remove Burkholderia mallei (RKI) and only use Burkholderia mallei (FLI). Make sure to rename the folder to Burkholderia_mallei_cgMLST_alleles.

def main(args):
    """Main function"""
    check_kma_installed()
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    decompress_tar_archive(args.zip, args.output)
    allele_folders = os.listdir(args.output)
    for allele_folder in allele_folders:
        if not allele_folder.startswith('.'):
            path = args.output + '/' + allele_folder
            rename_alleles(path)
            concat_renamed_fasta_files(path)
            os.system('kma index -i {}/{}_complete.fasta -o {}/{}_complete'.format(path, allele_folder, path, allele_folder))
            os.system('rm {}/*.fasta'.format(path))

def check_kma_installed():
    try:
        # Attempt to execute the KMA command and get its version
        result = subprocess.run(["kma", "-v"], capture_output=True, text=True, check=True)
        if result.returncode == 0:
            return "KMA is installed and available."
        else:
            return "KMA is installed but there was an error executing it."
    except FileNotFoundError:
        return "KMA is not installed or not available in PATH."


def clean(path):
    # Remove hidden files
    files = os.listdir(path)
    for file in files:
        sub_files = os.listdir(path + '/' + file)
        for sub_file in sub_files:
            if sub_file.startswith('._'):
                os.remove(path + '/' + file + '/' + sub_file)
            if sub_file.startswith('rename'):
                os.remove(path + '/' + file + '/' + sub_file)
            if 'complete' in sub_file:
                os.remove(path + '/' + file + '/' + sub_file)


def concat_renamed_fasta_files(path):
    species_name = path.split('/')[-1]
    os.system('cat {}/renamed_* > {}/{}_complete.fasta'.format(path, path, species_name))

def rename_alleles(path):
    files = os.listdir(path)

    for file in files:
        file_name = file.split('.')[0]
        if not file.startswith('.'):
            with open(path + '/renamed_' + file, 'w') as outfile:
                with open(path + '/' + file, 'r') as f:
                    print ('renamed: ' + path + '/' + file)
                    for line in f:
                        if line.startswith('>'):
                            number = line.strip()[1:]
                            print('>{}_{}'.format(file_name, number), file=outfile)
                        else:
                            print(line.strip(), file=outfile)


def decompress_tar_archive(path, output):
    """Decompress a tar archive"""
    os.system('tar -xvf {} -C {}'.format(path, output))


if __name__ == '__main__':
    # initialize the options parser
    parser = argparse.ArgumentParser('build_cg_db', add_help=False)

    parser.add_argument('--zip', action="store", type=str, dest='zip',
                        help='Zip archive of cgmlst db which have been downloaded.')
    parser.add_argument('--output', action="store", type=str, default='output', dest="output",
                        help="Output directory")
    parser.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    main(args)