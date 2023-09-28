import os
import sys

def produce_features(args):
    os.mkdir('output')
    os.system('kmc -k5 -cs1000000000 {} output/5mers .'.format(args.input[0]))
    os.system('kmc_dump output/5mers output/5mers.txt')
    os.system('kmc -k13 -cs1000000000 {} output/13mers .'.format(args.input[0]))
    os.system('kmc_dump output/13mers output/13mers.txt')
    os.system('kmc -k21 -cs1000000000 {} output/21mers .'.format(args.input[0]))
    os.system('kmc_dump output/21mers output/21mers.txt')
    return 'test'