import os
import sys

def produce_features(args):
    os.mkdir('output')
    os.system('kmc -k1 -cs1000000000000 {} output/1mers . > output/1mers_stats.txt'.format(args.input[0]))
    os.system('kmc_dump output/1mers output/1mers.txt')
    os.system('kmc -k5 -cs1000000000 {} output/5mers . > output/5mers_stats.txt'.format(args.input[0]))
    os.system('kmc_dump output/5mers output/5mers.txt')
    os.system('kmc -k13 -cs1000000000 {} output/13mers . > output/13mers_stats.txt'.format(args.input[0]))
    os.system('kmc_dump output/13mers output/13mers.txt')
    os.system('kmc -k21 -cs1000000000 {} output/21mers . > output/21mers_stats.txt'.format(args.input[0]))
    os.system('kmc_dump output/21mers output/21mers.txt')
    return 'test'