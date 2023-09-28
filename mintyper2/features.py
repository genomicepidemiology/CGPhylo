import os
import sys

def produce_features(args):
    os.mkdir('output')
    os.system('kmc -k 5 {} output/5mers .'.format(args.input))
    return 'test'