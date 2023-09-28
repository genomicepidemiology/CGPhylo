import os
import sys

def produce_features(args):
    os.system('kmc -k 5 {} 5mers .'.format(args.input))
    return 'test'