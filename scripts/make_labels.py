#!/usr/bin/env python
from sys import argv
import argparse

def parse_args(argv):
    desc = 'Translates assembly accession IDs to representative OTU IDs'
    copyright = 'Copyright (c) David Lund 2024.'
    parser = argparse.ArgumentParser(description=desc+'. '+copyright)
    parser.add_argument('--list', '-l', required=True,
                        help='List of the desired length.')
    parser.add_argument('--word', '-w', required=True,
                        help='Word to repeat')
    parser.add_argument('--header', required=True,
                        help='Header of output file')
    parser.add_argument('--output', '-o', required=True,
                        help='Name of output file.')
    
    arguments = parser.parse_args()
    
    return arguments

def make_labels(list, name):
    labels = []
    for i in range(len(list)-1):
        labels.append(name + "\n")

    return labels

def main():
    arguments = parse_args(argv)

    with open(arguments.list, 'r') as file:
        list = file.readlines()

    labels = make_labels(list, arguments.word)

    stringToExport = str(arguments.header) + "\n"

    for row in labels:
        stringToExport += str(row)

    with open(arguments.output, 'w') as file:
        file.write(stringToExport)

if __name__ == '__main__':
    main()
