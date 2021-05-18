#!/usr/bin/env python3

from pysam import VariantFile
import argparse


def main(in_variants, out_variants,  i, n, mode):
    variants = VariantFile(in_variants)

    if mode == 'v':
        mode = ''
    out = VariantFile(out_variants, 'w' + mode, header=variants.header)

    for j, rec in enumerate(variants.fetch()):
        if (j % n) + 1 == i:
            out.write(rec)
    out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('in_variants', help='input variant file (.vcf, .vcf.gz, .bcf)')
    parser.add_argument('i', type=int, help='split index')
    parser.add_argument('n', type=int, help='number of splits')
    parser.add_argument('--out', default='-', help='output variant file, default to stdout')
    parser.add_argument('--mode', default='v', help='output variant type (v - VCF, b - BCF, u - uncompressed BCF)')
    args = parser.parse_args()
    if args.mode not in ['v', 'b', 'u']:
        print('Error: Mode must be one of "v", "b" or "u')
        exit(1)
    main(args.in_variants, args.out, args.i, args.n, args.mode)
