#!/usr/bin/env python3

from pysam import VariantFile
import argparse
import subprocess
from subprocess import Popen


class VcfWriteProc:
    def __init__(self, pref, header, as_bcf):
        fn = (pref + '.bcf') if as_bcf else (pref + '.vcf.gz')
        mode = '-Ob' if as_bcf else '-Oz'
        self.Popen = Popen(['bcftools', 'view', mode, '--no-version', '-o', fn],
                           stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        self.variantFile = VariantFile(self.Popen.stdin, 'wu', header=header)

    def close(self):
        self.variantFile.close()
        self.Popen.communicate()
        return self.Popen.poll()

    def write(self, rec):
        self.variantFile.write(rec)


def main(pref, chunk_size, as_bcf):
    n = 0
    i = 1
    vf_in = VariantFile('-')
    vf_out = VcfWriteProc("{}-{}".format(pref, i), vf_in.header, as_bcf)
    for rec in vf_in:
        if n == chunk_size:
            vf_out.close()
            n = 0
            i += 1
            vf_out = VcfWriteProc("{}-{}".format(pref, i), vf_in.header, as_bcf)
        vf_out.write(rec)
        n += 1
    vf_out.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pref', default='output', help='output variant files prefix')
    parser.add_argument('--chunk-size', required=True, type=int, help='ideal (max) size of chunks')
    parser.add_argument('--bcf', action='store_true', help='write bcf instead of vcf.gz')
    args = parser.parse_args()
    main(args.pref, args.chunk_size, args.bcf)
