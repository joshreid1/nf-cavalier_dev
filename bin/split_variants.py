#!/usr/bin/env python3

from pysam import VariantFile
import argparse
import subprocess
from subprocess import Popen, check_call


class VcfWriteProc:
    def __init__(self, pref, header, as_bcf):
        self.as_bcf = as_bcf
        self.fn = (pref + '.bcf') if as_bcf else (pref + '.vcf.gz')
        mode = '-Ob' if as_bcf else '-Oz'
        self.Popen = Popen(['bcftools', 'view', mode, '--no-version', '-o', self.fn],
                           stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        self.variantFile = VariantFile(self.Popen.stdin, 'wu', header=header)

    def index(self):
        args = ['bcftools', 'index', self.fn] if self.as_bcf else ['bcftools', 'index', '-t', self.fn]
        return check_call(args, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)

    def close(self, index=False):
        self.variantFile.close()
        self.Popen.communicate()
        if index:
            self.Popen.poll()
            return self.index()
        return self.Popen.poll()

    def write(self, rec):
        self.variantFile.write(rec)


def main(pref, chunk_size, as_bcf, do_index):
    n = 0
    i = 1
    ix = '01'
    vf_in = VariantFile('-', 'rb')
    vf_out = VcfWriteProc("{}-{}".format(pref, ix), vf_in.header, as_bcf)
    for rec in vf_in:
        if n == chunk_size:
            vf_out.close(index=do_index)
            n = 0
            i += 1
            ix = ('0' + str(i)) if i < 10 else str(i)
            vf_out = VcfWriteProc("{}-{}".format(pref, ix), vf_in.header, as_bcf)
        vf_out.write(rec)
        n += 1
    vf_out.close(index=do_index)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pref', default='output', help='output variant files prefix')
    parser.add_argument('--chunk-size', required=True, type=int, help='ideal (max) size of chunks')
    parser.add_argument('--bcf', action='store_true', help='write bcf instead of vcf.gz')
    parser.add_argument('--index', action='store_true', help='index output files')
    args = parser.parse_args()
    main(args.pref, args.chunk_size, args.bcf, args.index)
