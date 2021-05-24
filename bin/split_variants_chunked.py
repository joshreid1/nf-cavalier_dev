#!/usr/bin/env python3

from pysam import VariantFile
import argparse
import subprocess
from subprocess import Popen
from queue import Queue
from threading import Thread
import math


class VcfReadThread(Thread):
    def __init__(self, fn, q):
        Thread.__init__(self, daemon=True)
        self.Popen = Popen(['bcftools', 'view', '-Ou', '--no-version', fn],
                           stdin=subprocess.DEVNULL, stdout=subprocess.PIPE)
        self.variantFile = VariantFile(self.Popen.stdout)
        self.q = q

    def run(self):
        for j, rec in enumerate(self.variantFile):
            self.q.put((j, rec))
        self.q.put((-1, -1))
        self.variantFile.close()
        self.Popen.communicate()
        self.Popen.poll()


class VcfWriteProc:
    def __init__(self, fn, header):
        self.Popen = Popen(['bcftools', 'view', '-Oz', '--no-version', '-o', fn],
                           stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        self.variantFile = VariantFile(self.Popen.stdin, 'wu', header=header)

    def close(self):
        self.variantFile.close()
        self.Popen.communicate()
        return self.Popen.poll()

    def write(self, rec):
        self.variantFile.write(rec)


def nrec(fn):
    popen = Popen(['bcftools', 'index', '--nrecords', fn],
                  stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    return int(popen.stdout.readline())


def main(input, out_pref, n_chunk, chunk_size):
    if n_chunk != -1:
        chunk_size = int(math.ceil(nrec(input) / n_chunk))
    q = Queue(10)
    reader = VcfReadThread(input, q)
    reader.start()
    n = 0
    i = 1
    writer = VcfWriteProc("{}-{}.vcf.gz".format(out_pref, i), reader.variantFile.header)
    while True:
        j, rec = q.get()
        if j == -1:
            break
        if n == chunk_size:
            writer.close()
            n = 0
            i += 1
            writer = VcfWriteProc("{}-{}.vcf.gz".format(out_pref, i), reader.variantFile.header)
        writer.write(rec)
        q.task_done()
        n += 1
    writer.close()
    print('Done. Wrote {} records over {} files.'.format(n, i))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input variant file (.vcf, .vcf.gz, .bcf)')
    parser.add_argument('--out', default='output', help='output variant files prefix')
    parser.add_argument('--n-chunk', default=-1, type=int, help='number of chunks')
    parser.add_argument('--chunk-size', default=-1, type=int, help='size of chunks')
    args = parser.parse_args()
    if args.n_chunk == -1 and args.chunk_size == -1:
        print('Error: must specify one of --n-chunk or --chunk-size')
        exit(1)
    main(args.input, args.out, args.n_chunk, args.chunk_size)
