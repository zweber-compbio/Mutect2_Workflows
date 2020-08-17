# File Name: prepare_multivcf.py
# Created By: ZWeber
# Created On: 2020-05-30
# Purpose: prepares annotations from a Multisample VCF created by Mutect2
# LICENSE: software not developed by Z.Weber retains original licenses.
#   while software developed by Z.Weber can be used under the conditions
#   outlined by the GNU-GPLv3.0

# external libraries
import argparse
import gzip

# Main Execution
if __name__ == "__main__":
    # command line  argument capture
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", type=str)
    parser.add_argument("-o","--output", type=str)
    args = parser.parse_args()

    # if gzipped, unzip
    isgzipped = (args.input[-3:] == ".gz")

    if (isgzipped):
        invcf = gzip.open(args.input, "r")
        outvcf = open(args.output, "w")

        for line in invcf:
            outvcf.write(line.decode())
        invcf.close()
        outvcf.close()
    else:
        invcf = open(args.input, "r")
        outvcf = open(args.output, "w")

        for line in invcf:
            outvcf.write(line)
        invcf.close()
        outvcf.close()
