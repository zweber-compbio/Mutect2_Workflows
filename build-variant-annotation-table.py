# File Name:
# Created By:
# Created On:
# Purpose:
# License:

import argparse
import pandas as pd

# stream vep annotation file and split into data dict, header, and data lines
# write out to file a data dictionary, and return a pandas dataframe
def parse_vep_input(file_name, output_prefix):
    output_filename = output_prefix + "-datadict.txt"
    datadictobj = open(output_filename, 'w')
    vepobj = open(file_name, 'r')
    tabvep = []
    for line in vepobj:
        if line[0:2] == "##":
            towrite = line.lstrip("#")
            datadictobj.write(towrite)
        elif line[0] == "#":
            header = line.lstrip("#").split('\t')
        else:
            tabvep.append(line.split('\t'))
    return pd.DataFrame(tabvep, columns=header)

# stream sample names file, expanding into proper set of column names
def parse_sample_names(sample_names, patterns):
    s_names, output_colnames = [], []
    with open(sample_names, 'r') as fobj:
        for line in fobj: s_names.append(line.rstrip())
    for s in s_names:
        output_colnames.extend(["{}.{}".format(s,p) for p in patterns])
    return output_colnames, s_names

# stream vcf info file, converting to properly labeled pandas dataframe
def parse_var_input(file_name, sample_names, colnames):
    return pd.read_csv(file_name, sep="\t", names=colnames)

# compute allele_frequency related statistics and resort columns
def adjust_allele_frequency_columns(jointdf, sample_names):
    cols_to_add = {}
    for s in sample_names:
        adepthcolnames = s + ".AD"
        refacolnames = s + ".REF"
        altacolnames = s + ".ALT"
        dectcolnames = s + ".DETECTED"

        splitalleledepths = jointdf[adepthcolnames].str.split(',', expand=True)
        detected = [int(a) > 0 for a in splitalleledepths[1]]
        cols_to_add[refacolnames] = splitalleledepths[0]
        cols_to_add[altacolnames] = splitalleledepths[1]
        cols_to_add[dectcolnames] = detected

        del jointdf[adepthcolnames]

    newcolumns = pd.DataFrame(cols_to_add)
    return pd.concat([jointdf, newcolumns], axis=1)


if __name__ == "__main__":
    # setup CLI argument intake
    parser = argparse.ArgumentParser()
    parser.add_argument("--tab", type=str)
    parser.add_argument("--vcf-info", type=str)
    parser.add_argument("--sample-names", type=str)
    parser.add_argument("--output-filename-prefix",type=str)
    args = parser.parse_args()

    # parse the VEP and VCF info table to generate pandas dataframes
    sample_cnames, sample_names = parse_sample_names(args.sample_names, ["DP","AD"])
    vcfcnames = ["FILTER","TLOD","NLOD","NALOD"]; vcfcnames.extend(sample_cnames)
    vepdata = parse_vep_input(args.tab, args.output_filename_prefix)
    vardata = parse_var_input(args.vcf_info, args.sample_names, vcfcnames)
    joined_df = pd.concat([vepdata, vardata], axis=1)

    # create and sort columns related to allele frequency statistics
    alleleadj_df = adjust_allele_frequency_columns(joined_df, sample_names)

    # write out post-processed file
    writeoutname = "{}-annotation-table.txt".format(args.output_filename_prefix)
    alleleadj_df.to_csv(writeoutname, sep="\t", index=False)
