import os
import sys
import math
import shutil
import argparse
import itertools
import yaml
import pandas as pd
import numpy as np
import mosaic.io as mio
from pathlib import Path

def get_filtered_mutations(manually_annotated_snvs):
    if os.stat(manually_annotated_snvs).st_size != 0:
        snvs = pd.read_csv(manually_annotated_snvs, sep="\t", comment="#")
        snvs.replace(np.nan, "", inplace=True)
        # remove the ones labeled as artifact
        snvs = snvs.loc[~snvs["annotation"].str.lower().str.contains("artifact")]
        germline_mutations = snvs[(snvs["annotation"] == "germline_HET")][
            "condensed_format"
        ].to_list()
        total_mutations = snvs[
            (snvs["annotation"] != "likely_artifact")
            & (snvs["annotation"] != "germline_HOM")
        ]["condensed_format"].to_list()
        somatic_mutations = snvs[snvs["annotation"] == "bulk_somatic"][
            "condensed_format"
        ].to_list()
    else:
        raise ValueError("No manually annotated SNVs found")

    germline_mutations = list(set(germline_mutations))
    somatic_mutations = list(set(somatic_mutations))
    total_mutations = list(set(total_mutations))

    print(
        f"""
    germline: {len(germline_mutations)}
    somatic: {len(somatic_mutations)}
    total: {len(total_mutations)}"""
    )

    return total_mutations, germline_mutations, somatic_mutations

def generate_condor_input(sample_name, cn_assignment_df, args, bin_thres=0.5):
    merged_cn_assignment_df = cn_assignment_df.copy()
    print(merged_cn_assignment_df)
    if "Unnamed: 0" in merged_cn_assignment_df.columns:
        merged_cn_assignment_df = merged_cn_assignment_df.rename(
            columns={"Unnamed: 0": "sample"}
        )
    merged_cn_assignment_df["cell_barcode"] = (
        merged_cn_assignment_df["sample"].astype(str)
        + ":"
        + merged_cn_assignment_df["cell_barcode"].astype(str)
    )
    merged_cn_assignment_df.drop(["sample"], axis=1, inplace=True)

    df_total_samples = []
    df_alt_samples = []
    common_mutations = []
    germline_mutations = []
    somatic_mutations = []
    # embed()
    for file in (Path(args.l) / sample_name).iterdir():
        file = file.name
        if str(file).startswith(sample_name):
            if str(file).endswith(".h5"):
                print(f"Found H5: {file}")
                hdf5_f = Path(args.l) / sample_name / file
                sample_obj = mio.load(hdf5_f)
                sample_obj.dna.genotype_variants(
                    min_dp=8,
                    min_alt_read=3,
                    assign_low_conf_genotype=True,
                )
                df_alt_snv = sample_obj.dna.get_attribute(
                    "alt_read_count", constraint="row"
                )
                df_total_snv = sample_obj.dna.get_attribute("DP", constraint="row")

                print(df_alt_snv.shape)
                cmuts, gmuts, smuts = get_filtered_mutations(args.snvs)
                common_mutations = cmuts
                germline_mutations = gmuts
                somatic_mutations = smuts

                df_alt_snv.reset_index(inplace=True)
                df_alt_snv = df_alt_snv.rename(columns={"index": "cell_barcode"})

                df_alt_samples.append(df_alt_snv)

                df_total_snv.reset_index(inplace=True)
                df_total_snv = df_total_snv.rename(columns={"index": "cell_barcode"})

                df_total_samples.append(df_total_snv)

    def mut_replace(x):
        x = x.replace(":", "_").replace("/", "_").split("_")
        x[2], x[3] = x[3], x[2]
        return "_".join(x)

    # common_mutations = list(map(mut_replace, common_mutations))
    # germline_mutations_list = list(map(mut_replace, germline_mutations))
    germline_mutations_list = list(germline_mutations)
    # somatic_mutations_list = list(map(mut_replace, somatic_mutations))
    somatic_mutations_list = list(somatic_mutations)
    # print(len(common_mutations))

    df_total = pd.concat(df_total_samples, ignore_index=True)
    df_total = df_total.rename(columns={"DP": "cell_barcode"})
    # df_total = df_total.rename(
    #     columns={
    #         c: mut_replace(c) for c in df_total.columns if c not in ["cell_barcode"]
    #     }
    # )

    df_total = df_total[["cell_barcode"] + common_mutations]
    df_total = df_total.set_index("cell_barcode")
    df_total = df_total.fillna(0)

    df_alt = pd.concat(df_alt_samples)
    df_alt = df_alt.rename(columns={"alt_read_count": "cell_barcode"})
    # df_alt = df_alt.rename(
    #     columns={c: mut_replace(c) for c in df_alt.columns if c not in ["cell_barcode"]}
    # )

    df_alt = df_alt[["cell_barcode"] + common_mutations]
    df_alt = df_alt.set_index("cell_barcode")
    df_alt = df_alt.fillna(0)

    print(df_total.shape)
    print(df_alt.shape)

    df_character_mat = df_total.copy()
    df_character_mat = df_alt.divide(df_total)

    print(df_character_mat.index)
    print(merged_cn_assignment_df["cell_barcode"])

    def rename_barcode(s):
        return s.split(":")[1] + "-" + s.split(":")[0]

    merged_cn_assignment_df["cell_barcode"] = merged_cn_assignment_df[
        "cell_barcode"
    ].apply(rename_barcode)
    print(merged_cn_assignment_df["clone_id"].value_counts())
    print(
        set([c.split("-")[2] for c in df_character_mat.index.tolist()]),
        merged_cn_assignment_df.shape,
    )

    df_character_mat = pd.merge(
        df_character_mat,
        merged_cn_assignment_df,
        left_on=df_character_mat.index,
        right_on="cell_barcode",
        how="left",
    )
    df_character_mat = df_character_mat.set_index("cell_barcode")
    print(df_character_mat["clone_id"].value_counts())

    df_character_mat.rename(columns={"clone_id": "cluster_id"}, inplace=True)

    l_ids = list(df_character_mat["cluster_id"].unique())
    print(l_ids)
    l_ids.sort()
    l_ids_dict = {}
    index = 0
    for i in l_ids:
        l_ids_dict[i] = index
        index += 1

    df_character_mat["cluster_id"] = df_character_mat["cluster_id"].replace(l_ids_dict)

    df_character_vaf = df_character_mat.copy()
    df_character_mat[df_character_mat.columns[:-1]] = df_character_mat[
        df_character_mat.columns[:-1]
    ].applymap(lambda x: -1 if pd.isna(x) else 0 if x < bin_thres else 1)

    return (
        df_total,
        df_alt,
        df_character_mat,
        df_character_vaf,
        germline_mutations_list,
        somatic_mutations_list,
    )


def main(args):
    dataset = args.d
    cn_assignment_df = pd.read_csv(args.i)
    (
        df_total,
        df_alt,
        df_character_mat,
        df_character_vaf,
        germline_mutations,
        somatic_mutations,
    ) = generate_condor_input(dataset, cn_assignment_df, args)
    with open(args.g, "w") as fp:
        for item in germline_mutations:
            fp.write("%s\n" % item)

    with open(args.s, "w") as fp:
        for item in somatic_mutations:
            fp.write("%s\n" % item)

    df_character_vaf.to_csv(args.v, index_label="")
    df_character_mat.to_csv(args.m, index_label="")
    df_alt.to_csv(args.a, index_label="")
    df_total.to_csv(args.t, index_label="")

    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", type=str, help="dataset name")
    parser.add_argument("-l", type=str, help="input path for hdf5 and CRAVAT files")
    parser.add_argument(
        "-i", type=str, help="input path to refined cell assignments csv file"
    )
    # parser.add_argument(
    #     "-c", type=str, help="input path to snv selection config parameter file"
    # )
    parser.add_argument(
        "-snvs", type=str, default=None, help="input path to manually annotated SNVs"
    )
    parser.add_argument("-v", type=str, help="output path for ConDoR VAF matrix")
    parser.add_argument("-m", type=str, help="output path for ConDoR binary matrix")
    parser.add_argument(
        "-a", type=str, help="output path for ConDoR alternate readcount matrix"
    )
    parser.add_argument(
        "-t", type=str, help="output path for ConDoR total readcount matrix "
    )
    parser.add_argument(
        "-g", type=str, help="output path for ConDoR germline mutation list"
    )
    parser.add_argument("-s", type=str, help="output path for somatic mutation list")

    args = parser.parse_args(None if sys.argv[1:] else ["-h"])
    main(args)
