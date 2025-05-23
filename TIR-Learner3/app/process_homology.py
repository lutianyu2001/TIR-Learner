#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *

BLAST_HEADER_FULL_COVERAGE = ("qacc", "sacc", "length", "pident", "gaps", "mismatch",
                              "qstart", "qend", "sstart", "send", "evalue", "qcovhsp")

BLAST_HEADER_EIGHTY_SIMILARITY = ("qseqid", "sseqid", "length", "pident", "gaps", "mismatch",
                                  "qstart", "qend", "sstart", "send", "evalue", "qcovhsp")

BLAST_TYPE = {"length": int, "gaps": int, "mismatch": int,
              "qstart": int, "qend": int, "sstart": int, "send": int}


def _process_homology_full_coverage(genome_name: str, species: str, TIR_type: str) -> Optional[pd.DataFrame]:
    blast = f"{genome_name}{SPLITER}blast{SPLITER}{species}_{TIR_type}_RefLib"
    df = None
    if os.path.exists(blast) and os.path.getsize(blast) != 0:
        df = pd.read_csv(blast, sep='\t', header=None, names=BLAST_HEADER_FULL_COVERAGE, dtype=BLAST_TYPE, engine='c',
                         memory_map=True)
        df = df.loc[(df["qcovhsp"] == 100) & (df["pident"] >= 80)].reset_index(drop=True)
        df = df.sort_values(["sacc", "sstart", "send", "qcovhsp", "pident"],
                            ascending=[True, True, True, True, True], ignore_index=True)
        df = df.drop_duplicates(["sacc", "sstart", "send"], keep="last", ignore_index=True)
        df.insert(0, "TIR_type", TIR_type)
    return df


def _process_homology_eighty_similarity(file_name: str, species: str, TIR_type: str, flag_verbose: bool = True) -> Optional[pd.DataFrame]:
    blast = f"{file_name}{SPLITER}blast{SPLITER}{species}_{TIR_type}_RefLib"
    df = None
    if os.path.exists(blast) and os.path.getsize(blast) != 0:
        df = pd.read_csv(blast, sep='\t', header=None, names=BLAST_HEADER_EIGHTY_SIMILARITY, dtype=BLAST_TYPE,
                         engine='c', memory_map=True)
        df = df.loc[(df["qcovhsp"] >= 80) & (df["pident"] >= 80)].reset_index(drop=True)

        df["sseqid"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: x["qseqid"].split(":")[0], axis=1)
        df["sstart"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: int(x["qseqid"].split(":")[1]), axis=1)
        df["send"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: int(x["qseqid"].split(":")[2]), axis=1)

        df = df.sort_values(["sseqid", "sstart", "send", "qcovhsp", "pident"],
                            ascending=[True, True, True, True, True], ignore_index=True)
        df = df.drop_duplicates(["sseqid", "sstart", "send"], keep="last", ignore_index=True)
        df.insert(0, "TIR_type", TIR_type)
    return df


def _process_result(df_list: List[pd.DataFrame], species: str) -> pd.DataFrame:
    try:
        df = pd.concat(df_list, ignore_index=True).iloc[:, [0, 1, 2, 9, 10]].copy()
    except ValueError:
        # print(f"""
        # [ERROR] No sequence is found similar to the TIR database of {species}!
        # You may have specified the wrong species. Please double-check or set species=others and rerun TIR-Learner.
        # """)
        # sys.exit(-1)
        raise SystemExit(f"""
        [ERROR] No sequence is found similar to the TIR database of {species}!
        You may have specified the wrong species. Please double-check or set species=others and rerun TIR-Learner. 
        """)
    df = df.set_axis(["TIR_type", "id", "seqid", "sstart", "send"], axis=1)
    return df


def select_full_coverage(TIRLearner_instance) -> pd.DataFrame:
    print("Module 1, Step 2: Select 100% coverage entries from blast results")
    mp_args_list = [(TIRLearner_instance.genome_name, TIRLearner_instance.species, TIR_type)
                    for TIR_type in TIR_SUPERFAMILIES]
    with mp.Pool(int(TIRLearner_instance.processors)) as pool:
        df_list = pool.starmap(_process_homology_full_coverage, mp_args_list)
    # subprocess.Popen(["rm", "-f", f"*{spliter}blast{spliter}*"])  # remove blast files
    if not TIRLearner_instance.flag_debug:
        subprocess.Popen(["find", ".", "-name", f"*{SPLITER}blast{SPLITER}*", "-delete"])
    # subprocess.Popen(f"rm -f *{spliter}blast{spliter}*", shell=True)
    return _process_result(df_list, TIRLearner_instance.species)


def select_eighty_similarity(TIRLearner_instance) -> pd.DataFrame:
    print("Module 2, Step 7: Select 80% similar entries from blast results")
    mp_args_list = [(TIRLearner_instance.processed_de_novo_result_file_name,
                     TIRLearner_instance.species,
                     TIR_type,
                     TIRLearner_instance.flag_verbose)
                    for TIR_type in TIR_SUPERFAMILIES]
    with mp.Pool(int(TIRLearner_instance.processors)) as pool:
        df_list = pool.starmap(_process_homology_eighty_similarity, mp_args_list)
    if not TIRLearner_instance.flag_debug:
        subprocess.Popen(["find", ".", "-name", f"*{SPLITER}blast{SPLITER}*", "-delete"])
    return _process_result(df_list, TIRLearner_instance.species)
