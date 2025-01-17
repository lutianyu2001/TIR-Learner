from const import *

from run_TIRvish import process_fasta, retrieve_split_sequence_offset


def GRF(genome_file, genome_name, processors, TIR_length, GRF_path):
    GRF_bin_path = os.path.join(GRF_path, "grf-main")
    GRF_result_dir_name = f"{genome_name}_GRFmite"
    grf = (f"\"{GRF_bin_path}\" -i \"{genome_file}\" -o \"{GRF_result_dir_name}\" -c 1 -t {processors} -p 20 "
           f"--min_space 10 --max_space {int(TIR_length)} --max_indel 0 --min_tr 10 "
           f"--min_spacer_len 10 --max_spacer_len {int(TIR_length)}")
    shell_filter = r" | grep -vE 'start:|end:|print:|calculate|^$'"
    subprocess.Popen(grf + shell_filter, shell=True).wait()


def GRF_mp(genome_file_path, genome_name, processors, TIR_length, GRF_path):
    GRF_working_dir = os.path.dirname(genome_file_path)
    genome_file_name = os.path.basename(genome_file_path)
    os.chdir(GRF_working_dir)
    GRF(genome_file_name, genome_name, processors, TIR_length, GRF_path)
    os.chdir("../")


def processors_allocation(processors):
    if processors <= 16:
        num_process = int(math.sqrt(processors))
    else:
        # Divide by sqrt(3) ≈ 1.732 for larger processor counts to favor more threads per process
        num_process = int(math.sqrt(processors/math.sqrt(3)))
    num_thread_per_process = processors // num_process
    total_num_threads = num_process * num_thread_per_process
    return num_process, num_thread_per_process, total_num_threads


def run_GRF_native(genome_file, genome_name, processors, TIR_length, flag_debug, GRF_path):
    print("  Step 1/2: Executing GRF in native mode\n")
    GRF(genome_file, genome_name, processors, TIR_length, GRF_path)
    print("  Step 2/2: Getting GRF result")
    return get_GRF_result_df_native(genome_name, flag_debug)


def run_GRF_py_para(genome_file, genome_name, TIR_length, processors, flag_debug, GRF_path, seq_num,
                    fasta_files_path_list):
    os.makedirs(f"{SPLIT_FASTA_TAG}_mp", exist_ok=True)  # TODO revise to include mode detection function
    os.chdir(f"./{SPLIT_FASTA_TAG}_mp")

    print("  Step 1/3: Checking processed FASTA files")
    if (len(fasta_files_path_list) == 0 or
            any(not os.path.exists(f) or os.path.getsize(f) == 0 for f in fasta_files_path_list)):
        print("    Processed FASTA files not found / invalid. Re-process FASTA files.")
        fasta_files_path_list.extend(process_fasta(genome_file, TIRVISH_SPLIT_SEQ_LEN, TIRVISH_OVERLAP_SEQ_LEN))

    print("  Step 2/3: Executing GRF with python multiprocessing")
    num_process, num_thread_per_process, _ = processors_allocation(processors)

    mp_args_list = [(file_path, genome_name, num_thread_per_process, TIR_length, GRF_path) for
                    file_path in fasta_files_path_list]

    print()
    with mp.Pool(num_process) as pool:
        pool.starmap(GRF_mp, mp_args_list)
    print()

    print("  Step 3/3: Getting GRF result")
    return get_GRF_result_df_para(fasta_files_path_list, genome_name, flag_debug,
                                  TIRVISH_SPLIT_SEQ_LEN, TIRVISH_OVERLAP_SEQ_LEN)


# TODO serial to parallel
# def get_GRF_result_df_boost(fasta_files_path_list, genome_name, flag_debug, split_seq_len, overlap_seq_len):
#     GRF_result_dir_name = f"{genome_name}_GRFmite"
#     GRF_result_dir_list = [os.path.join(os.path.dirname(file), GRF_result_dir_name) for file in
#                            fasta_files_path_list]
#
#     df_list = []
#     for f in GRF_result_dir_list:
#         try:
#             df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
#                             for rec in SeqIO.parse(os.path.join(f, "candidate.fasta"), "fasta")]
#             df_in = pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype({"len": int})
#
#             id_pattern = r"^(\w+)_split_([\w.]+of\d+):(\d+):(\d+):(\w+):(\w+)$"
#
#             def revise_id(match):
#                 seqid, segment_position, start, end, TIR_pattern, TSD = match.groups()
#                 offset = retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len)
#                 return f"{seqid}:{int(start) + offset}:{int(end) + offset}:{TIR_pattern}:{TSD}"
#
#             df_in["id"] = df_in["id"].str.replace(id_pattern, revise_id, regex=True)
#
#             df_list.append(df_in.copy())
#         except FileNotFoundError:
#             continue
#         if not flag_debug:
#             subprocess.Popen(["rm", "-rf", f])
#     return pd.concat(df_list).drop_duplicates(ignore_index=True).sort_values("id", ignore_index=True)

# TODO parallelize this
def get_GRF_result_df_para(fasta_files_path_list, genome_name, flag_debug, split_seq_len, overlap_seq_len):
    GRF_result_dir_name = f"{genome_name}_GRFmite"
    GRF_result_dir_list = [os.path.join(os.path.dirname(file), GRF_result_dir_name) for file in
                           fasta_files_path_list]

    df_list = []
    for f in GRF_result_dir_list:
        try:
            df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
                            for rec in SeqIO.parse(os.path.join(f, "candidate.fasta"), "fasta")]
            df_in = pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype({"len": int})

            id_pattern = r"^(\w+)_split_([\w.]+of\d+):(\d+):(\d+):(\w+):(\w+)$"

            def revise_id(match):
                seqid, segment_position, start, end, TIR_pattern, TSD = match.groups()
                offset = retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len)
                return f"{seqid}:{int(start) + offset}:{int(end) + offset}:{TIR_pattern}:{TSD}"

            df_in["id"] = df_in["id"].str.replace(id_pattern, revise_id, regex=True)

            df_list.append(df_in.copy())
        except FileNotFoundError:
            continue
        if not flag_debug:
            subprocess.Popen(["rm", "-rf", f])

    if len(df_list) == 0:
        return None

    return pd.concat(df_list).drop_duplicates(ignore_index=True).sort_values("id", ignore_index=True)


def get_GRF_result_df_native(genome_name, flag_debug):
    GRF_result_dir_name = f"{genome_name}_GRFmite"
    GRF_result_file_path = os.path.join(GRF_result_dir_name, "candidate.fasta")

    # df_data_dict = {"id": [rec.id for rec in SeqIO.parse(GRF_result_file_path, "fasta")],
    # df_data_dict = {"id": list(SeqIO.index(GRF_result_file_path, "fasta")),
    #                 "seq": [str(rec.seq) for rec in SeqIO.parse(GRF_result_file_path, "fasta")],
    #                 "len": [len(rec) for rec in SeqIO.parse(GRF_result_file_path, "fasta")]}
    try:
        df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
                        for rec in SeqIO.parse(GRF_result_file_path, "fasta")]
    except FileNotFoundError:
        return None

    if not flag_debug:
        subprocess.Popen(["rm", "-rf", GRF_result_dir_name])
    return pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype({"len": int})


def execute(TIRLearner_instance):
    genome_file = TIRLearner_instance.genome_file_path
    genome_name = TIRLearner_instance.genome_name
    TIR_length = TIRLearner_instance.TIR_length
    processors = TIRLearner_instance.processors
    para_mode = TIRLearner_instance.para_mode
    flag_debug = TIRLearner_instance.flag_debug
    GRF_path = TIRLearner_instance.GRF_path
    additional_args = TIRLearner_instance.additional_args
    seq_num = TIRLearner_instance.genome_file_stat["num"]
    fasta_files_path_list = TIRLearner_instance.split_fasta_files_path_list

    if NO_PARALLEL in additional_args:
        return run_GRF_native(genome_file, genome_name, TIR_length, processors, flag_debug, GRF_path)
    elif para_mode == "gnup":
        raise NotImplementedError()
    else:
        return run_GRF_py_para(genome_file, genome_name, TIR_length, processors, flag_debug, GRF_path,
                               seq_num, fasta_files_path_list)
