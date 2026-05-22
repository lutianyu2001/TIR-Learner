#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *

# DTA:8
# DTC:2/3
# DTH:3(twa)
# DTM:7-10
# DTT:2(TA)

TSD = {"DTA": [8],
       "DTC": [3, 2],
       "DTH": [3],
       "DTM": [10, 9, 8, 7],
       "DTT": [2]}


def _compare(tir1: str, tir2: str) -> int:
    d = 0
    for i in range(0, len(tir1)):
        if tir1[i] != tir2[i]:
            d += 1
    return d


def _sliding_window(seq1: str, seq2: str, TSD_length: int) -> Tuple[list[str], list[str]]:
    set1 = []
    set2 = []
    for i in range(0, len(seq1) - TSD_length + 1):
        set1.append(seq1[i:i + TSD_length])
        set2.append(seq2[i:i + TSD_length])
    return set1, set2


def _conserved(superfamily: str, s1: str) -> bool:
    no_motif = ("DTX", "NonTIR")
    if superfamily in no_motif:
        return True
    else:
        motif1 = " "
        pattern = " "
        if superfamily == "DTA":
            # YARNG
            motif1 = s1[0:5]
            pattern = "[CT]A[AG][ATGC]G"
        if superfamily == "DTC":
            # CMCWR
            motif1 = s1[0:5]
            pattern = "CACT[AG]"
        if superfamily == "DTH":
            motif1 = s1[0:4]
            pattern = "G[GA][GC]C"
        if superfamily == "DTM":
            motif1 = s1[0:1]
            pattern = "[GC]"
        if superfamily == "DTT":
            motif1 = s1[0:10]
            pattern = "CT[ATCG][ATCG]CTC[ATCG][ATCG]T"
        if superfamily == "DTE":
            # GGNRM
            motif1 = s1[0:5]
            pattern = "GG[ATCG][AG][AC]"
        if superfamily == "DTR":
            # CACWATG
            motif1 = s1[0:7]
            pattern = "CAC[AT]ATG"
        if superfamily == "DTP":
            # CANRG
            motif1 = s1[0:5]
            pattern = "CA[ATGC][AG]G"
        motif2 = str(Seq(motif1).reverse_complement())
        z1 = bool(re.match(pattern, motif1))
        z2 = bool(re.match(pattern, motif2))
        if z1 or z2:
            return True
        return False


def _get_difference(set1: List[str], set2: List[str]) -> Dict[str, int]:
    tsd_diff = {}
    for i in range(0, len(set1)):
        for j in range(0, len(set2)):
            name = str(i) + ":" + str(j)
            diff = _compare(set1[i], set2[j])
            tsd_diff[name] = diff
    return tsd_diff


def _conserved_DTH(set1: List[str], TSD_dffset: Dict[str, int], l: int) -> bool:
    for i in TSD_dffset:
        if TSD_dffset[i] < l * 0.2:
            s1 = set1[int(i.split(":")[0])]
            if s1 in ("TTA", "TAA"):
                return True
    return False


def _conserved_DTT(set1: List[str], TSD_dffset: Dict[str, int], l: int) -> bool:
    for i in TSD_dffset:
        if TSD_dffset[i] < l * 0.2:
            s1 = set1[int(i.split(":")[0])]
            if s1[0:2] == "TA":
                return True
    return False


def _is_TSD(TSD_dffset: Dict[str, int], l: int) -> bool:
    for i in TSD_dffset:
        if TSD_dffset[i] < l * 0.2:
            return True
    return False


def _check_TIR(x: pd.Series) -> Union[int, float]:
    TIR_MIN_LEN = 10

    family = x.iloc[0]  # Use index since the column name for TIR superfamily may differ
    s = x["seq"][200:-200]
    len_s = len(s)

    if len_s <= 200:
        l_list = list(range(TIR_MIN_LEN, int(len_s / 2)))
    else:
        l_list = list(range(TIR_MIN_LEN, 100))

    for l in l_list:
        s1 = s[0:l]
        s2_ = s[-l:]
        s2 = str(Seq(s2_).reverse_complement())
        d = _compare(s1, s2)
        if d < l * 0.2 and _conserved(family, s1):
            return l
    return np.nan  # np.nan is float


def _check_TSD(x: pd.Series) -> Union[int, float]:
    family = x.iloc[0]  # Use index since the column name for TIR superfamily may differ
    s = x["seq"]
    l = TSD[family]

    for i in l:
        s1 = s[200 - i:200]
        last20 = s[-200:]
        s2 = last20[0:i]
        set1, set2 = _sliding_window(s1, s2, i)
        dff = _get_difference(set1, set2)
        if family == "DTH" and _conserved_DTH(set1, dff, i):
            return i
        elif family == "DTT" and _conserved_DTT(set1, dff, i):
            return i
        elif family != "DTH" and family != "DTT" and _is_TSD(dff, i):
            return i
    return np.nan  # np.nan is float


def _get_TIR(x: pd.Series) -> pd.Series:
    s = x["seq"][200:-200]
    TIR_len = x["TIR_len"]
    return pd.Series([s[0:TIR_len], s[-TIR_len:]])


def _get_TSD(x: pd.Series) -> pd.Series:
    s = x["seq"]
    TSD_len = x["TSD_len"]

    s1tsd = s[200 - TSD_len:200]
    last200 = s[-200:]
    s2tsd = last200[0:TSD_len]
    set1, set2 = _sliding_window(s1tsd, s2tsd, TSD_len)
    tsd_dffset = _get_difference(set1, set2)
    for i in tsd_dffset:
        if tsd_dffset[i] < TSD_len * 0.2:
            seq1 = set1[int(i.split(":")[0])]
            seq2 = set2[int(i.split(":")[1])]
            return pd.Series([seq1, seq2])
    return pd.Series([np.nan * 2])


def _TIR_TSD_percent(seq1: str, seq2: str) -> float:
    d = _compare(seq1, seq2)
    l = len(seq1)
    p = (l - d) / l
    p *= 100
    p = round(p, 2)
    return p


def _process_result(df_in: pd.DataFrame, module: str) -> pd.DataFrame:
    df = df_in.copy()
    df["source"] = module
    df = df.loc[:, ["seqid", "source", "TIR_type", "sstart", "send",
                    "TIR1", "TIR2", "TIR_percent", "TSD1", "TSD2", "TSD_percent", "len"]]
    df = df.rename(columns={"TIR_type": "type"})
    return df


def execute(TIRLearner_instance, module: str) -> pd.DataFrame:
    df = TIRLearner_instance["base"].copy()
    df["len"] = df["end"] - df["start"] + 1
    df = df[df["len"] >= 450].reset_index(drop=True)  # ?

    print("  Step 1/6: Checking TIR")
    df["TIR_len"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(_check_TIR, axis=1)
    print("  Step 2/6: Checking TSD")
    df["TSD_len"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(_check_TSD, axis=1)
    df = df.dropna(ignore_index=True)
    df = df.astype({"TIR_len": int, "TSD_len": int})

    if df.shape[0] == 0:
        return df

    print("  Step 3/6: Retrieving TIR")
    df[["TIR1", "TIR2"]] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(_get_TIR, axis=1)
    print("  Step 4/6: Calculating TIR percentage")
    df["TIR_percent"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: _TIR_TSD_percent(x["TIR1"], str(Seq(x["TIR2"]).reverse_complement())), axis=1)
    print("  Step 5/6: Retrieving TSD")
    df[["TSD1", "TSD2"]] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(_get_TSD, axis=1)
    print("  Step 6/6: Calculating TSD percentage")
    df["TSD_percent"] = df.swifter.progress_bar(TIRLearner_instance.flag_verbose).apply(
        lambda x: _TIR_TSD_percent(x["TSD1"], x["TSD2"]), axis=1)

    df["len"] = df["len"] - 400  # ?
    return _process_result(df, module)
