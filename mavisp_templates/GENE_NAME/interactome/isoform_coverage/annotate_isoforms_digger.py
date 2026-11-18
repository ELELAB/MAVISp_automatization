import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import requests

# ============================================================
# HARD-CODED PATHS (Homo sapiens only)
# ============================================================
DIGGER_HUMAN_DIR = Path("/data/databases/digger/data/data/Homo sapiens[human]")
DIGGER_DOMAIN_MAP_FILE = Path("/data/databases/digger/domain_to_exon/domain_mapped_to_exons.csv")
DIGGER_DDI_FILE = DIGGER_HUMAN_DIR / "predicted_ddi_ppi.tsv"   # domain_1, domain_2, class (ENTREZ/PFAM)
DIGGER_GENE_INFO_FILE = DIGGER_HUMAN_DIR / "gene_info.csv"     # ENST -> NCBI gene ID

UNIPROT_BASE = "https://rest.uniprot.org"

# ============================================================
# Generic helpers
# ============================================================
def split_uniprot_isoform(acc: str) -> Tuple[str, Optional[str]]:
    m = re.match(r"^([A-Z0-9]+?)(?:-(\d+))?$", str(acc).strip())
    if not m:
        raise ValueError(f"Invalid UniProt accession: {acc}")
    return m.group(1), m.group(2)

def read_table_any(path: Path) -> pd.DataFrame:
    p = str(path)
    if p.endswith(".csv"):
        return pd.read_csv(p, low_memory=False)
    if p.endswith(".tsv") or p.endswith(".txt"):
        return pd.read_csv(p, sep="\t", low_memory=False)
    try:
        return pd.read_csv(p, low_memory=False)
    except Exception:
        return pd.read_csv(p, sep="\t", low_memory=False)

def pick_col(df: pd.DataFrame, candidates: List[str], required: bool = True) -> Optional[str]:
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    if required:
        raise KeyError(f"Missing expected column. Tried {candidates}. Available: {list(df.columns)}")
    return None

def normalize_enst(x: str) -> str:
    x = str(x).strip()
    if not x or x.lower() == "nan":
        return ""
    return x.split(".")[0]

def normalize_pfam(x: str) -> str:
    x = str(x).strip()
    if not x or x.lower() == "nan":
        return ""
    return x.split(".")[0]

def parse_digger_domain_token(token: str) -> Tuple[str, str]:
    """
    DIGGER predicted_ddi_ppi.tsv stores tokens like '7157/PF00870'
    Returns (entrez_id, pfam_id)
    """
    s = str(token).strip()
    if not s or s.lower() == "nan":
        return "", ""
    if "/" not in s:
        return "", normalize_pfam(s)
    entrez, pfam = s.split("/", 1)
    return entrez.strip(), normalize_pfam(pfam)

# ============================================================
# UniProt helpers (sequence authority + ENST xrefs)
# ============================================================
def uniprot_entry_json(acc_or_iso: str) -> dict:
    url = f"{UNIPROT_BASE}/uniprotkb/{acc_or_iso}.json"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return r.json()

def uniprot_fasta_seq(acc_or_iso: str) -> str:
    url = f"{UNIPROT_BASE}/uniprotkb/{acc_or_iso}.fasta"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    return "".join([ln.strip() for ln in r.text.splitlines() if ln and not ln.startswith(">")])

def extract_ensembl_transcripts_from_uniprot_entry(entry_json: dict) -> List[str]:
    ensts = set()
    xrefs = entry_json.get("uniProtKBCrossReferences", []) or entry_json.get("dbReferences", [])
    for xr in xrefs:
        db = xr.get("database") or xr.get("type")
        if str(db).lower() != "ensembl":
            continue

        xid = str(xr.get("id", "")).strip().split(".")[0]
        if xid.startswith("ENST"):
            ensts.add(xid)

        for p in xr.get("properties", []) or []:
            val = str(p.get("value", "")).strip().split(".")[0]
            if val.startswith("ENST"):
                ensts.add(val)

    return sorted(ensts)

def get_uniprot_to_enst_from_uniprot_api(uniprot_ids: List[str]) -> pd.DataFrame:
    rows = []
    for uid in sorted(set(uniprot_ids)):
        try:
            entry = uniprot_entry_json(uid)
            ensts = extract_ensembl_transcripts_from_uniprot_entry(entry)
            if ensts:
                rows.append({"uniprot_id": uid, "enst": ";".join(ensts)})
        except Exception:
            continue

    if not rows:
        return pd.DataFrame(columns=["uniprot_id", "enst"])
    return pd.DataFrame(rows)

# ============================================================
# DIGGER loaders
# ============================================================
def load_digger_domain_map_human() -> pd.DataFrame:
    d = read_table_any(DIGGER_DOMAIN_MAP_FILE)
    enst_col = pick_col(d, ["Transcript stable ID"])
    pfam_col = pick_col(d, ["Pfam ID"])

    out = pd.DataFrame({
        "enst": d[enst_col].astype(str).map(normalize_enst),
        "pfam": d[pfam_col].astype(str).map(normalize_pfam),
    })
    out = out[(out["enst"] != "") & (out["pfam"] != "")]
    return out.drop_duplicates()

def load_enst_to_entrez_from_gene_info() -> Dict[str, str]:
    g = read_table_any(DIGGER_GENE_INFO_FILE)
    enst_col = pick_col(g, ["Transcript stable ID"])
    entrez_col = pick_col(g, ["NCBI gene ID"])

    mapping = {}
    for _, r in g.iterrows():
        enst = normalize_enst(r[enst_col])
        gid = str(r[entrez_col]).strip()
        if not enst or not gid or gid.lower() == "nan":
            continue
        gid_clean = gid.split(".")[0]
        if gid_clean.isdigit():
            mapping.setdefault(enst, gid_clean)
    return mapping

def load_digger_ddi_rules_oriented_human() -> Dict[Tuple[str, str], Set[Tuple[str, str]]]:
    """
    Oriented DDI rules:
      (entrezA, entrezB) -> set of (pfamA, pfamB)

    This preserves direction so we can report DIGGER-style target-pfam -> partner-pfam DDIs.
    Also stores the reverse orientation.
    """
    d = read_table_any(DIGGER_DDI_FILE)
    c1 = pick_col(d, ["domain_1"])
    c2 = pick_col(d, ["domain_2"])

    rules: Dict[Tuple[str, str], Set[Tuple[str, str]]] = {}
    for _, r in d.iterrows():
        e1, p1 = parse_digger_domain_token(r[c1])
        e2, p2 = parse_digger_domain_token(r[c2])
        if not e1 or not e2 or not p1 or not p2:
            continue

        rules.setdefault((e1, e2), set()).add((p1, p2))
        rules.setdefault((e2, e1), set()).add((p2, p1))

    return rules

# ============================================================
# Core logic helpers
# ============================================================
def domain_set_for_enst(dom_map: pd.DataFrame, enst: str) -> Set[str]:
    return set(dom_map.loc[dom_map["enst"] == enst, "pfam"].dropna().astype(str))

def choose_best_enst(enst_candidates: List[str], dom_map: pd.DataFrame, enst2entrez: Dict[str, str]) -> str:
    """
    Pick the most useful ENST:
    1) present in domain map and has Entrez
    2) present in domain map
    3) has Entrez
    4) first available
    """
    if not enst_candidates:
        return ""

    domain_ensts = set(dom_map["enst"].unique())

    for e in enst_candidates:
        if e in domain_ensts and e in enst2entrez:
            return e
    for e in enst_candidates:
        if e in domain_ensts:
            return e
    for e in enst_candidates:
        if e in enst2entrez:
            return e
    return enst_candidates[0]

def classify_confidence(missing_pct: float, has_digger_edge: bool) -> str:
    """
    Practical DIGGER-like confidence approximation.
    """
    if not has_digger_edge:
        return "Uncertain"
    if missing_pct == 0:
        return "Original"
    return "Low"

# ============================================================
# Main
# ============================================================
def annotate_isoform_digger_human(
    aggregated_csv: str,
    target_isoform_uniprot: str,
    output_csv: str,
) -> pd.DataFrame:
    for p in [DIGGER_DOMAIN_MAP_FILE, DIGGER_DDI_FILE, DIGGER_GENE_INFO_FILE]:
        if not p.exists():
            raise FileNotFoundError(f"Missing required DIGGER file: {p}")

    df = pd.read_csv(aggregated_csv, low_memory=False)
    for c in ["Target_Uniprot_AC", "Target_protein", "Interactor_UniProt_AC", "Interactor"]:
        if c not in df.columns:
            raise KeyError(f"Input CSV missing required column: {c}")

    target_base, iso_suffix = split_uniprot_isoform(target_isoform_uniprot)
    if iso_suffix is None:
        raise ValueError("Please provide target isoform like P04637-4")
    if target_base not in set(df["Target_Uniprot_AC"].astype(str).unique()):
        raise ValueError(f"Target base {target_base} not found in Target_Uniprot_AC column")

    # Load DIGGER evidence
    dom_map = load_digger_domain_map_human()
    enst2entrez = load_enst_to_entrez_from_gene_info()
    ddi_rules_oriented = load_digger_ddi_rules_oriented_human()

    # UniProt target metadata
    target_base_entry = uniprot_entry_json(target_base)
    target_ref_len = int(target_base_entry.get("sequence", {}).get("length", 0))
    target_iso_len = len(uniprot_fasta_seq(target_isoform_uniprot))

    # Build UniProt->ENST mapping from UniProt xrefs (target + all partner base IDs)
    interactor_bases = sorted(set(df["Interactor_UniProt_AC"].astype(str)))
    ids_to_map = sorted({target_base, target_isoform_uniprot, *interactor_bases})
    uniprot_enst = get_uniprot_to_enst_from_uniprot_api(ids_to_map)

    # ---- Target ENST mapping (smart selection) ----
    target_iso_enst_candidates = []
    target_ref_enst_candidates = []

    row_target_iso = uniprot_enst[uniprot_enst["uniprot_id"] == target_isoform_uniprot]
    if not row_target_iso.empty and str(row_target_iso.iloc[0]["enst"]).strip():
        target_iso_enst_candidates = [x for x in str(row_target_iso.iloc[0]["enst"]).split(";") if x.startswith("ENST")]

    row_target_base = uniprot_enst[uniprot_enst["uniprot_id"] == target_base]
    if not row_target_base.empty and str(row_target_base.iloc[0]["enst"]).strip():
        target_ref_enst_candidates = [x for x in str(row_target_base.iloc[0]["enst"]).split(";") if x.startswith("ENST")]

    target_ref_enst = choose_best_enst(target_ref_enst_candidates, dom_map, enst2entrez)
    target_iso_enst = choose_best_enst(target_iso_enst_candidates, dom_map, enst2entrez)

    if not target_iso_enst and target_ref_enst:
        target_iso_enst = target_ref_enst

    # Target domains / target Entrez
    target_domains = domain_set_for_enst(dom_map, target_iso_enst) if target_iso_enst else set()
    ref_domains = domain_set_for_enst(dom_map, target_ref_enst) if target_ref_enst else set()
    target_entrez = enst2entrez.get(target_iso_enst, "") if target_iso_enst else ""

    if not target_iso_enst or not target_ref_enst:
        target_domain_status = "no_digger_mapping"
        lost_domains, gained_domains = [], []
    else:
        lost_domains = sorted(ref_domains - target_domains)
        gained_domains = sorted(target_domains - ref_domains)
        if lost_domains and gained_domains:
            target_domain_status = "domain_rewired"
        elif lost_domains:
            target_domain_status = "domain_loss"
        elif gained_domains:
            target_domain_status = "domain_gain"
        else:
            target_domain_status = "no_domain_change_detected"

    # Replace target accession in the final output with the isoform accession
    df["Target_Uniprot_AC"] = target_isoform_uniprot

    # Fill target output columns
    df["Target_Reference_Seq_Length"] = target_ref_len
    df["Target_Isoform_Seq_Length"] = target_iso_len
    df["Target_DIGGER_Mapped_ENST"] = target_iso_enst
    df["Target_DIGGER_Domain_Status"] = target_domain_status
    df["Target_DIGGER_Lost_Domains"] = ";".join(lost_domains)
    df["Target_DIGGER_Gained_Domains"] = ";".join(gained_domains)

    # DIGGER-style partner-level output columns
    df["DIGGER_Partner_NCBI_Gene_ID"] = ""
    df["DIGGER_Retained_DDIs"] = ""
    df["DIGGER_Missing_DDIs"] = ""
    df["DIGGER_Missing_DDI_Percent"] = pd.NA
    df["DIGGER_PPI_Effect"] = "Uncertain"
    df["DIGGER_Confidence"] = "Uncertain"

    # Row-wise DIGGER-style partner effect
    for idx, row in df.iterrows():
        inter_base = str(row["Interactor_UniProt_AC"]).strip()

        if not target_iso_enst or not target_entrez:
            continue

        # Partner ENST from UniProt base accession
        rmap = uniprot_enst[uniprot_enst["uniprot_id"] == inter_base]
        if rmap.empty:
            continue

        inter_enst_candidates = [x for x in str(rmap.iloc[0]["enst"]).split(";") if x.startswith("ENST")]
        inter_enst = choose_best_enst(inter_enst_candidates, dom_map, enst2entrez)
        if not inter_enst:
            continue

        inter_entrez = enst2entrez.get(inter_enst, "")
        if not inter_entrez:
            continue

        pair_key = (target_entrez, inter_entrez)  # oriented target -> partner
        oriented_ddis = ddi_rules_oriented.get(pair_key, set())

        if not oriented_ddis:
            # Mentha/STRING edge exists, but DIGGER has no gene-specific DDI rule
            df.at[idx, "DIGGER_PPI_Effect"] = "Uncertain"
            df.at[idx, "DIGGER_Confidence"] = "Uncertain"
            continue

        retained = []
        missing = []

        # DIGGER-style split based on target isoform target-side domain presence
        # A DDI is missing if target-side Pfam is absent in isoform but present in reference
        # A DDI is retained if target-side Pfam is present in isoform
        for tpf, ppf in sorted(oriented_ddis):
            ddi_label = f"{tpf}-{ppf}"

            if tpf in target_domains:
                retained.append(ddi_label)
            elif tpf in ref_domains:
                missing.append(ddi_label)
            else:
                # Target-side domain not in ref or isoform (rare/noisy rule relative to chosen ENST)
                # ignore for percent calculation
                pass

        retained = sorted(set(retained))
        missing = sorted(set(missing))

        total_considered = len(retained) + len(missing)
        if total_considered == 0:
            missing_pct = pd.NA
            effect = "Uncertain"
            confidence = "Uncertain"
        else:
            missing_pct_val = round((len(missing) / total_considered) * 100, 2)
            missing_pct = missing_pct_val
            effect = "Affected" if len(missing) > 0 else "Retained"
            confidence = classify_confidence(missing_pct_val, has_digger_edge=True)

        df.at[idx, "DIGGER_Partner_NCBI_Gene_ID"] = inter_entrez
        df.at[idx, "DIGGER_Retained_DDIs"] = ", ".join(retained)
        df.at[idx, "DIGGER_Missing_DDIs"] = ", ".join(missing)
        df.at[idx, "DIGGER_Missing_DDI_Percent"] = missing_pct
        df.at[idx, "DIGGER_PPI_Effect"] = effect
        df.at[idx, "DIGGER_Confidence"] = confidence

    df.to_csv(output_csv, index=False)
    return df

# ============================================================
# CLI
# ============================================================
if __name__ == "__main__":
    import sys

    if len(sys.argv) != 4:
        print("Usage: python annotate_isoforms_digger.py <aggregated_csv> <target_isoform_uniprot> <output_csv>")
        raise SystemExit(1)

    aggregated_csv = sys.argv[1]
    target_isoform_uniprot = sys.argv[2]
    output_csv = sys.argv[3]

    out = annotate_isoform_digger_human(
        aggregated_csv=aggregated_csv,
        target_isoform_uniprot=target_isoform_uniprot,
        output_csv=output_csv,
    )

    print(f"[INFO] Saved: {output_csv}")
