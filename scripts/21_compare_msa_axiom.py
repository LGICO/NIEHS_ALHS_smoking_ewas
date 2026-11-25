#!/usr/bin/env python3
"""Cross-match per-sample MSA VCFs against a multisample Axiom VCF.

The script uses rsID where available, falling back to chromosome/position. It
aligns REF/ALT alleles (including strand complements), skips ambiguous or
non-SNV records, and compares hard-call dosages. Optional MSA results,
sample-sheet, and bridge tables establish the expected cross-platform identity.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


BASES = frozenset("ACGT")
COMP = str.maketrans("ACGT", "TGCA")


def opn(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def chrom_norm(value: str) -> str:
    c = re.sub(r"^chr", "", value, flags=re.I).upper()
    return "MT" if c == "M" else c


def site_keys(chrom: str, pos: str, rsid: str, match_mode: str):
    keys = []
    if match_mode in {"auto", "rsid"}:
        keys.extend(("id", x) for x in dict.fromkeys(re.findall(r"rs\d+", rsid)))
    if match_mode in {"auto", "position"}:
        keys.append(("pos", chrom_norm(chrom), int(pos)))
    return keys


def dosage(gt: str):
    gt = gt.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return None
    alleles = re.split(r"[/|]", gt)
    if not alleles or any(a not in {"0", "1"} for a in alleles):
        return None
    return sum(a == "1" for a in alleles)


def parse_info(info: str):
    out = {}
    for item in info.split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
    return out


def read_table(path: Path):
    """Read comma-, tab-, or simple whitespace-delimited metadata."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip()]
    if not lines:
        raise ValueError(f"Empty metadata file: {path}")
    header = lines[0]
    if "\t" in header:
        return list(csv.DictReader(lines, delimiter="\t"))
    if "," in header:
        return list(csv.DictReader(lines, delimiter=","))
    names = re.split(r"\s+", header.strip())
    rows = []
    for line_no, line in enumerate(lines[1:], 2):
        values = re.split(r"\s+", line.strip())
        if len(values) != len(names):
            raise ValueError(
                f"Whitespace table {path}, line {line_no}: expected "
                f"{len(names)} fields, found {len(values)}"
            )
        rows.append(dict(zip(names, values)))
    return rows


def require_columns(rows, columns, label):
    if not rows:
        raise ValueError(f"No data rows in {label}")
    missing = [col for col in columns if col not in rows[0]]
    if missing:
        raise ValueError(f"Missing column(s) in {label}: {', '.join(missing)}")


def one_to_one_map(rows, key_col, value_col, label):
    mapping = {}
    for row in rows:
        key, value = row[key_col].strip(), row[value_col].strip()
        if not key or not value:
            continue
        if key in mapping and mapping[key] != value:
            raise ValueError(
                f"Non-unique {label}: {key_col}={key!r} maps to multiple "
                f"{value_col} values"
            )
        mapping[key] = value
    return mapping


def canonical_sample_id(value: str) -> str:
    # The bridge uses e.g. a550484-4227346-090615-570_C05 while the Axiom
    # VCF uses a550484.4227346.090615.570_C05.
    return value.strip().lower().replace("-", ".")


def build_axiom_id_resolver(sample_ids):
    exact = set(sample_ids)
    canonical = {}
    for sample_id in sample_ids:
        canonical.setdefault(canonical_sample_id(sample_id), []).append(sample_id)

    def resolve(value):
        value = value.strip()
        if value in exact:
            return value, "EXACT"
        hits = canonical.get(canonical_sample_id(value), [])
        if len(hits) == 1:
            return hits[0], "PUNCTUATION_NORMALIZED"
        if len(hits) > 1:
            return "", "AMBIGUOUS_AXIOM_ID"
        return "", "AXIOM_ID_NOT_FOUND"

    return resolve


def resolve_vcf_path(value: str, results_path: Path, msa_vcf_dir=None) -> Path:
    path = Path(value.strip())
    if path.is_absolute():
        return path
    candidates = [
        path,
        results_path.parent / path,
        results_path.parent.parent / path,
        results_path.parent / "msa_vcf" / path.name,
    ]
    if msa_vcf_dir:
        candidates.append(msa_vcf_dir / path.name)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return path


def valid_snv(ref: str, alt: str) -> bool:
    return len(ref) == len(alt) == 1 and ref in BASES and alt in BASES


def read_axiom(path: Path, match_mode: str):
    samples = []
    calls = {}
    with opn(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10 or "," in f[4] or not valid_snv(f[3].upper(), f[4].upper()):
                continue
            fmt = f[8].split(":")
            if "GT" not in fmt:
                continue
            gt_i = fmt.index("GT")
            info = parse_info(f[7])
            # affy2vcf commonly keeps AX-* in ID and stores the dbSNP name in
            # INFO/DBSNP_RS_ID.
            id_text = ";".join((
                f[2], info.get("DBSNP_RS_ID", ""),
                info.get("dbsnp_rs_id", ""),
            ))
            rec = (f[3].upper(), f[4].upper(), {})
            for sample, field in zip(samples, f[9:]):
                parts = field.split(":")
                d = dosage(parts[gt_i] if gt_i < len(parts) else ".")
                if d is not None:
                    rec[2][sample] = d
            for key in site_keys(f[0], f[1], id_text, match_mode):
                calls.setdefault(key, rec)
    if not samples:
        raise ValueError("No samples found in Axiom VCF header")
    return samples, calls


def read_msa(path: Path, min_gs: float, match_mode: str,
             include_implicit_snps: bool = False, allowed_variants=None):
    calls = {}
    with opn(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8 or "," in f[4]:
                continue
            ref, alt = f[3].upper(), f[4].upper()
            if not valid_snv(ref, alt):
                continue
            info = parse_info(f[7])
            probe_id = info.get("Probe_ID", "")
            # Explicit MSA SNP probes are named rs... in Probe_ID. SeSAMe also
            # emits implicit SNP calls from cg probes affected by variants;
            # those are excluded by default for cross-platform identity QC.
            if not include_implicit_snps and not re.match(r"^rs\d+", probe_id):
                continue
            # Depending on SeSAMe/annotation version, the rsID can occur in
            # Probe_ID or in the separate rs_ID annotation.
            id_text = ";".join((
                f[2], probe_id, info.get("rs_ID", ""), info.get("RS_ID", "")
            ))
            d = dosage(info.get("GT", "."))
            try:
                gs = float(info.get("GS", f[5] if f[5] != "." else "nan"))
            except ValueError:
                gs = math.nan
            if d is None or (not math.isnan(gs) and gs <= min_gs):
                continue
            keys = site_keys(f[0], f[1], id_text, match_mode)
            if allowed_variants is not None:
                keys = [
                    key for key in keys
                    if (key[1] if key[0] == "id" else f"{key[1]}:{key[2]}")
                    in allowed_variants
                ]
            for key in keys:
                old = calls.get(key)
                if old is None or (not math.isnan(gs) and gs > old[3]):
                    calls[key] = (ref, alt, d, gs)
    return calls


def orient(msa_ref, msa_alt, ax_ref, ax_alt):
    if (msa_ref, msa_alt) == (ax_ref, ax_alt):
        return False, "same"
    if (msa_ref, msa_alt) == (ax_alt, ax_ref):
        return True, "swapped"
    cr, ca = msa_ref.translate(COMP), msa_alt.translate(COMP)
    if (cr, ca) == (ax_ref, ax_alt):
        return False, "complement"
    if (cr, ca) == (ax_alt, ax_ref):
        return True, "complement_swapped"
    return None, "incompatible"


def aligned_sites(msa_calls, axiom_calls):
    used = set()
    for key, (mr, ma, md, _gs) in msa_calls.items():
        ar = axiom_calls.get(key)
        if ar is None:
            continue
        # The same record may be indexed by both rsID and position.
        identity = id(ar)
        if identity in used:
            continue
        used.add(identity)
        flip, _status = orient(mr, ma, ar[0], ar[1])
        if flip is None:
            continue
        aligned_md = 2 - md if flip else md
        yield key, aligned_md, ar[2], _status


def compare_one(msa_calls, axiom_samples, axiom_calls):
    counts = {s: [0, 0] for s in axiom_samples}  # compared, concordant
    for _key, aligned_md, axiom_site_calls, _orientation in aligned_sites(
        msa_calls, axiom_calls
    ):
        for sample, ad in axiom_site_calls.items():
            counts[sample][0] += 1
            counts[sample][1] += int(aligned_md == ad)
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axiom-vcf", required=True, type=Path)
    ap.add_argument("--msa-vcf-dir", type=Path)
    ap.add_argument("--msa-results", type=Path)
    ap.add_argument("--msa-sample-sheet", type=Path)
    ap.add_argument("--bridge-file", type=Path)
    ap.add_argument(
        "--fallback-sample-file", type=Path,
        help=(
            "Optional table for MSA IDs absent from the main sample sheet; "
            "for example, missing_ev.csv containing Sample_ID and LHID."
        ),
    )
    ap.add_argument("--msa-results-id-col", default="msa_id")
    ap.add_argument(
        "--sheet-msa-id-col", default="auto",
        help="Sample-sheet column matching msa_id; auto checks Sample_Name and Sample_ID.",
    )
    ap.add_argument("--sheet-bridge-id-col", default="IID")
    ap.add_argument("--bridge-sheet-id-col", default="IID")
    ap.add_argument("--bridge-axiom-id-col", default="BestArrayID")
    ap.add_argument("--fallback-msa-id-col", default="Sample_ID")
    ap.add_argument("--fallback-bridge-id-col", default="LHID")
    ap.add_argument("--bridge-fallback-id-col", default="LHID")
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--min-gs", type=float, default=20.0)
    ap.add_argument("--min-sites", type=int, default=50)
    ap.add_argument("--match-threshold", type=float, default=0.90)
    ap.add_argument("--min-margin", type=float, default=0.05)
    ap.add_argument(
        "--include-implicit-msa-snps", action="store_true",
        help="Also use cg-probe-derived SNP calls; excluded by default.",
    )
    ap.add_argument(
        "--variant-allowlist", type=Path,
        help="One rsID (or CHROM:POS in position mode) per line to retain.",
    )
    ap.add_argument(
        "--match-mode", choices=("auto", "rsid", "position"), default="auto",
        help="Use rsID, coordinate, or both. Use rsid when genome builds differ.",
    )
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    allowed_variants = None
    if args.variant_allowlist:
        with args.variant_allowlist.open() as handle:
            allowed_variants = {
                line.strip().split()[0] for line in handle
                if line.strip() and not line.startswith("#")
            }
        if not allowed_variants:
            raise ValueError("--variant-allowlist contains no variants")
        print(f"Loaded {len(allowed_variants)} allowed variants")

    axiom_samples, axiom_calls = read_axiom(args.axiom_vcf, args.match_mode)
    metadata_args = (args.msa_results, args.msa_sample_sheet, args.bridge_file)
    if any(metadata_args) and not all(metadata_args):
        raise SystemExit(
            "--msa-results, --msa-sample-sheet, and --bridge-file must be "
            "provided together"
        )

    msa_inputs = []
    if all(metadata_args):
        result_rows = read_table(args.msa_results)
        sheet_rows = read_table(args.msa_sample_sheet)
        bridge_rows = read_table(args.bridge_file)
        require_columns(
            result_rows, (args.msa_results_id_col, "vcf_path"), "MSA results"
        )
        sheet_msa_id_col = args.sheet_msa_id_col
        if sheet_msa_id_col == "auto":
            result_ids = {
                row[args.msa_results_id_col].strip() for row in result_rows
                if row[args.msa_results_id_col].strip()
            }
            candidates = [
                col for col in ("Sample_Name", "Sample_ID")
                if col in sheet_rows[0]
            ]
            if not candidates:
                raise ValueError(
                    "Automatic sample-sheet ID selection requires Sample_Name or Sample_ID"
                )
            overlap = {
                col: len(result_ids & {row[col].strip() for row in sheet_rows})
                for col in candidates
            }
            sheet_msa_id_col = max(candidates, key=lambda col: overlap[col])
            print(
                f"Using sample-sheet ID column {sheet_msa_id_col!r}: "
                f"matched {overlap[sheet_msa_id_col]}/{len(result_ids)} MSA result IDs"
            )
        require_columns(
            sheet_rows, (sheet_msa_id_col, args.sheet_bridge_id_col),
            "MSA sample sheet",
        )
        require_columns(
            bridge_rows, (args.bridge_sheet_id_col, args.bridge_axiom_id_col),
            "bridge file",
        )
        sheet_map = one_to_one_map(
            sheet_rows, sheet_msa_id_col, args.sheet_bridge_id_col,
            "sample-sheet join",
        )
        bridge_map = one_to_one_map(
            bridge_rows, args.bridge_sheet_id_col, args.bridge_axiom_id_col,
            "bridge join",
        )
        fallback_map = {}
        bridge_fallback_map = {}
        if args.fallback_sample_file:
            fallback_rows = read_table(args.fallback_sample_file)
            require_columns(
                fallback_rows,
                (args.fallback_msa_id_col, args.fallback_bridge_id_col),
                "fallback sample file",
            )
            require_columns(
                bridge_rows,
                (args.bridge_fallback_id_col, args.bridge_axiom_id_col),
                "bridge file fallback join",
            )
            fallback_map = one_to_one_map(
                fallback_rows,
                args.fallback_msa_id_col,
                args.fallback_bridge_id_col,
                "fallback sample join",
            )
            bridge_fallback_map = one_to_one_map(
                bridge_rows,
                args.bridge_fallback_id_col,
                args.bridge_axiom_id_col,
                "bridge fallback join",
            )
        resolve_axiom_id = build_axiom_id_resolver(axiom_samples)
        seen_msa_ids = set()
        for row in result_rows:
            msa_id = row[args.msa_results_id_col].strip()
            if not msa_id:
                continue
            if msa_id in seen_msa_ids:
                raise ValueError(f"Duplicate MSA result ID: {msa_id!r}")
            seen_msa_ids.add(msa_id)
            sheet_iid = sheet_map.get(msa_id, "")
            fallback_bridge_id = fallback_map.get(msa_id, "")
            bridge_axiom_id = bridge_map.get(sheet_iid, "") if sheet_iid else ""
            linkage_method = "SAMPLE_SHEET_IID" if bridge_axiom_id else ""
            if not bridge_axiom_id and fallback_bridge_id:
                bridge_axiom_id = bridge_fallback_map.get(fallback_bridge_id, "")
                if bridge_axiom_id:
                    linkage_method = "FALLBACK_LHID"
            if bridge_axiom_id:
                expected_axiom_id, mapping_status = resolve_axiom_id(bridge_axiom_id)
            elif sheet_iid:
                expected_axiom_id, mapping_status = "", "IID_NOT_IN_BRIDGE"
            elif fallback_bridge_id:
                expected_axiom_id, mapping_status = "", "FALLBACK_ID_NOT_IN_BRIDGE"
            else:
                expected_axiom_id, mapping_status = "", "MSA_ID_NOT_IN_SAMPLE_SHEET"
            vcf_path = resolve_vcf_path(
                row["vcf_path"], args.msa_results, args.msa_vcf_dir
            )
            msa_inputs.append((
                vcf_path, msa_id, row, sheet_iid, fallback_bridge_id,
                linkage_method, bridge_axiom_id, expected_axiom_id,
                mapping_status,
            ))
    else:
        if not args.msa_vcf_dir:
            raise SystemExit(
                "Provide --msa-vcf-dir, or provide the three metadata files"
            )
        msa_paths = sorted(list(args.msa_vcf_dir.glob("*.vcf")) +
                           list(args.msa_vcf_dir.glob("*.vcf.gz")))
        msa_inputs = [
            (
                path, re.sub(r"\.vcf(?:\.gz)?$", "", path.name), {}, "", "",
                "", "", "", "NOT_PROVIDED",
            )
            for path in msa_paths
        ]
    if not msa_inputs:
        raise SystemExit("No MSA VCF inputs found")

    all_rows, best_rows = [], []
    expected_confusion = Counter()
    expected_orientations = Counter()
    expected_loci = defaultdict(lambda: [0, 0])
    for (
        path, msa_id, result_meta, sheet_iid, fallback_bridge_id,
        linkage_method, bridge_axiom_id, expected_axiom_id, mapping_status,
    ) in msa_inputs:
        if not path.is_file():
            raise FileNotFoundError(f"MSA VCF not found for {msa_id!r}: {path}")
        msa_calls = read_msa(
            path, args.min_gs, args.match_mode,
            include_implicit_snps=args.include_implicit_msa_snps,
            allowed_variants=allowed_variants,
        )
        counts = compare_one(msa_calls, axiom_samples, axiom_calls)
        if expected_axiom_id:
            for key, md, ax_site_calls, orientation in aligned_sites(msa_calls, axiom_calls):
                ad = ax_site_calls.get(expected_axiom_id)
                if ad is None:
                    continue
                expected_confusion[(md, ad)] += 1
                expected_orientations[orientation] += 1
                locus = key[1] if key[0] == "id" else f"{key[1]}:{key[2]}"
                expected_loci[locus][0] += 1
                expected_loci[locus][1] += int(md == ad)
        ranked = []
        for axiom_id, (n, same) in counts.items():
            concordance = same / n if n else math.nan
            ranked.append((concordance, n, same, axiom_id))
            all_rows.append({
                "msa_id": msa_id, "expected_axiom_id": expected_axiom_id,
                "axiom_id": axiom_id,
                "n_compared": n, "n_concordant": same,
                "concordance": "NA" if math.isnan(concordance) else f"{concordance:.6f}",
            })
        ranked.sort(key=lambda x: (-(-1 if math.isnan(x[0]) else x[0]), -x[1], x[3]))
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else (math.nan, 0, 0, "")
        margin = best[0] - second[0] if not any(math.isnan(x) for x in (best[0], second[0])) else math.nan
        passed = (best[1] >= args.min_sites and best[0] >= args.match_threshold and
                  (math.isnan(margin) or margin >= args.min_margin))
        genetic_status = "PASS" if passed else (
            "TOO_FEW_SITES" if best[1] < args.min_sites else
            "LOW_CONCORDANCE" if best[0] < args.match_threshold else "AMBIGUOUS_BEST_MATCH"
        )
        expected_match = None if not expected_axiom_id else best[3] == expected_axiom_id
        sex_matches = result_meta.get("sex_concordant", "").strip().upper()
        if mapping_status not in {"EXACT", "PUNCTUATION_NORMALIZED", "NOT_PROVIDED"}:
            overall_status = "MAPPING_UNRESOLVED"
        elif genetic_status != "PASS":
            overall_status = genetic_status
        elif expected_match is False:
            overall_status = "UNEXPECTED_MATCH"
        elif sex_matches == "FALSE":
            overall_status = "SEX_DISCORDANCE"
        else:
            overall_status = "PASS"
        best_rows.append({
            "msa_id": msa_id,
            "sample_sheet_iid": sheet_iid,
            "bridge_axiom_id": bridge_axiom_id,
            "expected_axiom_id": expected_axiom_id,
            "mapping_status": mapping_status,
            "best_axiom_id": best[3],
            "expected_match": "NA" if expected_match is None else str(expected_match).upper(),
            "recorded_sex": result_meta.get("recorded_sex", ""),
            "msa_inferred_sex": result_meta.get("msa_inferred_sex", ""),
            "sex_concordant": result_meta.get("sex_concordant", ""),
            "n_compared": best[1],
            "best_concordance": "NA" if math.isnan(best[0]) else f"{best[0]:.6f}",
            "second_axiom_id": second[3],
            "second_concordance": "NA" if math.isnan(second[0]) else f"{second[0]:.6f}",
            "margin": "NA" if math.isnan(margin) else f"{margin:.6f}",
            "genetic_status": genetic_status,
            "status": overall_status,
            "fallback_bridge_id": fallback_bridge_id,
            "linkage_method": linkage_method,
        })

    def write_tsv(name, rows):
        with (args.output_dir / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    write_tsv("all_pair_concordance.tsv", all_rows)
    write_tsv("best_matches.tsv", best_rows)
    if expected_confusion:
        confusion_rows = []
        for md in range(3):
            for ad in range(3):
                n = expected_confusion[(md, ad)]
                confusion_rows.append({
                    "msa_alt_dosage": md,
                    "axiom_alt_dosage": ad,
                    "n": n,
                    "concordant": str(md == ad).upper(),
                })
        write_tsv("expected_pair_genotype_confusion.tsv", confusion_rows)
        locus_rows = []
        for locus, (n, same) in expected_loci.items():
            locus_rows.append({
                "variant": locus,
                "n_expected_pairs": n,
                "n_concordant": same,
                "concordance": f"{same / n:.6f}" if n else "NA",
            })
        locus_rows.sort(key=lambda row: (float(row["concordance"]), -row["n_expected_pairs"], row["variant"]))
        write_tsv("expected_pair_snp_qc.tsv", locus_rows)
        orientation_rows = [
            {"orientation": orientation, "n": n}
            for orientation, n in sorted(expected_orientations.items())
        ]
        write_tsv("expected_pair_orientation_counts.tsv", orientation_rows)


if __name__ == "__main__":
    main()
