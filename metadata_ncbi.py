#!/usr/bin/env python3
"""
Fetch metadata from NCBI for a list of protein accessions (e.g. WP_226588834.1).

Per accession it collects:
  - Identical Protein Group (IPG): organisms, strains, assemblies, genomic location
  - Taxonomy lineage (via the assembly's TaxID)
  - BioSample attributes (isolation_source, host, environment, geo location...)
  - Genomic neighborhood FASTA (for antiSMASH / thiopeptide detection)

Usage:
    pip install biopython pandas
    python fetch_ncbi_metadata.py metadata_noactinos.csv numero_access_node_id metadata_enriched.csv
(all three arguments are optional; these are the defaults)
"""


import csv
import io
import os
import sys
import time
from xml.etree import ElementTree as ET

import pandas as pd
from Bio import Entrez

Entrez.email = "erikdavid.f104@gmail.com"   # REQUIRED by NCBI: put your real email
Entrez.api_key = None                      # optional: paste NCBI API key for 10 req/s
DELAY = 0.11 if Entrez.api_key else 0.34   # stay under NCBI rate limits
WINDOW = 25_000                            # bp on each side of the gene
NEIGHBORHOOD_DIR = "neighborhoods"

BIOSAMPLE_KEYS = [
    "isolation_source", "host", "env_broad_scale", "env_local_scale",
    "env_medium", "geo_loc_name", "collection_date", "strain", "lat_lon",
]

_tax_cache, _asm_cache, _bs_cache = {}, {}, {}


def pause():
    time.sleep(DELAY)


def get_ipg(acc):
    """Identical Protein Group report: every genome record carrying this protein."""
    h = Entrez.efetch(db="protein", id=acc, rettype="ipg", retmode="text")
    txt = h.read()
    h.close()
    pause()
    if isinstance(txt, bytes):
        txt = txt.decode()
    rows = list(csv.reader(io.StringIO(txt), delimiter="\t"))
    header = next((r for r in rows if r and r[0] == "Id"), None)
    if header is None:
        return []
    out = []
    for r in rows[rows.index(header) + 1:]:
        if len(r) >= len(header):
            out.append(dict(zip(header, r)))
    return out


def get_assembly(asm_acc):
    """Assembly accession -> BioSample accession + TaxID."""
    if asm_acc in _asm_cache:
        return _asm_cache[asm_acc]
    s = Entrez.read(Entrez.esearch(db="assembly", term=asm_acc))
    pause()
    if not s["IdList"]:
        _asm_cache[asm_acc] = {}
        return {}
    summ = Entrez.read(Entrez.esummary(db="assembly", id=s["IdList"][0]))
    pause()
    d = summ["DocumentSummarySet"]["DocumentSummary"][0]
    res = {
        "biosample": d.get("BioSampleAccn", ""),
        "taxid": str(d.get("Taxid", "")),
        "assembly_level": d.get("AssemblyStatus", ""),        # Complete Genome / Chromosome / Scaffold / Contig
        "refseq_category": d.get("RefSeq_category", ""),      # reference genome / representative genome / na
        "from_type_material": d.get("FromType", ""),          # non-empty = derived from a type strain
        "assembly_name": d.get("AssemblyName", ""),
    }
    _asm_cache[asm_acc] = res
    return res


def get_taxonomy(taxid):
    if not taxid:
        return {}
    if taxid in _tax_cache:
        return _tax_cache[taxid]
    rec = Entrez.read(Entrez.efetch(db="taxonomy", id=taxid, retmode="xml"))[0]
    pause()
    lineage = {x["Rank"]: x["ScientificName"] for x in rec.get("LineageEx", [])}
    res = {
        "scientific_name": rec.get("ScientificName", ""),
        "phylum": lineage.get("phylum", ""),
        "class": lineage.get("class", ""),
        "order": lineage.get("order", ""),
        "family": lineage.get("family", ""),
        "genus": lineage.get("genus", ""),
        "full_lineage": rec.get("Lineage", ""),
    }
    _tax_cache[taxid] = res
    return res


def get_biosample(sample_acc):
    """BioSample attributes: this is where habitat/isolation source lives."""
    if not sample_acc:
        return {}
    if sample_acc in _bs_cache:
        return _bs_cache[sample_acc]
    s = Entrez.read(Entrez.esearch(db="biosample", term=sample_acc))
    pause()
    if not s["IdList"]:
        _bs_cache[sample_acc] = {}
        return {}
    h = Entrez.efetch(db="biosample", id=s["IdList"][0], retmode="xml")
    xml = h.read()
    h.close()
    pause()
    root = ET.fromstring(xml)
    attrs = {}
    for a in root.iter("Attribute"):
        key = (a.get("harmonized_name") or a.get("attribute_name") or "").lower()
        if key in BIOSAMPLE_KEYS and a.text:
            attrs.setdefault(key, a.text.strip())
    title = root.find(".//Description/Title")
    attrs["biosample_title"] = title.text if title is not None else ""
    _bs_cache[sample_acc] = attrs
    return attrs


def save_neighborhood(acc, nuc_acc, start, stop):
    """Download genomic region around the gene, for antiSMASH."""
    os.makedirs(NEIGHBORHOOD_DIR, exist_ok=True)
    path = os.path.join(NEIGHBORHOOD_DIR, acc.replace(".", "_") + ".fasta")
    if os.path.exists(path):
        return path
    lo, hi = sorted((int(start), int(stop)))
    h = Entrez.efetch(
        db="nuccore", id=nuc_acc, rettype="fasta", retmode="text",
        seq_start=max(1, lo - WINDOW), seq_stop=hi + WINDOW,
    )
    with open(path, "w") as fh:
        fh.write(h.read())
    h.close()
    pause()
    return path


def process(acc):
    result = {"accession": acc}
    ipg = get_ipg(acc)
    if not ipg:
        result["note"] = "no IPG record found"
        return result

    # Summary across ALL genomes that carry this identical protein
    result["n_genome_records"] = len(ipg)
    result["organisms_all"] = "; ".join(sorted({r.get("Organism", "") for r in ipg if r.get("Organism")}))
    result["protein_name"] = ipg[0].get("Protein Name", "")
    pname = result["protein_name"].lower()
    result["protein_is_hypothetical"] = (not pname) or ("hypothetical" in pname) or ("uncharacterized" in pname)
    result["protein_db"] = "RefSeq" if acc.startswith("WP_") else "GenBank/INSDC"

    # Pick the best representative: prefer RefSeq rows that have an assembly
    best = next((r for r in ipg if r.get("Source") == "RefSeq" and r.get("Assembly")), None) \
        or next((r for r in ipg if r.get("Assembly")), ipg[0])
    result["organism"] = best.get("Organism", "")
    result["strain"] = best.get("Strain", "")
    result["assembly"] = best.get("Assembly", "")
    result["nuccore"] = best.get("Nucleotide Accession", "")

    if result["assembly"]:
        asm = get_assembly(result["assembly"])
        result["biosample"] = asm.get("biosample", "")
        result.update(get_taxonomy(asm.get("taxid", "")))
        result.update(get_biosample(asm.get("biosample", "")))

    if result["nuccore"] and best.get("Start") and best.get("Stop"):
        result["neighborhood_fasta"] = save_neighborhood(
            acc, result["nuccore"], best["Start"], best["Stop"])
    return result


def main():
    args = sys.argv[1:]
    in_csv = args[0] if len(args) > 0 else "03_Results/conjunto_noactinomicetos_YcaO/metadata_noactinos.csv"
    col = args[1] if len(args) > 1 else "numero_access_node_id"
    out_csv = args[2] if len(args) > 2 else "metadata_enriched.csv"
    df = pd.read_csv(in_csv)
    accs = df[col].dropna().astype(str).str.strip().unique()

    rows = []
    for i, acc in enumerate(accs, 1):
        print(f"[{i}/{len(accs)}] {acc}")
        try:
            rows.append(process(acc))
        except Exception as e:  # keep going on failures
            rows.append({"accession": acc, "note": f"ERROR: {e}"})
        if i % 25 == 0:  # checkpoint
            pd.DataFrame(rows).to_csv(out_csv, index=False)

    meta = pd.DataFrame(rows)
    merged = df.merge(meta, how="left", left_on=col, right_on="accession")
    if "especie" in merged.columns:
        sp = merged["especie"].fillna("")
        org = merged.get("organism", pd.Series("", index=merged.index)).fillna("")
        merged["species_is_named"] = ~sp.str.contains(r"\bsp\.|unclassified|bacterium", case=False, regex=True)
        merged["species_matches_ncbi"] = [o.startswith(" ".join(s_.split()[:2])) if o else None
                                          for s_, o in zip(sp, org)]
    merged.to_csv(out_csv, index=False)
    print(f"Done. Wrote {out_csv}")


if __name__ == "__main__":
    main()