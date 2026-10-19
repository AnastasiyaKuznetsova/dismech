#!/usr/bin/env python3
"""Download an NPC literature corpus from PubMed and PubMed Central.

What it saves (in --out, default ./corpus):
  query.json            the exact query, date and counts (provenance)
  abstracts.jsonl       one record per paper: pmid, title, abstract, journal, year,
                        publication types, MeSH terms, DOI, PMCID and flags
                        (review, case report, clinical trial, retracted)
  fulltext/PMC*.json.gz full text from the PMC open-access subset (BioC JSON)
  fulltext_txt/PMC*.txt the same full text as plain text with section headings
  fulltext_status.tsv   which PMCIDs were downloaded / not open access

Re-running is safe: papers and full texts already on disk are skipped, so you
can update the corpus later or resume after an interruption.

Examples
  python download_npc_corpus.py --email you@scripps.edu
  python download_npc_corpus.py --email you@scripps.edu --no-fulltext
  python download_npc_corpus.py --email you@scripps.edu --query '"Gaucher Disease"[MeSH]' --out corpus_gaucher

Set NCBI_API_KEY in your environment to go faster (10 instead of 3 requests/s).
Only the Python standard library is needed.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIOC = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/{}/unicode"

DEFAULT_QUERY = (
    '("Niemann-Pick Disease, Type C"[MeSH] OR "Niemann-Pick C1 Protein"[MeSH] '
    'OR "Niemann-Pick type C"[tiab] OR "Niemann-Pick disease type C"[tiab] OR "Niemann-Pick C"[tiab] '
    'OR "Niemann-Pick type C1"[tiab] OR "Niemann-Pick disease type C1"[tiab] OR "Niemann-Pick C1"[tiab] '
    'OR "Niemann-Pick type C2"[tiab] OR "Niemann-Pick disease type C2"[tiab] OR "Niemann-Pick C2"[tiab] '
    'OR "NPC intracellular cholesterol transporter"[tiab] '
    'OR NPC1[tiab] OR NPC2[tiab])'
)
# Note: NPC1[tiab] also returns Ebola/filovirus entry papers (NPC1 is the viral
# receptor). They are kept and can be filtered later with the MeSH terms.


class Client:
    def __init__(self, email: str, api_key: str | None):
        self.email, self.api_key = email, api_key
        self.delay = 0.11 if api_key else 0.35
        self._last = 0.0

    def _wait(self):
        gap = time.time() - self._last
        if gap < self.delay:
            time.sleep(self.delay - gap)
        self._last = time.time()

    def get(self, url: str, params: dict | None = None, data: dict | None = None,
            retries: int = 5) -> bytes | None:
        params = dict(params or {})
        if url.startswith(EUTILS):
            params.update(tool="npc_corpus_downloader", email=self.email)
            if self.api_key:
                params["api_key"] = self.api_key
        full = url + ("?" + urllib.parse.urlencode(params) if params and data is None else "")
        body = None
        if data is not None:
            body = urllib.parse.urlencode({**params, **data}).encode()
        for attempt in range(retries):
            self._wait()
            try:
                req = urllib.request.Request(full, data=body,
                                             headers={"User-Agent": f"npc_corpus_downloader ({self.email})"})
                with urllib.request.urlopen(req, timeout=60) as r:
                    return r.read()
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return None
                if e.code in (429, 500, 502, 503, 504):
                    time.sleep(2 ** attempt)
                    continue
                raise
            except (urllib.error.URLError, TimeoutError):
                time.sleep(2 ** attempt)
        raise RuntimeError(f"Failed after {retries} attempts: {url}")


# ----------------------------------------------------------------------------
# Search
# ----------------------------------------------------------------------------
def search_pmids(c: Client, query: str, start_year: int, end_year: int) -> list[str]:
    """ESearch caps at 10,000 ids per query, so search one year at a time."""
    pmids: list[str] = []
    for year in range(start_year, end_year + 1):
        params = {"db": "pubmed", "term": query, "retmax": 9999, "retmode": "json",
                  "datetype": "pdat", "mindate": str(year), "maxdate": str(year)}
        res = json.loads(c.get(EUTILS + "esearch.fcgi", params))["esearchresult"]
        ids, count = res.get("idlist", []), int(res.get("count", 0))
        if count > 9999:
            print(f"  warning: {year} has {count} hits; only 9999 retrieved. Narrow the query.", file=sys.stderr)
        pmids += ids
        if ids:
            print(f"  {year}: {len(ids)}")
    return list(dict.fromkeys(pmids))  # dedupe, keep order


# ----------------------------------------------------------------------------
# Fetch + parse PubMed XML
# ----------------------------------------------------------------------------
def text_of(node) -> str:
    return " ".join("".join(node.itertext()).split()) if node is not None else ""


def parse_article(art: ET.Element) -> dict | None:
    mc = art.find("MedlineCitation")
    if mc is None:
        return None
    a = mc.find("Article")
    pmid = text_of(mc.find("PMID"))

    parts = []
    for ab in a.findall("Abstract/AbstractText") if a is not None else []:
        label = ab.get("Label")
        t = text_of(ab)
        if t:
            parts.append(f"{label}: {t}" if label else t)

    year = text_of(a.find("Journal/JournalIssue/PubDate/Year")) if a is not None else ""
    if not year and a is not None:
        md = text_of(a.find("Journal/JournalIssue/PubDate/MedlineDate"))
        year = md[:4] if md[:4].isdigit() else text_of(a.find("ArticleDate/Year"))

    pub_types = [text_of(p) for p in mc.findall("Article/PublicationTypeList/PublicationType")]
    mesh = [text_of(m) for m in mc.findall("MeshHeadingList/MeshHeading/DescriptorName")]
    ids = {i.get("IdType"): text_of(i) for i in art.findall("PubmedData/ArticleIdList/ArticleId")}
    retracted = ("Retracted Publication" in pub_types) or any(
        cc.get("RefType") == "RetractionIn" for cc in mc.findall("CommentsCorrectionsList/CommentsCorrections"))
    lower = [p.lower() for p in pub_types]

    return {
        "pmid": pmid,
        "title": text_of(a.find("ArticleTitle")) if a is not None else "",
        "abstract": "\n".join(parts),
        "journal": text_of(a.find("Journal/Title")) if a is not None else "",
        "year": int(year) if year.isdigit() else None,
        "publication_types": pub_types,
        "mesh_terms": mesh,
        "doi": ids.get("doi"),
        "pmcid": ids.get("pmc"),
        "is_review": any("review" in p for p in lower),
        "is_case_report": "case reports" in lower,
        "is_clinical_trial": any(p.startswith("clinical trial") for p in lower),
        "is_retracted": retracted,
        "has_abstract": bool(parts),
    }


def fetch_records(c: Client, pmids: list[str], batch: int = 200):
    for i in range(0, len(pmids), batch):
        chunk = pmids[i:i + batch]
        xml = c.get(EUTILS + "efetch.fcgi", params={"db": "pubmed", "retmode": "xml"},
                    data={"id": ",".join(chunk)})
        root = ET.fromstring(xml)
        for art in root.findall("PubmedArticle"):
            rec = parse_article(art)
            if rec:
                yield rec
        print(f"  fetched {min(i + batch, len(pmids))}/{len(pmids)}")


# ----------------------------------------------------------------------------
# Full text (PMC open-access subset, BioC JSON)
# ----------------------------------------------------------------------------
def bioc_to_text(doc_json) -> str:
    collections = doc_json if isinstance(doc_json, list) else [doc_json]
    lines, last_section = [], None
    for col in collections:
        for doc in col.get("documents", []):
            for p in doc.get("passages", []):
                inf = p.get("infons", {})
                section = inf.get("section_type", "")
                if section in ("REF", "COMP_INT", "AUTH_CONT", "ACK_FUND"):
                    continue
                if section != last_section:
                    lines.append(f"\n## {section}")
                    last_section = section
                t = (p.get("text") or "").strip()
                if t:
                    lines.append(t)
    return "\n".join(lines).strip() + "\n"


def fetch_fulltext(c: Client, records: list[dict], out: Path):
    ft_dir, txt_dir = out / "fulltext", out / "fulltext_txt"
    ft_dir.mkdir(exist_ok=True)
    txt_dir.mkdir(exist_ok=True)
    status_path = out / "fulltext_status.tsv"
    done = {}
    if status_path.exists():
        for line in status_path.read_text().splitlines()[1:]:
            pmcid, st = line.split("\t")[:2]
            done[pmcid] = st
    todo = [r for r in records if r.get("pmcid") and r["pmcid"] not in done]
    print(f"Full text: {len(todo)} PMC articles to try ({len(done)} already checked)")
    new_file = not status_path.exists()
    with open(status_path, "a") as st:
        if new_file:
            st.write("pmcid\tstatus\tpmid\n")
        for k, r in enumerate(todo, 1):
            raw = c.get(BIOC.format(r["pmcid"]))
            ok = False
            if raw and raw.strip().startswith((b"[", b"{")):
                try:
                    j = json.loads(raw)
                    text = bioc_to_text(j)
                    if len(text) > 500:
                        with gzip.open(ft_dir / f"{r['pmcid']}.json.gz", "wb") as fh:
                            fh.write(raw)
                        (txt_dir / f"{r['pmcid']}.txt").write_text(text)
                        ok = True
                except json.JSONDecodeError:
                    pass
            st.write(f"{r['pmcid']}\t{'downloaded' if ok else 'not_open_access'}\t{r['pmid']}\n")
            st.flush()
            if k % 50 == 0:
                print(f"  {k}/{len(todo)}")


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True, help="NCBI asks for a contact email")
    ap.add_argument("--query", default=DEFAULT_QUERY)
    ap.add_argument("--out", default="corpus")
    ap.add_argument("--start-year", type=int, default=1950)
    ap.add_argument("--end-year", type=int, default=dt.date.today().year)
    ap.add_argument("--no-fulltext", action="store_true", help="abstracts only")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    c = Client(args.email, os.environ.get("NCBI_API_KEY"))

    print("Searching PubMed...")
    pmids = search_pmids(c, args.query, args.start_year, args.end_year)
    print(f"{len(pmids)} PMIDs match the query")

    abs_path = out / "abstracts.jsonl"
    have = set()
    if abs_path.exists():
        with open(abs_path) as fh:
            have = {json.loads(l)["pmid"] for l in fh if l.strip()}
    new = [p for p in pmids if p not in have]
    print(f"Fetching {len(new)} new records ({len(have)} already on disk)")
    with open(abs_path, "a") as fh:
        for rec in fetch_records(c, new):
            fh.write(json.dumps(rec, ensure_ascii=False) + "\n")

    with open(abs_path) as fh:
        records = [json.loads(l) for l in fh if l.strip()]

    if not args.no_fulltext:
        fetch_fulltext(c, records, out)

    summary = {
        "query": args.query,
        "years": [args.start_year, args.end_year],
        "run_date": dt.datetime.now().isoformat(timespec="seconds"),
        "n_pmids_matching": len(pmids),
        "n_records": len(records),
        "n_with_abstract": sum(r["has_abstract"] for r in records),
        "n_reviews": sum(r["is_review"] for r in records),
        "n_case_reports": sum(r["is_case_report"] for r in records),
        "n_clinical_trials": sum(r["is_clinical_trial"] for r in records),
        "n_retracted": sum(r["is_retracted"] for r in records),
        "n_with_pmcid": sum(bool(r.get("pmcid")) for r in records),
        "n_fulltext_downloaded": len(list((out / "fulltext_txt").glob("*.txt"))) if (out / "fulltext_txt").exists() else 0,
    }
    (out / "query.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    if summary["n_retracted"]:
        print("Retracted papers are flagged (is_retracted). Never use them as evidence.")


if __name__ == "__main__":
    main()
