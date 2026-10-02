#!/usr/bin/env python3
"""Verifica meccanica delle voci citate nella sezione contro Crossref.

Uso:
    python verifica_doi.py CARTELLA_DI_LAVORO

Legge citazioni.txt e i .bib indicati in meta.json, e scrive report_doi.md.
Per ogni chiave controlla:
  - che esista nel .bib;
  - se ha un DOI, che il DOI esista e che titolo, anno e primo autore coincidano;
  - se non ha DOI, cerca la voce su Crossref e propone il DOI piu' probabile.
Libri, tesi e report spesso non sono su Crossref: in quel caso la voce e'
marcata "da verificare a mano", non come errore.

Variabile d'ambiente opzionale CROSSREF_MAILTO: un indirizzo email da passare
a Crossref per usare il loro "polite pool" (risposte piu' rapide).
"""
import difflib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

API = "https://api.crossref.org/works"


def parse_bib(text):
    entries = {}
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", text):
        typ, key = m.group(1).lower(), m.group(2)
        if typ in ("comment", "string", "preamble"):
            continue
        depth, i = 1, m.end()
        while i < len(text) and depth:
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            i += 1
        body = text[m.end():i - 1]
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*", body):
            name, j = fm.group(1).lower(), fm.end()
            if j >= len(body):
                continue
            if body[j] == "{":
                d, k = 0, j
                while k < len(body):
                    d += {"{": 1, "}": -1}.get(body[k], 0)
                    if d == 0:
                        break
                    k += 1
                val = body[j + 1:k]
            elif body[j] == '"':
                k = body.find('"', j + 1)
                val = body[j + 1:k]
            else:
                val = re.match(r"[^,\n]*", body[j:]).group(0)
            fields.setdefault(name, val.strip())
        entries[key] = {"type": typ, **fields}
    return entries


def clean(s):
    s = re.sub(r"\\[a-zA-Z]+\s*", "", s or "")
    s = re.sub(r"[{}$\\]", "", s)
    s = re.sub(r"[^\w\s]", " ", s.lower())
    return re.sub(r"\s+", " ", s).strip()


def first_author(authors):
    a = (authors or "").split(" and ")[0].strip()
    a = a.split(",")[0] if "," in a else (a.split()[-1] if a.split() else "")
    return clean(a)


def fetch(url):
    mail = os.environ.get("CROSSREF_MAILTO")
    if mail:
        url += ("&" if "?" in url else "?") + "mailto=" + urllib.parse.quote(mail)
    req = urllib.request.Request(url, headers={"User-Agent": "revisiona-tesi/1.0"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=20) as r:
                return json.load(r)["message"]
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if e.code in (429, 503):
                time.sleep(2 * (attempt + 1))
                continue
            raise
    raise RuntimeError("Crossref non risponde")


def compare(entry, item):
    issues = []
    ct = clean((item.get("title") or [""])[0])
    sim = difflib.SequenceMatcher(None, clean(entry.get("title")), ct).ratio()
    if sim < 0.85:
        issues.append(f"titolo diverso (somiglianza {sim:.2f}): Crossref = \"{(item.get('title') or [''])[0]}\"")
    parts = (item.get("issued") or {}).get("date-parts") or [[None]]
    year = parts[0][0]
    if year and entry.get("year") and str(year) != re.sub(r"\D", "", entry["year"])[:4]:
        issues.append(f"anno diverso: .bib {entry['year']}, Crossref {year}")
    authors = item.get("author") or []
    if authors and entry.get("author"):
        fam = clean(authors[0].get("family", ""))
        if fam and fam != first_author(entry["author"]):
            issues.append(f"primo autore diverso: .bib {first_author(entry['author'])}, Crossref {fam}")
    return issues, sim


def main():
    wd = sys.argv[1] if len(sys.argv) > 1 else "."
    meta = json.load(open(os.path.join(wd, "meta.json"), encoding="utf-8"))
    keys = [k.strip() for k in open(os.path.join(wd, "citazioni.txt"), encoding="utf-8") if k.strip()]

    entries = {}
    for b in meta.get("bib", []):
        try:
            entries.update(parse_bib(open(b, encoding="utf-8", errors="replace").read()))
        except FileNotFoundError:
            print(f"ATTENZIONE: .bib non trovato: {b}", file=sys.stderr)

    rows, details = [], []
    for key in keys:
        e = entries.get(key)
        if e is None:
            rows.append((key, "ERRORE", "chiave assente dal .bib"))
            continue
        doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", e.get("doi", "")).strip()
        try:
            if doi:
                item = fetch(f"{API}/{urllib.parse.quote(doi)}")
                if item is None:
                    rows.append((key, "ERRORE", f"DOI inesistente su Crossref: {doi}"))
                    continue
                issues, _ = compare(e, item)
                if issues:
                    rows.append((key, "CONTROLLARE", "; ".join(issues)))
                else:
                    rows.append((key, "OK", f"DOI {doi}"))
            else:
                q = urllib.parse.quote(f"{e.get('title', '')} {first_author(e.get('author'))}")
                res = fetch(f"{API}?query.bibliographic={q}&rows=3")
                best, best_sim = None, 0
                for it in (res or {}).get("items", []):
                    _, sim = compare(e, it)
                    if sim > best_sim:
                        best, best_sim = it, sim
                if best and best_sim >= 0.9:
                    issues, _ = compare(e, best)
                    note = f"DOI mancante nel .bib; Crossref propone {best.get('DOI')}"
                    if issues:
                        note += " (" + "; ".join(issues) + ")"
                    rows.append((key, "AGGIUNGERE DOI", note))
                else:
                    rows.append((key, "A MANO", f"nessun DOI e nessuna corrispondenza sicura su Crossref ({e['type']})"))
        except Exception as ex:  # rete assente, timeout, ecc.
            rows.append((key, "NON VERIFICATO", f"errore di rete: {ex}"))
        time.sleep(0.2)

    order = {"ERRORE": 0, "CONTROLLARE": 1, "AGGIUNGERE DOI": 2, "NON VERIFICATO": 3, "A MANO": 4, "OK": 5}
    rows.sort(key=lambda r: order[r[1]])
    with open(os.path.join(wd, "report_doi.md"), "w", encoding="utf-8") as f:
        f.write("# Verifica bibliografica automatica (Crossref)\n\n")
        f.write("Controllo solo formale: esistenza delle voci e coerenza dei metadati. "
                "Non dice nulla su cosa affermano le fonti.\n\n")
        f.write("| Chiave | Esito | Dettagli |\n|---|---|---|\n")
        for k, s, d in rows:
            f.write(f"| `{k}` | {s} | {d.replace('|', '/')} |\n")
    counts = {}
    for _, s, _ in rows:
        counts[s] = counts.get(s, 0) + 1
    print("Voci controllate:", len(rows), "-", ", ".join(f"{k}: {v}" for k, v in counts.items()))


if __name__ == "__main__":
    main()
