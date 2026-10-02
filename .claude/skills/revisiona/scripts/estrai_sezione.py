#!/usr/bin/env python3
"""Estrae un chapter/section/subsection da una tesi LaTeX (anche divisa in piu'
file con \\input/\\include) e prepara la cartella di lavoro per la revisione.

Uso:
    python estrai_sezione.py MAIN.tex "Titolo o parte del titolo" [--out DIR]
    python estrai_sezione.py MAIN.tex "label:sec:funnel"

Produce in DIR:
    target.tex        la sezione selezionata, copiata identica
    thesis_flat.tex   la tesi intera appiattita (contesto per gli agenti)
    outline.md        indice della tesi con file e riga di ogni titolo
    citazioni.txt     chiavi citate nella sezione
    meta.json         posizione della sezione nei sorgenti e file .bib
"""
import argparse
import datetime
import json
import os
import re
import sys
import unicodedata

LEVELS = {"part": 0, "chapter": 1, "section": 2, "subsection": 3, "subsubsection": 4}
HEAD_RE = re.compile(r"\\(part|chapter|section|subsection|subsubsection)\*?\s*(\[[^\]]*\])?\s*\{")
INPUT_RE = re.compile(r"\\(input|include|subfile)\s*\{([^}]+)\}")
STOP_RE = re.compile(r"\\(end\{document\}|appendix\b|backmatter\b|bibliography\{|printbibliography)")
CITE_RE = re.compile(r"\\[a-zA-Z]*cite[a-zA-Z]*\*?\s*(?:\[[^\]]*\]\s*){0,2}\{([^}]+)\}")
LABEL_RE = re.compile(r"\\label\{([^}]+)\}")
BIB_RE = re.compile(r"\\(?:bibliography|addbibresource)\s*(?:\[[^\]]*\])?\{([^}]+)\}")


def strip_comment(line):
    """Rimuove il commento LaTeX (% non preceduto da backslash)."""
    out = []
    i = 0
    while i < len(line):
        c = line[i]
        if c == "\\" and i + 1 < len(line):
            out.append(line[i:i + 2])
            i += 2
            continue
        if c == "%":
            break
        out.append(c)
        i += 1
    return "".join(out)


def resolve(path, base):
    p = os.path.join(base, path.strip())
    if not os.path.splitext(p)[1]:
        p += ".tex"
    return os.path.normpath(p)


def flatten(path, base, seen=None):
    """Restituisce una lista di (testo_riga, file, numero_riga)."""
    seen = seen or set()
    path = os.path.normpath(path)
    if path in seen:
        return []
    seen.add(path)
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except FileNotFoundError:
        print(f"ATTENZIONE: file incluso non trovato: {path}", file=sys.stderr)
        return []
    out = []
    for n, line in enumerate(lines, 1):
        code = strip_comment(line)
        m = INPUT_RE.search(code)
        if m:
            before = code[:m.start()].strip()
            if before:
                out.append((line, path, n))
            out.extend(flatten(resolve(m.group(2), base), base, seen))
        else:
            out.append((line, path, n))
    return out


def braced(text, start):
    """Contenuto fra graffe bilanciate a partire da text[start] == '{'."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:i]
    return text[start + 1:]


def norm(s):
    s = unicodedata.normalize("NFKD", s)
    s = re.sub(r"\\[a-zA-Z]+\*?", " ", s)
    s = re.sub(r"[{}$\\~^_]", " ", s)
    return re.sub(r"\s+", " ", s).strip().lower()


def find_headings(flat):
    heads = []
    for idx, (line, path, n) in enumerate(flat):
        code = strip_comment(line)
        for m in HEAD_RE.finditer(code):
            title = braced(code, m.end() - 1)
            label = None
            lm = LABEL_RE.search(code[m.end():])
            if not lm and idx + 1 < len(flat):
                nxt = strip_comment(flat[idx + 1][0]).strip()
                if nxt.startswith("\\label"):
                    lm = LABEL_RE.match(nxt)
            if lm:
                label = lm.group(1)
            heads.append({"idx": idx, "level": LEVELS[m.group(1)], "cmd": m.group(1),
                          "title": title.strip(), "file": path, "line": n, "label": label})
    return heads


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("main")
    ap.add_argument("target")
    ap.add_argument("--out")
    ap.add_argument("--struttura", action="store_true",
                    help="cartella revisione/struttura-<slug>-<data> (revisione strutturale)")
    a = ap.parse_args()

    main_path = os.path.abspath(a.main)
    base = os.path.dirname(main_path)
    flat = flatten(main_path, base)
    heads = find_headings(flat)
    if not heads:
        sys.exit("Nessun titolo (\\chapter, \\section, ...) trovato.")

    q = a.target.strip()
    if q.lower().startswith("label:"):
        lab = q[6:].strip()
        cands = [h for h in heads if h["label"] == lab]
    else:
        nq = norm(q)
        cands = [h for h in heads if norm(h["title"]) == nq] or \
                [h for h in heads if nq in norm(h["title"])]
    if not cands:
        print("Nessuna sezione corrisponde. Titoli disponibili:", file=sys.stderr)
        for h in heads:
            print("  " * h["level"] + f"- {h['title']}", file=sys.stderr)
        sys.exit(2)
    if len(cands) > 1:
        print("Piu' sezioni corrispondono, specifica meglio (o usa label:...):", file=sys.stderr)
        for h in cands:
            print(f"  - [{h['cmd']}] {h['title']}  ({h['file']}:{h['line']}, label={h['label']})",
                  file=sys.stderr)
        sys.exit(3)
    t = cands[0]

    end = len(flat)
    for h in heads:
        if h["idx"] > t["idx"] and h["level"] <= t["level"]:
            end = h["idx"]
            break
    for i in range(t["idx"] + 1, end):
        if STOP_RE.search(strip_comment(flat[i][0])):
            end = i
            break
    seg = flat[t["idx"]:end]
    while len(seg) > 1 and not seg[-1][0].strip():
        seg.pop()

    files = sorted({p for _, p, _ in seg})
    single = len(files) == 1 and all(seg[k + 1][2] == seg[k][2] + 1 for k in range(len(seg) - 1))

    slug = re.sub(r"[^a-z0-9]+", "-", norm(t["title"]))[:40].strip("-") or "sezione"
    prefix = "struttura-" if a.struttura else ""
    out = a.out or os.path.join(base, "revisione", f"{prefix}{slug}-{datetime.date.today():%Y%m%d}")
    os.makedirs(out, exist_ok=True)

    target_text = "\n".join(l for l, _, _ in seg) + "\n"
    with open(os.path.join(out, "target.tex"), "w", encoding="utf-8") as f:
        f.write(target_text)
    with open(os.path.join(out, "thesis_flat.tex"), "w", encoding="utf-8") as f:
        f.write("\n".join(l for l, _, _ in flat) + "\n")

    with open(os.path.join(out, "outline.md"), "w", encoding="utf-8") as f:
        f.write("# Indice della tesi\n\n")
        for h in heads:
            mark = "  <== SEZIONE IN REVISIONE" if h is t else ""
            rel = os.path.relpath(h["file"], base)
            lab = f" `{h['label']}`" if h["label"] else ""
            f.write("  " * h["level"] + f"- [{h['cmd']}] {h['title']}{lab} ({rel}:{h['line']}){mark}\n")

    keys = []
    code_only = "\n".join(strip_comment(l) for l, _, _ in seg)
    for m in CITE_RE.finditer(code_only):
        for k in m.group(1).split(","):
            k = k.strip()
            if k and k not in keys:
                keys.append(k)
    with open(os.path.join(out, "citazioni.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(keys) + ("\n" if keys else ""))

    bibs = []
    for l, _, _ in flat:
        for m in BIB_RE.finditer(strip_comment(l)):
            for b in m.group(1).split(","):
                # LaTeX risolve i percorsi rispetto alla cartella del main, non del file incluso
                bp = os.path.join(base, b.strip())
                if not bp.endswith(".bib"):
                    bp += ".bib"
                bp = os.path.normpath(bp)
                if bp not in bibs:
                    bibs.append(bp)

    chap = next((h for h in reversed(heads) if h["idx"] <= t["idx"] and h["level"] <= 1), t)
    meta = {
        "titolo": t["title"], "comando": t["cmd"], "label": t["label"],
        "capitolo": chap["title"],
        "slug_capitolo": re.sub(r"[^a-z0-9]+", "-", norm(chap["title"]))[:40].strip("-"),
        "file_sorgente": files[0] if single else files,
        "riga_inizio": seg[0][2], "riga_fine": seg[-1][2],
        "reinserimento_automatico": single,
        "bib": bibs, "n_citazioni": len(keys), "main": main_path,
    }
    with open(os.path.join(out, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    print(f"Sezione: [{t['cmd']}] {t['title']}")
    print(f"Sorgente: {', '.join(os.path.relpath(x, base) for x in files)} "
          f"righe {seg[0][2]}-{seg[-1][2]} ({len(seg)} righe)")
    print(f"Citazioni nella sezione: {len(keys)}; file .bib: {len(bibs)}")
    if not single:
        print("NOTA: la sezione attraversa piu' file; il reinserimento andra' fatto a mano.")
    print(f"Cartella di lavoro: {out}")


if __name__ == "__main__":
    main()
