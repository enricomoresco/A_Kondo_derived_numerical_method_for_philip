#!/usr/bin/env python3
"""Confronta la sezione rivista con l'originale e, se richiesto, la reinserisce
nel file sorgente della tesi.

Uso:
    python applica_revisione.py CARTELLA_DI_LAVORO            # crea solo diff.html
    python applica_revisione.py CARTELLA_DI_LAVORO --applica  # reinserisce nella tesi

Con --applica il file sorgente viene prima copiato in <file>.bak, poi le righe
della sezione vengono sostituite con target_revised.tex. Si rifiuta di procedere
se il sorgente e' cambiato rispetto a quando la sezione e' stata estratta.
"""
import difflib
import json
import os
import shutil
import sys


def main():
    wd = sys.argv[1] if len(sys.argv) > 1 else "."
    apply = "--applica" in sys.argv
    meta = json.load(open(os.path.join(wd, "meta.json"), encoding="utf-8"))
    orig = open(os.path.join(wd, "target.tex"), encoding="utf-8").read().splitlines()
    new = open(os.path.join(wd, "target_revised.tex"), encoding="utf-8").read().splitlines()

    html = difflib.HtmlDiff(wrapcolumn=90).make_file(
        orig, new, "originale", "rivista", context=True, numlines=2)
    with open(os.path.join(wd, "diff.html"), "w", encoding="utf-8") as f:
        f.write(html)
    ratio = difflib.SequenceMatcher(None, "\n".join(orig), "\n".join(new)).ratio()
    print(f"diff.html creato. Somiglianza con l'originale: {ratio:.0%}")

    if not apply:
        return
    if not meta.get("reinserimento_automatico"):
        sys.exit("La sezione attraversa piu' file: reinserimento automatico non possibile. "
                 "Copia target_revised.tex a mano.")
    src = meta["file_sorgente"]
    try:
        lines = open(src, encoding="utf-8").read().splitlines()
    except UnicodeDecodeError:
        sys.exit(f"{src} non e' in UTF-8: non sovrascrivo nulla. Reinserisci a mano.")
    a, b = meta["riga_inizio"] - 1, meta["riga_fine"]
    if lines[a:b] != orig:
        sys.exit(f"{src} e' cambiato dopo l'estrazione: non sovrascrivo nulla. "
                 "Riestrai la sezione o reinserisci a mano.")
    eol = "\r\n" if b"\r\n" in open(src, "rb").read() else "\n"
    shutil.copy2(src, src + ".bak")
    lines[a:b] = new
    with open(src, "w", encoding="utf-8", newline="") as f:
        f.write(eol.join(lines) + eol)
    print(f"Sezione reinserita in {src} (backup: {os.path.basename(src)}.bak)")


if __name__ == "__main__":
    main()
