---
name: revisiona
description: Revisione multi-agente della tesi LaTeX. "/revisiona struttura <capitolo>" fa la revisione strutturale di un capitolo intero; "/revisiona <section o subsection>" revisiona e riscrive una sezione con editor, controllo letteratura, controllo stile AI e coordinatore, con approvazione dell'autore prima di ogni modifica.
---
# Revisione della tesi

Argomenti: `$ARGUMENTS`

- `struttura <titolo>` → modalità STRUTTURA su un capitolo.
- `<titolo>` → modalità SEZIONE su una section o subsection.
- Il titolo può essere una parte del titolo o `label:<label>`.
- `--main <file.tex>` facoltativo; altrimenti usa il file nella cartella corrente che contiene `\documentclass`.
- `--pdf <cartella>` facoltativo: cartella con i PDF delle fonti, da passare all'agente letteratura.
- `--auto` facoltativo: salta le approvazioni e applica le proposte raccomandate SÌ.

Gli script sono in `.claude/skills/revisiona/scripts/`. L'utente lavora su Windows: lancia gli script con `python`.
In tutti i prompt agli agenti indica il percorso assoluto della cartella di lavoro (`WD`).
Parla con l'utente in italiano, in modo sintetico.

## Modalità STRUTTURA

1. `python .claude/skills/revisiona/scripts/estrai_sezione.py <main> "<titolo>" --struttura`
   Se lo script elenca più corrispondenze o nessuna, mostra l'elenco all'utente e chiedi quale intende.
2. Agente `editor`, modalità STRUTTURA, su `WD`.
3. Mostra all'utente una sintesi di `report_struttura.md` (decisione per ogni section e i tre interventi più importanti) e chiedi quali decisioni approva. Salva la sua risposta in `WD/decisioni_struttura.md`.
4. Fine. Ricorda che le singole section si rivedono con `/revisiona <section>`.

## Modalità SEZIONE

1. **Estrazione.** `python .claude/skills/revisiona/scripts/estrai_sezione.py <main> "<titolo>"`. Gestisci ambiguità come sopra.
   Se `meta.json` ha comando `chapter`, suggerisci di fare prima la modalità STRUTTURA e di procedere per section.
2. **Contesto strutturale.** Cerca in `revisione/` la cartella `struttura-<slug_capitolo>-*` più recente. Se esiste, passa agli agenti il percorso del suo `report_struttura.md` e di `decisioni_struttura.md`.
3. **Editor.** Agente `editor`, modalità SEZIONE.
4. **Approvazione 1** (salta con `--auto`). Mostra le decisioni TAGLIARE / SPOSTARE / AMPLIARE di `report_edit.md` e chiedi conferma. Salva in `WD/decisioni_edit.md`.
5. **Verifica DOI.** `python .claude/skills/revisiona/scripts/verifica_doi.py <WD>`. Se l'esito è NON VERIFICATO per errori di rete, dillo all'utente e prosegui.
6. **Letteratura e stile, in parallelo.** Lancia nello stesso messaggio gli agenti `letteratura` (con `--pdf` se fornito) e `stile-ai`, indicando `decisioni_edit.md`.
7. **Proposte.** Agente `coordinatore`, FASE PROPOSTE.
8. **Approvazione 2** (con `--auto`, approva tutte le SÌ). Mostra `proposte.md`: prima le questioni per l'autore, poi la tabella. Chiedi quali approvare (accetta risposte come "tutte le SÌ tranne 4 e 9, la 12 sì"). Registra la scelta in `WD/approvate.md`.
9. **Applicazione.** Agente `coordinatore`, FASE APPLICAZIONE, con l'elenco di `approvate.md`.
10. **Diff.** `python .claude/skills/revisiona/scripts/applica_revisione.py <WD>`. Riporta la percentuale di somiglianza e segnala `diff.html`, `changelog.md` ed eventuale `bib_correzioni.md`.
11. **Reinserimento.** Chiedi all'utente se reinserire la sezione nella tesi. Solo con un sì esplicito: `python .claude/skills/revisiona/scripts/applica_revisione.py <WD> --applica` (crea un backup `.bak`). Se `reinserimento_automatico` è false, spiega che va copiata a mano da `target_revised.tex`.

## Regole
- Non modificare mai i file `.tex` o `.bib` della tesi se non al passo 11.
- Se un agente fallisce o produce un report vuoto, dillo e chiedi se rilanciarlo; non proseguire come se avesse funzionato.
- A fine lavoro, una riga di sintesi: quante proposte applicate, quante rifiutate, quanti TODO e citazioni da verificare restano.
