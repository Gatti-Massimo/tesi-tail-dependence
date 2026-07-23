# Istruzioni Codex per questa tesi

## Scopo del progetto

Questo repository serve ora a capire i tre paper sulla tail dependence e, più avanti, a sviluppare codice R e analisi. La struttura deve restare minima.

## Regole non negoziabili

- I tre PDF in `papers/` sono sorgenti originali: non modificarli, rinominarli, annotarli o sovrascriverli.
- Gli appunti personali sono esclusivamente file LaTeX in `appunti/`, uno per paper: `schmidt.tex`, `durante.tex`, `hua_joe.tex`.
- Gli appunti devono spiegare i paper in modo didattico per Massimo; non sono testo da inviare al relatore né capitoli definitivi della tesi.
- Non creare o compilare PDF senza una richiesta esplicita dell'utente.
- Non creare commit, non fare push e non modificare la storia Git senza una richiesta esplicita dell'utente.
- Mantenere la struttura minimale: creare nuove cartelle solo quando servono davvero.
- Il codice R va in `codice/`; i grafici generati vanno in `grafici/`.
- Non aggiungere README, workflow, bibliografie o file organizzativi salvo richiesta esplicita.

## Instradamento dei task

Quando il task corrisponde chiaramente a una riga, il coordinatore deve usare il subagent indicato. Per task composti, seguire l'ordine indicato e non far modificare gli stessi file a più agenti contemporaneamente.

| Task richiesto | Subagent | Azione autorizzata |
|---|---|---|
| Leggere, estrarre o confrontare uno o più paper | `paper_reader` | sola lettura; restituisce una sintesi strutturata al coordinatore |
| Controllare formule, definizioni o interpretazioni contro i PDF | `formula_reviewer` | sola lettura; segnala errori e ambiguità, non modifica file |
| Scrivere o aggiornare gli appunti di un paper | `latex_notes_writer` | modifica solo il file `.tex` del paper in `appunti/`; non compila PDF |
| Scrivere codice R, simulazioni o stime | `r_analyst` | modifica solo `codice/` e, se richiesto, produce file in `grafici/` |
| Controllare diff, creare commit o push | `git_publisher` | opera solo dopo autorizzazione esplicita e sul perimetro confermato |

### Sequenze consigliate

- Comprensione di un paper: `paper_reader` → eventualmente `formula_reviewer` → `latex_notes_writer`.
- Modifica di appunti esistenti: `formula_reviewer` → `latex_notes_writer` se l'utente chiede anche la correzione.
- Analisi empirica: `r_analyst`; chiedere `formula_reviewer` solo per verificare la traduzione matematica del paper.
- Git: `git_publisher` soltanto quando l'utente dice esplicitamente di fare commit o push.

## Comportamento del coordinatore

- Prima di delegare, chiarire quale file o paper è in scope.
- Restituire all'utente il risultato del subagent in italiano, distinguendo fatti dal paper, spiegazioni intuitive e interpretazioni.
- Non delegare task banali come `git status`, elenco file o una singola spiegazione già contenuta nel contesto.
- Preferire un solo subagent ben scelto; usare più agenti solo per fasi indipendenti o per una verifica mirata.
- Se un subagent propone una modifica fuori dal proprio perimetro, fermarlo e riportare il problema al coordinatore.
