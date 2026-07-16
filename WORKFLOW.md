# Workflow di tesi

## 1. Principi

1. **GitHub è la fonte ufficiale.** Non mantenere versioni definitive separate su Desktop, email o Overleaf.
2. **Codice prima dei risultati manuali.** Figure e tabelle devono essere prodotte da script R quando possibile.
3. **Un commit è un'unità logica.** Deve essere possibile capire, dal nome, cosa cambia e perché.
4. **La tesi deve compilare prima di un commit LaTeX importante.**
5. **Gli script devono essere rieseguibili da una sessione R vuota.**

## 2. Giornata di lavoro tipica

All'inizio:

1. aprire `Tesi.Rproj`;
2. controllare la scheda Git di RStudio;
3. eseguire `Pull` prima di modificare file, soprattutto se è stato usato Overleaf;
4. lavorare su un solo obiettivo coerente.

Prima di concludere:

1. salvare tutti i file;
2. rieseguire gli script modificati da una sessione pulita;
3. compilare `thesis/Tesi_tail_dependence.tex` se sono cambiati file LaTeX;
4. esaminare il diff;
5. creare un commit con messaggio informativo;
6. eseguire `Push`.

## 3. Quando creare un commit

Creare un commit quando è completa una modifica autonoma, per esempio:

- completata la scheda di una sezione di un paper;
- aggiunta una definizione con spiegazione verificata;
- introdotto uno script R funzionante;
- aggiunta o corretta una figura riproducibile;
- completata una sottosezione LaTeX che compila;
- modificata la struttura del progetto o una dipendenza.

Non creare un commit:

- per ogni singola frase;
- con file non salvati o output non verificati;
- mescolando appunti, analisi e correzioni LaTeX non correlate;
- usando messaggi come `update`, `modifiche`, `test`, `varie` o `finale`.

## 4. Formato dei messaggi

Formato:

```text
tipo: descrizione breve e specifica
```

Tipi ammessi:

- `docs:` appunti, testo della tesi, bibliografia;
- `analysis:` nuova analisi statistica o modifica sostanziale;
- `feat:` nuova funzionalità riutilizzabile nel codice;
- `fix:` correzione di errore nel codice, nei risultati o nel testo;
- `data:` acquisizione, pulizia o documentazione dei dati;
- `refactor:` riorganizzazione senza cambiamento dei risultati;
- `chore:` configurazione, dipendenze e manutenzione;
- `results:` figure, tabelle o risultati validati.

Esempi corretti:

```text
docs: riassumi definizioni di tail dependence in Schmidt
analysis: confronta coefficienti di coda per copule gaussiane e t
fix: correggi convenzione tra rendimenti e perdite
results: aggiungi figura della tail concentration function
chore: inizializza ambiente riproducibile con renv
```

La descrizione deve:

- usare l'imperativo presente;
- iniziare con lettera minuscola;
- non terminare con un punto;
- restare idealmente entro 72 caratteri;
- descrivere il contenuto effettivo del commit.

Per modifiche non ovvie aggiungere un corpo, separato da una riga vuota, che spieghi motivazione, scelte e limiti.

## 5. Branch

Per il lavoro ordinario individuale si usa `main`, con commit piccoli e verificati.

Creare un branch soltanto per:

- esperimenti che potrebbero essere scartati;
- grandi ristrutturazioni;
- confronto tra modelli alternativi;
- modifiche rischiose al template LaTeX.

Nomi consigliati:

```text
analysis/copule-archimedee
experiment/tail-order-estimator
docs/template-universita
```

## 6. R e riproducibilità

- Non usare `setwd()` negli script: il progetto RStudio definisce la radice.
- Non affidarsi a oggetti presenti nel Global Environment.
- Non salvare `.RData` o `.Rhistory`.
- Usare percorsi relativi alla radice del progetto.
- Numerare gli script eseguiti in sequenza: `01_import.R`, `02_clean.R`, `03_estimate.R`.
- Mettere le funzioni riutilizzabili in `R/` e non duplicarle negli script.
- Registrare le versioni dei pacchetti con `renv`.
- Impostare un seed esplicito per simulazioni e bootstrap.

## 7. LaTeX e Overleaf

Workflow consigliato:

1. modificare i file `.tex` localmente con RStudio, Codex o un editor LaTeX;
2. compilare localmente quando TinyTeX/TeX Live sarà installato;
3. esaminare il PDF prodotto;
4. committare e inviare su GitHub;
5. usare Overleaf come compilatore/collaborazione opzionale, sincronizzato con GitHub.

Se si modifica da Overleaf:

1. sincronizzare Overleaf verso GitHub;
2. eseguire `Pull` in RStudio **prima** di continuare localmente;
3. evitare modifiche simultanee agli stessi file su Overleaf e in locale.

La bibliografia vive in `thesis/bibliography/references.bib`. Figure e tabelle definitive devono avere nomi stabili e descrittivi, senza spazi.

## 8. Ruolo di Codex

Prima di ogni commit Codex deve:

1. mostrare quali file rientrano nella modifica;
2. controllare il diff;
3. eseguire i controlli pertinenti;
4. proporre un messaggio coerente con questa convenzione;
5. evitare di includere file estranei o dati riservati;
6. comunicare hash e messaggio del commit dopo la creazione.
