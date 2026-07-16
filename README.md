# Tesi magistrale - Tail Dependence

Repository di lavoro per la tesi magistrale in Economia e Finanza sulla tail dependence.

## Fonte ufficiale

La copia ufficiale del progetto è questo repository GitHub. RStudio, Codex ed eventualmente Overleaf devono sincronizzarsi con questa cronologia, evitando copie parallele non versionate.

## Struttura

- `literature/papers/`: PDF forniti dal relatore, conservati solo localmente;
- `notes/Appunti_tail_dependence.md`: appunti progressivi sui paper e sull'argomento;
- `R/`: funzioni R riutilizzabili;
- `analysis/`: script di analisi eseguiti in ordine numerico;
- `data/`: istruzioni e metadati relativi ai dati;
- `output/`: figure, tabelle e risultati riproducibili;
- `thesis/`: sorgenti LaTeX della tesi;
- `WORKFLOW.md`: regole operative, Git e convenzioni dei commit.

I paper PDF sono conservati localmente e intenzionalmente esclusi da GitHub.

## Apertura in RStudio

Aprire `Tesi.Rproj`. Il progetto parte con un ambiente R vuoto: `.RData` e `.Rhistory` non vengono utilizzati come fonte dei risultati. Ogni risultato importante deve poter essere rigenerato dagli script.

## Riproducibilità R

Quando inizierà la parte empirica, inizializzare `renv` dalla console di RStudio:

```r
install.packages("renv")
renv::init()
```

Dopo avere installato o aggiornato pacchetti e verificato che l'analisi funzioni:

```r
renv::snapshot()
```

Il file `renv.lock` dovrà essere incluso nel commit insieme al codice che richiede le nuove dipendenze.

## Tesi LaTeX

Il file principale è `thesis/Tesi_tail_dependence.tex`. Il primo capitolo, `thesis/chapters/01-quadro-della-letteratura.tex`, introduce e collega i tre testi assegnati. La struttura sarà adattata al template ufficiale dell'università non appena disponibile.
