# Simulazioni Monte Carlo per la stima della GTCF

Qui ho raccolto i tre script R che ho usato nella tesi. Per ogni campione simulato stimo i parametri α e c della *Generalized Tail Concentration Function* in due modi: prima uso le distribuzioni marginali note per ottenere i PIT teorici, poi ripeto la stima usando soltanto i ranghi delle osservazioni. Confronto i risultati separatamente per la coda inferiore e per quella superiore.

## Come leggere il codice

Ho scritto uno script per ciascuna copula: [gaussiana](MC_Gaussiana.R), [Clayton](MC_Clayton.R) e [Gumbel](MC_Gumbel.R). Ogni script genera i propri dati e si può eseguire da solo. Ho mantenuto la stessa numerazione delle sezioni in tutti e tre, così è facile ritrovare uno stesso passaggio:

| Sezioni | Contenuto |
| --- | --- |
| 1–3 | Generazione del campione, costruzione dei PIT e dei ranghi, scelta delle soglie. |
| 4–5 | Calcolo delle frequenze di coda e stima di α e c mediante regressione log–log. |
| 6–8 | Applicazione ai sottocampioni ed esecuzione delle 1.000 repliche per ciascun valore del tau di Kendall. |
| 9–10 | Riepiloghi per singola base e confronto tra PIT e ranghi sulle repliche valide per entrambi i metodi. |
| 11–12 | Costruzione dei grafici e salvataggio dei risultati. |

Nel codice, la sigla `IID` indica i PIT teorici. Per seguire il confronto discusso nel capitolo 3, il passaggio centrale è la sezione 10: è qui che affianco le stime ottenute dalle stesse osservazioni nella stessa replica.

## Come eseguire una simulazione

Dalla console di R si possono installare i pacchetti necessari:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "scales"))
```

Per riprodurre, ad esempio, i risultati di Clayton, basta aprire un terminale nella cartella dei tre script ed eseguire:

```sh
Rscript --vanilla MC_Clayton.R
```

In RStudio conviene partire da una sessione vuota: ogni script cancella gli oggetti presenti nell’ambiente di lavoro e fissa il seme a 44. Una volta impostata come directory di lavoro la cartella degli script, si può eseguire l’intero file con **Source**. Per esportare i PDF serve il supporto Cairo di R, verificabile con `capabilities("cairo")`.

## Dove trovare i risultati

Per Clayton, lo script salva tabelle e grafici in `output/Clayton/Monte Carlo Clayton/`; per le altre copule cambia il nome della famiglia nel percorso. Ogni cartella contiene i riepiloghi in CSV, 14 grafici in PNG e PDF e un unico PDF che li raccoglie tutti.

Per il confronto tra i due metodi, il file da aprire è quello con `Confronto_Basi` nel nome. Per Clayton e Gumbel contiene entrambi i parametri. Per la gaussiana ho separato i risultati: `MonteCarlo_Confronto_Basi_Gaussiana.csv` riguarda α, mentre `MonteCarlo_Confronto_Basi_GTCF_Gaussiana.csv` riguarda c. Per quest’ultimo confronto riporto media e varianza: c descrive l’approssimazione sulla finestra di stima e non ha un valore asintotico positivo rispetto al quale calcolare il bias e l’MSE.
