# Stima della Generalized Tail Concentration Function

Codice R della tesi magistrale di **Massimo Gatti**, *Stima della Generalized Tail Concentration Function: un confronto Monte Carlo tra PIT teorici e ranghi campionari* (Economia e finanza, Università degli Studi di Milano-Bicocca; relatore: Prof. Lucio De Capitani).

Il lavoro confronta due modi di stimare la dipendenza di coda a partire dalle stesse osservazioni simulate: applicare le funzioni di ripartizione marginali vere (PIT teorici) oppure costruire pseudo-osservazioni dai ranghi campionari. Il Capitolo 2 descrive il metodo e il Capitolo 3 presenta i risultati.

## Quale script eseguire

I tre script sono autonomi: ciascuno genera i dati, stima i parametri nelle due code e salva i risultati della propria copula. Non richiedono dati esterni né un ordine di esecuzione.

| Script | Copula | Risultati |
| --- | --- | --- |
| [MC_Gaussiana.R](MC_Gaussiana.R) | Gaussiana | 4 CSV; per il coefficiente `c` il confronto riguarda media e varianza sulla finestra, senza bias o MSE rispetto a un valore asintotico positivo. |
| [MC_Clayton.R](MC_Clayton.R) | Clayton | 3 CSV con riepiloghi per parametro e confronto appaiato tra PIT teorici e ranghi. |
| [MC_Gumbel.R](MC_Gumbel.R) | Gumbel | 3 CSV con gli stessi tipi di riepilogo e confronto utilizzati per Clayton. |

## Esecuzione

Servono R e i pacchetti `dplyr`, `tidyr`, `ggplot2` e `scales`. Per installarli, dalla console R:

```r
install.packages(c("dplyr", "tidyr", "ggplot2", "scales"))
```

Dopo aver scaricato la repository, aprire un terminale nella cartella che contiene i tre file ed eseguire lo script desiderato, ad esempio:

```sh
Rscript --vanilla MC_Clayton.R
```

In alternativa, in RStudio aprire una nuova sessione, impostare la stessa cartella come directory di lavoro ed eseguire lo script dall'inizio alla fine. Ogni script inizia con `rm(list = ls())`. Per l'esportazione dei PDF viene usato il dispositivo grafico Cairo di R (verificabile con `capabilities("cairo")`).

## Disegno Monte Carlo

- Tre valori del tau di Kendall: 0,20, 1/3 e 0,50. Per ogni copula e valore di tau si eseguono 1.000 repliche. Ogni script fissa il seme a 44 all'avvio.
- In ogni replica si genera un campione di 100.000 coppie e se ne prendono le prime 5.000, 10.000, 20.000, 50.000 e 100.000. I ranghi vengono ricalcolati per ciascun sottocampione.
- Per ogni coda si usano 20 soglie equidistanti su una finestra dinamica con esponente 0,25. La regressione log–log delle frequenze empiriche fornisce le stime di α e c.
- I confronti tra le due basi utilizzano, per ciascun parametro, soltanto le repliche valide per entrambe. Nel codice e in alcuni nomi di file, `IID` indica i PIT teorici.

## Dove trovare gli output

Ogni script crea la propria cartella, relativa alla directory di lavoro:

- `output/Gaussiana/Monte Carlo Gaussiana/`
- `output/Clayton/Monte Carlo Clayton/`
- `output/Gumbel/Monte Carlo Gumbel/`

Oltre ai CSV, ogni script salva 14 grafici sia in PNG sia in PDF vettoriale e un PDF multipagina che li raccoglie. Una nuova esecuzione nella stessa cartella sovrascrive i file omonimi. Gli output si generano localmente e non sono inclusi nella repository.
