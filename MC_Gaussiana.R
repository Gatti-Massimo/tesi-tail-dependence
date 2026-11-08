# M.C. - Copula Gaussiana



# 0. Ambiente ------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(44)



# 1. Parametri e modello globale -----------------------------------------------


# 1.1 Disegno Monte Carlo

# Fisso tre livelli di dipendenza misurati attraverso il tau di Kendall
valori_tau <- c(0.20, 1/3, 0.50)

# Per ogni replica genero una sola volta il campione di dimensione massima,
# i campioni più piccoli saranno ottenuti selezionando le prime n osservazioni
numerosita_campioni <- c(5000L, 10000L, 20000L, 50000L, 100000L)

# Numero di repliche previste per ciascun valore di tau
numero_repliche <- 1000L


# 1.2 Parametro della copula e benchmark teorici

# La copula gaussiana ha gli stessi benchmark nelle due code
# rho = sin(pi * tau / 2), kappa = 2 / (1 + rho) e alpha = 1 / kappa
benchmark_teorici <- data.frame(tau = valori_tau) %>%
  mutate(
    rho = sin(pi * tau / 2),
    kappa_vero = 2 / (1 + rho),
    alpha_vero = 1 / kappa_vero
  )

# Per c_hat confronto le basi sulla stessa finestra, senza un c vero costante
# Il fattore lentamente variabile della coda gaussiana non è una costante positiva
# da usare come benchmark dell'intercetta stimata nella finestra.


# 1.3 Marginali normali

# Fisso media e varianza delle due variabili
mu <- c(X = 0.01, Y = 0.03)

varianze <- c(X = 0.001, Y = 0.003)

sigma <- sqrt(varianze)


# 1.4 Funzione per simulare il modello globale gaussiano

# La funzione restituisce le variabili osservabili X e Y
simula_modello_globale <- function(n, rho, mu, sigma) {

  # Genero due variabili normali standard indipendenti
  z_x <- rnorm(n)
  errore_indipendente <- rnorm(n)

  # Costruisco una seconda normale standard correlata con z_x
  # La correlazione tra z_x e z_y è uguale a rho
  z_y <- rho * z_x + sqrt(1 - rho^2) * errore_indipendente

  # Applico le marginali normali
  # La dipendenza deriva dalla copula mentre mu e sigma determinano le marginali
  data.frame(
    X = mu["X"] + sigma["X"] * z_x,
    Y = mu["Y"] + sigma["Y"] * z_y
  )
}



# 2. Costruzione delle basi informative ----------------------------------------


# 2.1 Trasformazione con marginali vere e metodo dei ranghi

# Dagli stessi valori osservabili X e Y costruisco due basi informative
#
# La base IID utilizza le distribuzioni marginali vere
#
# La base Ranghi utilizza invece le distribuzioni marginali empiriche
# e riproduce la situazione nella quale le marginali vere non sono conosciute
costruisci_basi_informative <- function(dati_xy, mu, sigma) {

  # Ricavo la numerosità del campione
  n <- nrow(dati_xy)

  dati_xy %>%
    mutate(

      # Applico le CDF marginali vere alle osservazioni X e Y
      # Le coppie sono indipendenti tra righe ma dipendenti al loro interno
      U_iid = pnorm(X, mean = mu["X"], sd = sigma["X"]),
      V_iid = pnorm(Y, mean = mu["Y"], sd = sigma["Y"]),

      # Trasformo X e Y in pseudo-osservazioni attraverso i ranghi
      U_ranghi = rank(X) / (n + 1),
      V_ranghi = rank(Y) / (n + 1)
    )
}



# 3. Finestre dinamiche e soglie -----------------------------------------------


# 3.1 Parametri della finestra

# epsilon determina la velocità con cui le soglie si avvicinano al bordo
epsilon <- 0.25

# Uso la stessa finestra nelle due code, data la simmetria della copula
# Con n = 10.000 ed epsilon = 0,25, x è compreso tra 0,01 e 0,08
# Nella upper x rappresenta la distanza dal bordo: conto i superamenti di 1 - x
h_meno <- 0.10
h_piu <- 0.80

# Mantengo lo stesso numero di soglie per tutte le numerosità
numero_soglie <- 20L


# 3.2 Funzione per costruire le soglie dinamiche

costruisci_soglie <- function(n, epsilon, h_meno, h_piu, numero_soglie) {

  # Calcolo gli estremi della finestra relativi alla numerosità n
  # All'aumentare di n entrambi gli estremi si avvicinano a zero
  x_min <- h_meno * n^(-epsilon)
  x_max <- h_piu * n^(-epsilon)

  # seq() distribuisce lo stesso numero di soglie nella nuova finestra
  # La distanza tra soglie consecutive diminuisce quando la finestra si restringe
  soglie_x <- seq(from = x_min, to = x_max, length.out = numero_soglie)

  return(soglie_x)
}



# 4. Probabilità empiriche di coda ---------------------------------------------


# 4.1 Funzione per calcolare le probabilità empiriche di coda

calcola_probabilita_coda <- function(U, V, soglie_x) {

  numerosita <- length(U)

  # Per ogni soglia conto le coppie nella coda inferiore
  conteggi_lower <- vapply(
    soglie_x,
    function(x) sum(U <= x & V <= x),
    integer(1)
  )

  # Nella coda superiore x rappresenta la distanza dal bordo
  conteggi_upper <- vapply(
    soglie_x,
    function(x) sum(U > 1 - x & V > 1 - x),
    integer(1)
  )

  # Restituisco soltanto le soglie e le probabilità utilizzate nelle regressioni
  data.frame(
    soglia_x = soglie_x,
    probabilita_lower = conteggi_lower / numerosita,
    probabilita_upper = conteggi_upper / numerosita
  )
}



# 5. Regressione log-log e stima dei parametri ---------------------------------


# 5.1 Funzione per stimare i parametri GTCF

stima_parametri_gtcf <- function(soglie, probabilita) {

  stima_non_valida <- data.frame(
    beta_hat = NA_real_,
    alpha_hat = NA_real_,
    c_hat = NA_real_,
    stima_valida = FALSE
  )

  # Escludo probabilità nulle e finestre con probabilità tutte uguali
  if (any(probabilita == 0) || length(unique(probabilita)) < 2L) {
    return(stima_non_valida)
  }

  # Stimo log(P(x)) = log(c) + beta * log(x)
  modello <- lm(log(probabilita) ~ log(soglie))

  intercetta <- coef(modello)[[1]]
  beta <- coef(modello)[[2]]

  # La pendenza deve essere positiva e finita
  if (!is.finite(beta) || beta <= 0) {
    return(stima_non_valida)
  }

  # Impongo soltanto il limite superiore alpha <= 1
  # Eventuali valori inferiori a 1/2 rimangono visibili
  alpha <- min(1, 1 / beta)

  # Stimo c dall'intercetta della regressione originale
  c_hat <- exp(intercetta)

  data.frame(
    beta_hat = beta,
    alpha_hat = alpha,
    c_hat = c_hat,
    stima_valida = TRUE
  )
}



# 6. Applicazione della procedura a un campione --------------------------------


# 6.1 Funzione per stimare le due code nelle due basi informative

analizza_campione <- function(dati_xy, mu, sigma, epsilon,
    h_meno, h_piu, numero_soglie) {

  # Ricavo la numerosità
  n <- nrow(dati_xy)

  # Costruisco PIT teorici e Ranghi dagli stessi dati X e Y
  #
  # da: [ Sezione 2 ]
  basi <- costruisci_basi_informative(dati_xy, mu, sigma)

  # Costruisco le soglie dinamiche comuni alle due code
  #
  # da: [ Sezione 3 ]
  soglie_x <- costruisci_soglie(n, epsilon, h_meno, h_piu, numero_soglie)

  # Calcolo le probabilità di coda per entrambe le basi
  #
  # da: [ Sezione 4 ]
  probabilita_iid <- calcola_probabilita_coda(basi$U_iid, basi$V_iid, soglie_x)
  probabilita_ranghi <- calcola_probabilita_coda(basi$U_ranghi, basi$V_ranghi, soglie_x)

  # Stimo alpha e c nelle quattro combinazioni di base e coda
  #
  # da: [ Sezione 5 ]
  stima_iid_lower <- stima_parametri_gtcf(probabilita_iid$soglia_x, probabilita_iid$probabilita_lower)
  stima_iid_upper <- stima_parametri_gtcf(probabilita_iid$soglia_x, probabilita_iid$probabilita_upper)
  stima_ranghi_lower <- stima_parametri_gtcf(probabilita_ranghi$soglia_x, probabilita_ranghi$probabilita_lower)
  stima_ranghi_upper <- stima_parametri_gtcf(probabilita_ranghi$soglia_x, probabilita_ranghi$probabilita_upper)

  # Unisco le quattro stime indicando numerosità, base informativa e coda
  risultati <- bind_rows(
    data.frame(n = n, base_informativa = "IID", coda = "lower", stima_iid_lower),
    data.frame(n = n, base_informativa = "IID", coda = "upper", stima_iid_upper),
    data.frame(n = n, base_informativa = "Ranghi", coda = "lower", stima_ranghi_lower),
    data.frame(n = n, base_informativa = "Ranghi", coda = "upper", stima_ranghi_upper)
  )

  return(risultati)
}



# 7. Funzione per una replica Monte Carlo --------------------------------------


# 7.1 Simulazione e analisi di una replica per un valore di tau

simula_replica_monte_carlo <- function(tau, replica, numerosita_campioni, mu, sigma, epsilon,
    h_meno, h_piu, numero_soglie) {

  # Ricavo da tau il parametro necessario
  rho <- sin(pi * tau / 2)

  # Genero una sola volta il campione di numerosità massima
  #
  # da: [ Sezione 1 ]
  dati_massimi <- simula_modello_globale(
    n = max(numerosita_campioni),
    rho = rho,
    mu = mu,
    sigma = sigma
  )

  # Analizzo i sottocampioni annidati per tutte le numerosità previste
  risultati_n <- lapply(numerosita_campioni, function(n) {

    # Prendo le prime n osservazioni così i campioni rimangono annidati
    dati_n <- dati_massimi[seq_len(n), ]

    # Analizzo il sottocampione e ricalcolo i ranghi sulle sue n osservazioni
    #
    # da: [ Sezione 6 ]
    stime_n <- analizza_campione(dati_xy = dati_n, mu = mu, sigma = sigma, epsilon = epsilon,
      h_meno = h_meno, h_piu = h_piu, numero_soglie = numero_soglie
    )

    # Aggiungo gli identificativi della replica e dello scenario
    data.frame(
      replica = replica,
      tau = tau,
      stime_n
    )
  })

  # Unisco i risultati delle diverse numerosità
  return(bind_rows(risultati_n))
}



# 8. Esecuzione del Monte Carlo ------------------------------------------------


# 8.1 Costruzione degli scenari Monte Carlo

# Creo una riga per ogni combinazione di tau e numero della replica -> 3x1000
scenari_monte_carlo <- expand.grid(
  tau = valori_tau,
  replica = seq_len(numero_repliche),
  KEEP.OUT.ATTRS = FALSE # Informazioni di struttura
)

# Ordino gli scenari prima per tau e poi per numero della replica
scenari_monte_carlo <- scenari_monte_carlo %>% arrange(tau, replica)


# 8.2 Esecuzione di tutti gli scenari

# Ogni chiamata produce le stime per le cinque numerosità, le due basi e le due code
risultati_monte_carlo <- bind_rows(

  lapply(seq_len(nrow(scenari_monte_carlo)),

         function(indice_scenario) {

            # Seleziono la singola riga con il valore di tau e il numero della replica da eseguire
            scenario <- scenari_monte_carlo[indice_scenario, ]

            # Genero il campione massimo e analizzo tutte le numerosità annidate
            #
            # da: [ Sezione 7 ]
            risultato_scenario <- simula_replica_monte_carlo(
              tau = scenario$tau,
              replica = scenario$replica,
              numerosita_campioni = numerosita_campioni,
              mu = mu,
              sigma = sigma,
              epsilon = epsilon,
              h_meno = h_meno,
              h_piu = h_piu,
              numero_soglie = numero_soglie
            )

            # Visualizzo lo status di avanzamento (QoL)
            if (scenario$replica %% 100L == 0L || scenario$replica == numero_repliche) {
              cat(
                "tau =", scenario$tau,
                "| repliche completate:", scenario$replica, "/", numero_repliche, "\n"
              )
            }

            # Restituisco una riga per ogni combinazione di numerosità, base e coda
            # Con cinque numerosità ottengo 20 righe per ciascuna replica e tau
            return(risultato_scenario)

          })
)



# 9. Confronto con il benchmark teorico ----------------------------------------


# 9.1 Aggiunta dei benchmark teorici e calcolo degli errori

# Confronto alpha_hat con il valore vero del proprio tau, uguale nelle due code
risultati_monte_carlo <- risultati_monte_carlo %>%

  left_join(benchmark_teorici, by = "tau") %>%

  mutate(
    errore_alpha = alpha_hat - alpha_vero,

    # Distinguo le stime utilizzabili per alpha e per c
    valida_alpha = coalesce(stima_valida, FALSE) & is.finite(alpha_hat),
    valida_c = coalesce(stima_valida, FALSE) & is.finite(c_hat) & c_hat > 0
  )


# 9.2 Bias, varianza e MSE di alpha

# Una riga per tau, numerosità, base informativa e coda
# Calcolo le performance sulle stime valide della singola base
sintesi_monte_carlo <- risultati_monte_carlo %>%

  group_by(tau, n, base_informativa, coda, alpha_vero) %>%

  summarise(
    numero_repliche = n(),
    numero_stime_valide = sum(valida_alpha),
    bias_alpha = mean(errore_alpha[valida_alpha]),
    varianza_alpha = mean((alpha_hat[valida_alpha] - mean(alpha_hat[valida_alpha]))^2),
    mse_alpha = mean(errore_alpha[valida_alpha]^2),
    .groups = "drop"
  ) %>%

  # Se mancano stime valide, le performance restano mancanti
  mutate(across(c(bias_alpha, varianza_alpha, mse_alpha),
                ~ replace(.x, is.nan(.x), NA_real_))) %>%

  arrange(tau, n, base_informativa, coda)


# 9.3 Media e varianza di c_hat

# Descrivo il livello stimato nella finestra senza confrontarlo con un c vero
sintesi_gtcf <- risultati_monte_carlo %>%

  group_by(tau, n, base_informativa, coda) %>%

  summarise(
    numero_repliche = n(),
    numero_stime_valide = sum(valida_c),
    media_c = mean(c_hat[valida_c]),
    varianza_c = mean((c_hat[valida_c] - mean(c_hat[valida_c]))^2),
    .groups = "drop"
  ) %>%

  mutate(across(c(media_c, varianza_c),
                ~ replace(.x, is.nan(.x), NA_real_))) %>%

  arrange(tau, n, base_informativa, coda)


# 10. Confronto tra basi informative -------------------------------------------


# 10.1 Affiancamento delle stime ottenute nella stessa replica

# Porto sulla stessa riga le stime relative alla stessa replica, tau, n e coda
repliche_affiancate <- risultati_monte_carlo %>%

  select(tau, n, replica, coda, alpha_vero, base_informativa,
         alpha_hat, c_hat, valida_alpha, valida_c) %>%

  pivot_wider(
    id_cols = c(tau, n, replica, coda, alpha_vero),
    names_from = base_informativa,
    values_from = c(alpha_hat, c_hat, valida_alpha, valida_c),
    names_glue = "{.value}_{base_informativa}"
  )


# 10.2 Costruzione della tabella finale di confronto

# Per ciascun parametro confronto le basi sulle medesime coppie valide
# Calcolo gli errori di alpha prima della selezione, rispetto al benchmark di ogni riga
confronto_basi <- repliche_affiancate %>%

  mutate(
    coppia_valida = coalesce(valida_alpha_IID, FALSE) &
      coalesce(valida_alpha_Ranghi, FALSE),

    coppia_valida_c = coalesce(valida_c_IID, FALSE) &
      coalesce(valida_c_Ranghi, FALSE),

    errore_alpha_IID = alpha_hat_IID - alpha_vero,
    errore_alpha_Ranghi = alpha_hat_Ranghi - alpha_vero
  ) %>%

  group_by(tau, n, coda, alpha_vero) %>%

  summarise(
    numero_repliche = n(),

    # 10.2.1 Performance appaiate di alpha

    numero_coppie_valide = sum(coppia_valida),
    bias_alpha_IID = mean(errore_alpha_IID[coppia_valida]),
    bias_alpha_Ranghi = mean(errore_alpha_Ranghi[coppia_valida]),

    # Il denominatore della varianza è il numero di coppie valide
    varianza_alpha_IID = mean((alpha_hat_IID[coppia_valida] - mean(alpha_hat_IID[coppia_valida]))^2),
    varianza_alpha_Ranghi = mean((alpha_hat_Ranghi[coppia_valida] - mean(alpha_hat_Ranghi[coppia_valida]))^2),

    mse_alpha_IID = mean(errore_alpha_IID[coppia_valida]^2),
    mse_alpha_Ranghi = mean(errore_alpha_Ranghi[coppia_valida]^2),

    # 10.2.2 Media e varianza appaiate di c_hat

    numero_coppie_valide_c = sum(coppia_valida_c),
    media_c_IID = mean(c_hat_IID[coppia_valida_c]),
    media_c_Ranghi = mean(c_hat_Ranghi[coppia_valida_c]),

    varianza_c_IID = mean((c_hat_IID[coppia_valida_c] - mean(c_hat_IID[coppia_valida_c]))^2),
    varianza_c_Ranghi = mean((c_hat_Ranghi[coppia_valida_c] - mean(c_hat_Ranghi[coppia_valida_c]))^2),

    .groups = "drop"
  ) %>%

  mutate(

    # 10.2.3 Misure comparative di alpha

    # Valori positivi favoriscono i Ranghi
    # Il vantaggio nel bias è normalizzato rispetto ad alpha vero
    vantaggio_bias_percentuale = 100 * (abs(bias_alpha_IID) - abs(bias_alpha_Ranghi)) / alpha_vero,

    riduzione_varianza = if_else(
      varianza_alpha_IID > 0,
      100 * (1 - varianza_alpha_Ranghi / varianza_alpha_IID),
      NA_real_
    ),

    riduzione_mse = if_else(
      mse_alpha_IID > 0,
      100 * (1 - mse_alpha_Ranghi / mse_alpha_IID),
      NA_real_
    ),

    # 10.2.4 Misure comparative di c_hat

    # Il segno indica quale base ha la media maggiore, non quale è più accurata
    differenza_media_c = media_c_Ranghi - media_c_IID,

    riduzione_varianza_c = if_else(
      varianza_c_IID > 0,
      100 * (1 - varianza_c_Ranghi / varianza_c_IID),
      NA_real_
    ),

    # Converto soltanto le performance non disponibili
    across(matches("^(bias_|media_|differenza_|varianza_|mse_|vantaggio_|riduzione_)"),
           ~ replace(.x, is.nan(.x), NA_real_))
  ) %>%

  arrange(tau, n, coda)



# 11. SISTEMA GRAFICO E GRAFICI MONTE CARLO ====================================


# 11.1 Palette e codifiche grafiche comuni

col_lower <- "#67A9CF"
col_upper <- "#EF8A62"
col_riferimento <- "#777777"
col_griglia <- "#E6E6E6"
col_testo <- "#222222"

colori_code <- c(
  "Coda inferiore" = col_lower,
  "Coda superiore" = col_upper
)

forme_code <- c(
  "Coda inferiore" = 16,
  "Coda superiore" = 17
)

tipi_linea_code <- c(
  "Coda inferiore" = "solid",
  "Coda superiore" = "dashed"
)

tipi_linea_basi <- c(
  "PIT teorici" = "solid",
  "Ranghi" = "dashed"
)

forme_basi <- c(
  "PIT teorici" = 16,
  "Ranghi" = 17
)


# 11.2 Tema comune per le figure della tesi

# Titoli, sottotitoli e note vengono riportati nelle didascalie LaTeX.
# Nei grafici rimangono assi, intestazioni dei pannelli e legende.
tema_tesi <- function(base_size = 11, griglia_x = FALSE) {
  
  theme_minimal(
    base_size = base_size,
    base_family = "sans"
  ) +
    theme(
      plot.background = element_rect(
        fill = "white",
        color = NA
      ),
      panel.background = element_rect(
        fill = "white",
        color = NA
      ),
      
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.grid.major.y = element_line(
        color = col_griglia,
        linewidth = 0.30
      ),
      
      panel.grid.major.x = if (griglia_x) {
        element_line(
          color = "#EFEFEF",
          linewidth = 0.25
        )
      } else {
        element_blank()
      },
      
      axis.line = element_line(
        color = col_riferimento,
        linewidth = 0.30
      ),
      axis.ticks = element_line(
        color = col_riferimento,
        linewidth = 0.30
      ),
      axis.ticks.length = grid::unit(2, "pt"),
      
      axis.text = element_text(
        size = 9,
        color = col_testo
      ),
      axis.text.x = element_text(
        size = 8.5,
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      
      axis.title = element_text(
        size = base_size * 0.95,
        color = col_testo
      ),
      axis.title.x = element_text(
        margin = margin(t = 6)
      ),
      axis.title.y = element_text(
        margin = margin(r = 6)
      ),
      
      strip.background = element_rect(
        fill = "#F4F4F4",
        color = NA
      ),
      strip.text = element_text(
        size = 9,
        face = "bold",
        color = col_testo,
        margin = margin(5, 3, 5, 3)
      ),
      
      panel.spacing.x = grid::unit(0.9, "lines"),
      panel.spacing.y = grid::unit(0.9, "lines"),
      
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(
        size = base_size * 0.90,
        color = col_testo
      ),
      legend.key.width = grid::unit(1.8, "lines"),
      legend.key.height = grid::unit(0.9, "lines"),
      legend.spacing.x = grid::unit(0.4, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(4, "pt"),
      
      plot.margin = margin(6, 8, 6, 6)
    )
}


# 11.3 Etichette comuni e preparazione dei dati grafici

etichette_n <- format(
  numerosita_campioni,
  big.mark = ".",
  decimal.mark = ",",
  scientific = FALSE,
  trim = TRUE
)

# Relazione specifica della copula gaussiana
etichette_tau <- paste0(
  "τ = ",
  formatC(
    valori_tau,
    format = "f",
    digits = 3,
    decimal.mark = ","
  ),
  " | ρ = ",
  formatC(
    sin(pi * valori_tau / 2),
    format = "f",
    digits = 3,
    decimal.mark = ","
  )
)


# 11.3.1 Preparazione delle distribuzioni di alpha

riferimenti_distribuzione <- risultati_monte_carlo %>%
  distinct(tau, coda, alpha_vero) %>%
  mutate(
    pannello_etichetta = paste0(
      "τ = ",
      formatC(tau, format = "f", digits = 3, decimal.mark = ","),
      " | α vero = ",
      formatC(alpha_vero, format = "f", digits = 3, decimal.mark = ",")
    )
  )

livelli_pannelli <- riferimenti_distribuzione %>%
  arrange(
    match(coda, c("lower", "upper")),
    match(tau, valori_tau)
  ) %>%
  pull(pannello_etichetta) %>%
  unique()

dati_grafico_distribuzione <- risultati_monte_carlo %>%
  filter(valida_alpha) %>%
  left_join(
    riferimenti_distribuzione %>%
      select(tau, coda, pannello_etichetta),
    by = c("tau", "coda")
  ) %>%
  mutate(
    numerosita_etichetta = factor(
      n,
      levels = numerosita_campioni,
      labels = etichette_n
    ),
    pannello_etichetta = factor(
      pannello_etichetta,
      levels = livelli_pannelli
    )
  )

# Limiti comuni a PIT e Ranghi per la stessa coda
limiti_scale_code <- dati_grafico_distribuzione %>%
  group_by(coda) %>%
  summarise(
    minimo = min(c(alpha_hat, alpha_vero), na.rm = TRUE),
    massimo = max(c(alpha_hat, alpha_vero), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    margine = pmax(0.02, 0.04 * (massimo - minimo)),
    limite_inferiore = minimo - margine,
    limite_superiore = massimo + margine
  )


# 11.3.2 Preparazione delle distribuzioni di c

# Per c non viene tracciato un benchmark teorico
# Le intestazioni riportano tau e rho
dati_grafico_distribuzione_c <- risultati_monte_carlo %>%
  filter(valida_c) %>%
  mutate(
    numerosita_etichetta = factor(
      n,
      levels = numerosita_campioni,
      labels = etichette_n
    ),
    pannello_etichetta = factor(
      tau,
      levels = valori_tau,
      labels = etichette_tau
    )
  )

limiti_scale_code_c <- dati_grafico_distribuzione_c %>%
  group_by(coda) %>%
  summarise(
    minimo = min(c_hat),
    massimo = max(c_hat),
    .groups = "drop"
  ) %>%
  mutate(
    margine = pmax(0.02, 0.04 * (massimo - minimo)),
    limite_inferiore = minimo - margine,
    limite_superiore = massimo + margine
  )


# 11.3.3 Preparazione delle misure appaiate

# Alpha: bias, varianza e MSE
# c: media e varianza
# I valori provengono dal confronto sulle repliche valide comuni
prepara_performance <- function(parametro) {
  
  confronto_basi %>%
    select(
      tau, n, coda,
      matches(
        paste0("^(media|bias|varianza|mse)_", parametro, "_(IID|Ranghi)$")
      )
    ) %>%
    pivot_longer(
      cols = -c(tau, n, coda),
      names_to = c("misura", "base_informativa"),
      names_pattern = paste0(
        "^(media|bias|varianza|mse)_", parametro, "_(IID|Ranghi)$"
      ),
      values_to = "valore"
    ) %>%
    mutate(
      numerosita_etichetta = factor(
        match(n, numerosita_campioni),
        levels = seq_along(numerosita_campioni),
        labels = etichette_n
      ),
      tau_etichetta = factor(
        tau,
        levels = valori_tau,
        labels = etichette_tau
      ),
      base_etichetta = factor(
        base_informativa,
        levels = c("IID", "Ranghi"),
        labels = c("PIT teorici", "Ranghi")
      ),
      misura = factor(
        misura,
        levels = c("media", "bias", "varianza", "mse"),
        labels = c("Media", "Bias", "Varianza", "MSE")
      )
    )
}

dati_grafico_performance_alpha <- prepara_performance("alpha")
dati_grafico_performance_c <- prepara_performance("c")


# 11.3.4 Preparazione del confronto relativo tra le basi

# Per alpha confronto bias assoluto, varianza e MSE
prepara_vantaggio <- function() {
  
  confronto_basi %>%
    transmute(
      tau, n, coda,
      bias_assoluto = vantaggio_bias_percentuale,
      varianza = riduzione_varianza,
      mse = riduzione_mse
    ) %>%
    pivot_longer(
      cols = c(bias_assoluto, varianza, mse),
      names_to = "misura",
      values_to = "vantaggio_percentuale"
    ) %>%
    mutate(
      numerosita_etichetta = factor(
        match(n, numerosita_campioni),
        levels = seq_along(numerosita_campioni),
        labels = etichette_n
      ),
      tau_etichetta = factor(
        tau,
        levels = valori_tau,
        labels = etichette_tau
      ),
      coda_etichetta = factor(
        coda,
        levels = c("lower", "upper"),
        labels = c("Coda inferiore", "Coda superiore")
      ),
      misura = factor(
        misura,
        levels = c("bias_assoluto", "varianza", "mse"),
        labels = c("Bias assoluto", "Varianza", "MSE")
      )
    )
}

dati_grafico_vantaggio_alpha <- prepara_vantaggio()

# Per c confronto soltanto la varianza
# La minore dispersione non dimostra una maggiore accuratezza
dati_grafico_vantaggio_c <- confronto_basi %>%
  transmute(
    numerosita_etichetta = factor(
      match(n, numerosita_campioni),
      levels = seq_along(numerosita_campioni),
      labels = etichette_n
    ),
    tau_etichetta = factor(
      tau,
      levels = valori_tau,
      labels = etichette_tau
    ),
    coda_etichetta = factor(
      coda,
      levels = c("lower", "upper"),
      labels = c("Coda inferiore", "Coda superiore")
    ),
    misura = factor("Varianza"),
    vantaggio_percentuale = riduzione_varianza_c
  )


# 11.4 Grafici della distribuzione delle stime di alpha

# Boxplot sulle stime valide della singola base
# La linea grigia tratteggiata indica alpha vero
crea_grafico_distribuzione_alpha <- function(base_scelta, coda_scelta) {
  
  dati_grafico <- dati_grafico_distribuzione %>%
    filter(
      base_informativa == base_scelta,
      coda == coda_scelta
    ) %>%
    droplevels()
  
  riferimenti_grafico <- dati_grafico %>%
    distinct(pannello_etichetta, alpha_vero)
  
  limiti <- limiti_scale_code %>%
    filter(coda == coda_scelta)
  
  colore_coda <- if (coda_scelta == "lower") col_lower else col_upper
  
  ggplot(
    dati_grafico,
    aes(x = numerosita_etichetta, y = alpha_hat)
  ) +
    geom_hline(
      data = riferimenti_grafico,
      aes(yintercept = alpha_vero),
      color = col_riferimento,
      linewidth = 0.50,
      linetype = "dashed"
    ) +
    geom_boxplot(
      width = 0.52,
      fill = colore_coda,
      color = "#333333",
      alpha = 0.94,
      linewidth = 0.40,
      outlier.shape = 16,
      outlier.color = "#666666",
      outlier.alpha = 0.38,
      outlier.size = 0.65
    ) +
    facet_wrap(
      vars(pannello_etichetta),
      nrow = 1,
      scales = "fixed"
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(
      limits = c(limiti$limite_inferiore, limiti$limite_superiore),
      breaks = scales::breaks_pretty(n = 6),
      labels = scales::label_number(
        accuracy = 0.01,
        decimal.mark = ","
      ),
      expand = expansion(mult = 0)
    ) +
    labs(
      x = "Numerosità campionaria",
      y = expression(hat(alpha))
    ) +
    tema_tesi(griglia_x = FALSE) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9)
    )
}

grafico_distribuzione_alpha_iid_lower <-
  crea_grafico_distribuzione_alpha("IID", "lower")

grafico_distribuzione_alpha_iid_upper <-
  crea_grafico_distribuzione_alpha("IID", "upper")

grafico_distribuzione_alpha_ranghi_lower <-
  crea_grafico_distribuzione_alpha("Ranghi", "lower")

grafico_distribuzione_alpha_ranghi_upper <-
  crea_grafico_distribuzione_alpha("Ranghi", "upper")


# 11.5 Grafici della distribuzione delle stime di c

# c descrive l'approssimazione sulla finestra di stima
# Non viene aggiunta una linea di benchmark teorico
crea_grafico_distribuzione_c <- function(base_scelta, coda_scelta) {
  
  dati_grafico <- dati_grafico_distribuzione_c %>%
    filter(
      base_informativa == base_scelta,
      coda == coda_scelta
    ) %>%
    droplevels()
  
  limiti <- limiti_scale_code_c %>%
    filter(coda == coda_scelta)
  
  colore_coda <- if (coda_scelta == "lower") col_lower else col_upper
  
  ggplot(
    dati_grafico,
    aes(x = numerosita_etichetta, y = c_hat)
  ) +
    geom_boxplot(
      width = 0.52,
      fill = colore_coda,
      color = "#333333",
      alpha = 0.94,
      linewidth = 0.40,
      outlier.shape = 16,
      outlier.color = "#666666",
      outlier.alpha = 0.38,
      outlier.size = 0.65
    ) +
    facet_wrap(
      vars(pannello_etichetta),
      nrow = 1,
      scales = "fixed"
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(
      limits = c(limiti$limite_inferiore, limiti$limite_superiore),
      breaks = scales::breaks_pretty(n = 6),
      labels = scales::label_number(
        accuracy = 0.01,
        decimal.mark = ","
      ),
      expand = expansion(mult = 0)
    ) +
    labs(
      x = "Numerosità campionaria",
      y = expression(hat(c))
    ) +
    tema_tesi(griglia_x = FALSE) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9)
    )
}

grafico_distribuzione_c_iid_lower <-
  crea_grafico_distribuzione_c("IID", "lower")

grafico_distribuzione_c_iid_upper <-
  crea_grafico_distribuzione_c("IID", "upper")

grafico_distribuzione_c_ranghi_lower <-
  crea_grafico_distribuzione_c("Ranghi", "lower")

grafico_distribuzione_c_ranghi_upper <-
  crea_grafico_distribuzione_c("Ranghi", "upper")


# 11.5.1 Scala logaritmica facoltativa per i boxplot di c

# Trasformo le coordinate dopo il calcolo delle statistiche del boxplot
# PIT e Ranghi della stessa coda mantengono limiti identici
applica_scala_log_c <- function(grafico, coda_scelta) {
  
  dati_coda <- dati_grafico_distribuzione_c %>%
    filter(coda == coda_scelta)
  
  limiti <- range(dati_coda$c_hat)
  limiti <- limiti * c(0.9, 1.1)
  
  grafico +
    scale_y_continuous(
      breaks = scales::breaks_log(n = 6),
      labels = scales::label_number(
        decimal.mark = ",",
        big.mark = " "
      ),
      expand = expansion(mult = 0)
    ) +
    coord_transform(
      y = "log10",
      ylim = limiti,
      expand = TRUE
    ) +
    labs(y = expression(hat(c) ~ "(scala logaritmica)"))
}

# Mantengo il logaritmo attivo in entrambe le code, come nel file originale
# Per disattivarlo in una coda, commentare le due assegnazioni corrispondenti
# e rieseguire l'intera sezione 11
grafico_distribuzione_c_iid_lower <-
  applica_scala_log_c(grafico_distribuzione_c_iid_lower, "lower")

grafico_distribuzione_c_ranghi_lower <-
  applica_scala_log_c(grafico_distribuzione_c_ranghi_lower, "lower")

grafico_distribuzione_c_iid_upper <-
  applica_scala_log_c(grafico_distribuzione_c_iid_upper, "upper")

grafico_distribuzione_c_ranghi_upper <-
  applica_scala_log_c(grafico_distribuzione_c_ranghi_upper, "upper")


# 11.6 Grafici delle misure appaiate di alpha e c

# Alpha: bias, varianza e MSE.
# c: media e varianza, senza bias o MSE rispetto a un benchmark.
# Ogni riga identifica la misura rappresentata:
# non aggiungo un titolo generico all'asse Y.
# Le scale verticali sono separate tra le misure
# e comuni ai tre valori di tau della stessa riga.
crea_grafico_performance <- function(
    dati,
    coda_scelta,
    parametro
) {
  
  colore_coda <- if (coda_scelta == "lower") {
    col_lower
  } else {
    col_upper
  }
  
  ggplot(
    filter(dati, coda == coda_scelta),
    aes(
      x = numerosita_etichetta,
      y = valore,
      linetype = base_etichetta,
      shape = base_etichetta,
      group = base_etichetta
    )
  ) +
    geom_hline(
      yintercept = 0,
      color = col_riferimento,
      linewidth = 0.40,
      linetype = "dotted"
    ) +
    geom_line(
      color = colore_coda,
      linewidth = 0.55
    ) +
    geom_point(
      color = colore_coda,
      size = 1.70,
      stroke = 0.50
    ) +
    facet_grid(
      rows = vars(misura),
      cols = vars(tau_etichetta),
      scales = "free_y"
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 4),
      labels = scales::label_number(
        decimal.mark = ","
      )
    ) +
    scale_linetype_manual(
      values = tipi_linea_basi
    ) +
    scale_shape_manual(
      values = forme_basi
    ) +
    labs(
      x = "Numerosità campionaria",
      y = NULL,
      linetype = "Base informativa",
      shape = "Base informativa"
    ) +
    tema_tesi(griglia_x = TRUE)
}


# 11.6.1 Performance di alpha

grafico_performance_alpha_lower <-
  crea_grafico_performance(
    dati_grafico_performance_alpha,
    "lower",
    "alpha"
  )

grafico_performance_alpha_upper <-
  crea_grafico_performance(
    dati_grafico_performance_alpha,
    "upper",
    "alpha"
  )


# 11.6.2 Media e varianza di c

grafico_performance_c_lower <-
  crea_grafico_performance(
    dati_grafico_performance_c,
    "lower",
    "c"
  )

grafico_performance_c_upper <-
  crea_grafico_performance(
    dati_grafico_performance_c,
    "upper",
    "c"
  )


# 11.7 Grafici del confronto relativo tra Ranghi e PIT teorici

# Per alpha confronto bias assoluto, varianza e MSE.
# Per c rappresento soltanto la riduzione della varianza.
# Valori positivi favoriscono i Ranghi.
# Per c il vantaggio riguarda la stabilità tra repliche,
# non dimostra una maggiore accuratezza.
crea_grafico_vantaggio <- function(dati, parametro) {
  
  etichetta_y <- if (parametro == "alpha") {
    "Vantaggio dei Ranghi (%)"
  } else {
    "Riduzione della varianza (%)"
  }
  
  ggplot(
    dati,
    aes(
      x = numerosita_etichetta,
      y = vantaggio_percentuale,
      color = coda_etichetta,
      shape = coda_etichetta,
      linetype = coda_etichetta,
      group = coda_etichetta
    )
  ) +
    geom_hline(
      yintercept = 0,
      color = col_riferimento,
      linewidth = 0.45,
      linetype = "dotted"
    ) +
    geom_line(
      linewidth = 0.55
    ) +
    geom_point(
      size = 1.75,
      stroke = 0.50
    ) +
    facet_grid(
      rows = vars(misura),
      cols = vars(tau_etichetta),
      scales = "free_y"
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.6)
    ) +
    scale_y_continuous(
      breaks = scales::breaks_pretty(n = 4),
      labels = scales::label_number(
        accuracy = 0.1,
        decimal.mark = ",",
        suffix = "%"
      ),
      expand = expansion(mult = c(0.07, 0.10))
    ) +
    scale_color_manual(
      values = colori_code
    ) +
    scale_shape_manual(
      values = forme_code
    ) +
    scale_linetype_manual(
      values = tipi_linea_code
    ) +
    labs(
      x = "Numerosità campionaria",
      y = etichetta_y,
      color = "Coda",
      shape = "Coda",
      linetype = "Coda"
    ) +
    tema_tesi(griglia_x = TRUE)
}


# 11.7.1 Confronto relativo per alpha

grafico_vantaggio_ranghi <-
  crea_grafico_vantaggio(
    dati_grafico_vantaggio_alpha,
    "alpha"
  )


# 11.7.2 Riduzione della varianza di c

grafico_vantaggio_ranghi_c <-
  crea_grafico_vantaggio(
    dati_grafico_vantaggio_c,
    "c"
  )



# 12. SALVATAGGIO DEGLI OUTPUT =================================================


# 12.1 Cartella e tabelle di output

cartella_output <- "output/Gaussiana/Monte Carlo Gaussiana"

dir.create(
  cartella_output,
  recursive = TRUE,
  showWarnings = FALSE
)

# Le tabelle interne rimangono numeriche e invariate
tabella_output_alpha <- sintesi_monte_carlo %>%
  mutate(
    epsilon = .env$epsilon,
    base_informativa = recode(base_informativa, "IID" = "PIT"),
    .before = 1
  )

tabella_output_c <- sintesi_gtcf %>%
  mutate(
    epsilon = .env$epsilon,
    base_informativa = recode(base_informativa, "IID" = "PIT"),
    .before = 1
  )

# Alpha e c vengono esportati in confronti separati
tabella_output_confronto <- confronto_basi %>%
  transmute(
    epsilon = .env$epsilon,
    parametro = "alpha",
    tau, n, coda,
    valore_vero = alpha_vero,
    numero_repliche,
    coppie_valide = numero_coppie_valide,
    bias_PIT = bias_alpha_IID,
    bias_Ranghi = bias_alpha_Ranghi,
    varianza_PIT = varianza_alpha_IID,
    varianza_Ranghi = varianza_alpha_Ranghi,
    MSE_PIT = mse_alpha_IID,
    MSE_Ranghi = mse_alpha_Ranghi,
    vantaggio_bias_pct_su_vero = vantaggio_bias_percentuale,
    riduzione_varianza_percentuale = riduzione_varianza,
    riduzione_MSE_percentuale = riduzione_mse
  )

tabella_output_confronto_c <- confronto_basi %>%
  transmute(
    epsilon = .env$epsilon,
    parametro = "c",
    tau, n, coda,
    numero_repliche,
    coppie_valide = numero_coppie_valide_c,
    media_PIT = media_c_IID,
    media_Ranghi = media_c_Ranghi,
    differenza_media_Ranghi_PIT = differenza_media_c,
    varianza_PIT = varianza_c_IID,
    varianza_Ranghi = varianza_c_Ranghi,
    riduzione_varianza_percentuale = riduzione_varianza_c
  )


# 12.2 Salvataggio dei CSV

write.csv(
  tabella_output_alpha,
  file.path(cartella_output, "MonteCarlo_Sintesi_Gaussiana.csv"),
  row.names = FALSE
)

write.csv(
  tabella_output_c,
  file.path(cartella_output, "MonteCarlo_Sintesi_GTCF_Gaussiana.csv"),
  row.names = FALSE
)

write.csv(
  tabella_output_confronto,
  file.path(cartella_output, "MonteCarlo_Confronto_Basi_Gaussiana.csv"),
  row.names = FALSE
)

write.csv(
  tabella_output_confronto_c,
  file.path(cartella_output, "MonteCarlo_Confronto_Basi_GTCF_Gaussiana.csv"),
  row.names = FALSE
)


# 12.3 Raccolta dei grafici

# L'ordine determina anche le pagine del PDF multipagina
lista_grafici <- list(
  Distribuzione_Alpha_IID_Lower = grafico_distribuzione_alpha_iid_lower,
  Distribuzione_Alpha_IID_Upper = grafico_distribuzione_alpha_iid_upper,
  Distribuzione_Alpha_Ranghi_Lower = grafico_distribuzione_alpha_ranghi_lower,
  Distribuzione_Alpha_Ranghi_Upper = grafico_distribuzione_alpha_ranghi_upper,
  
  Distribuzione_GTCF_IID_Lower = grafico_distribuzione_c_iid_lower,
  Distribuzione_GTCF_IID_Upper = grafico_distribuzione_c_iid_upper,
  Distribuzione_GTCF_Ranghi_Lower = grafico_distribuzione_c_ranghi_lower,
  Distribuzione_GTCF_Ranghi_Upper = grafico_distribuzione_c_ranghi_upper,
  
  Performance_Alpha_Lower = grafico_performance_alpha_lower,
  Performance_Alpha_Upper = grafico_performance_alpha_upper,
  Performance_GTCF_Lower = grafico_performance_c_lower,
  Performance_GTCF_Upper = grafico_performance_c_upper,
  
  Vantaggio_Ranghi = grafico_vantaggio_ranghi,
  Vantaggio_Ranghi_GTCF = grafico_vantaggio_ranghi_c
)


# 12.4 Salvataggio dei PNG e dei PDF singoli

# Larghezza comune di 18 cm
# L'altezza dipende dal numero di righe di pannelli:
# - boxplot: 9 cm
# - confronto relativo di c, una riga: 8 cm
# - media e varianza di c, due righe: 12 cm
# - performance e confronto di alpha, tre righe: 15 cm
dimensioni_grafico <- function(nome) {
  
  if (startsWith(nome, "Distribuzione_")) {
    
    c(larghezza = 18, altezza = 9)
    
  } else if (nome == "Vantaggio_Ranghi_GTCF") {
    
    c(larghezza = 18, altezza = 8)
    
  } else if (startsWith(nome, "Performance_GTCF")) {
    
    c(larghezza = 18, altezza = 12)
    
  } else {
    
    c(larghezza = 18, altezza = 15)
  }
}

for (nome in names(lista_grafici)) {
  
  dimensioni <- dimensioni_grafico(nome)
  
  nome_file <- paste0(
    "MonteCarlo_",
    nome,
    "_Gaussiana"
  )
  
  ggsave(
    filename = file.path(
      cartella_output,
      paste0(nome_file, ".png")
    ),
    plot = lista_grafici[[nome]],
    width = unname(dimensioni["larghezza"]),
    height = unname(dimensioni["altezza"]),
    units = "cm",
    dpi = 320,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      cartella_output,
      paste0(nome_file, ".pdf")
    ),
    plot = lista_grafici[[nome]],
    device = grDevices::cairo_pdf,
    width = unname(dimensioni["larghezza"]),
    height = unname(dimensioni["altezza"]),
    units = "cm",
    bg = "white"
  )
}


# 12.5 PDF multipagina di consultazione

# Il PDF complessivo usa pagine uniformi da 18 x 15 cm
# Per la tesi si utilizzano i PDF singoli, con le dimensioni dedicate
salva_pdf <- function() {
  
  grDevices::cairo_pdf(
    filename = file.path(
      cartella_output,
      "Grafici_MonteCarlo_Gaussiana_Tesi.pdf"
    ),
    width = 18 / 2.54,
    height = 15 / 2.54,
    onefile = TRUE
  )
  
  on.exit(grDevices::dev.off(), add = TRUE)
  
  for (grafico in lista_grafici) {
    print(grafico)
  }
}

salva_pdf()


# 12.6 Riepilogo in Console

# Formatto soltanto le copie destinate alla stampa, mantenendo numerici i CSV
formatta_decimali <- function(x, cifre = 4) {
  formatC(x, format = "f", digits = cifre, decimal.mark = ",")
}

formatta_n <- function(x) {
  format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE, trim = TRUE)
}


# Tabelle dei risultati: una riga per tau e n, con basi e code affiancate
prepara_riepilogo_console <- function(tabella, colonna) {

  tabella %>%
    transmute(
      tau, n,
      configurazione = paste(base_informativa, coda, sep = "_"),
      valore = .data[[colonna]]
    ) %>%
    pivot_wider(names_from = configurazione, values_from = valore) %>%
    rename(PIT_lower = IID_lower, PIT_upper = IID_upper) %>%
    select(tau, n, PIT_lower, Ranghi_lower, PIT_upper, Ranghi_upper) %>%
    arrange(match(tau, valori_tau), n) %>%
    mutate(
      tau = formatta_decimali(tau, 3),
      n = formatta_n(n),
      across(c(PIT_lower, Ranghi_lower, PIT_upper, Ranghi_upper), formatta_decimali)
    )
}

riepilogo_alpha_console <- prepara_riepilogo_console(sintesi_monte_carlo, "bias_alpha")
riepilogo_media_c_console <- prepara_riepilogo_console(sintesi_gtcf, "media_c")


# Confronti appaiati: conto le numerosità nelle quali i Ranghi migliorano i PIT
# Per il bias confronto i valori assoluti
# Le riduzioni mediane sono calcolate sulle numerosità, non sulle repliche
prepara_confronto_console <- function() {

  confronto_basi %>%
    transmute(
      tau, coda,
      bias = vantaggio_bias_percentuale,
      varianza = riduzione_varianza,
      mse = riduzione_mse
    ) %>%
    group_by(tau, coda) %>%
    summarise(
      bias_vinti = paste0(sum(is.finite(bias) & bias > 0), "/", sum(is.finite(bias))),
      var_vinte = paste0(sum(is.finite(varianza) & varianza > 0), "/", sum(is.finite(varianza))),
      rid_var_med = median(varianza[is.finite(varianza)]),
      mse_vinti = paste0(sum(is.finite(mse) & mse > 0), "/", sum(is.finite(mse))),
      rid_mse_med = median(mse[is.finite(mse)]),
      .groups = "drop"
    ) %>%
    arrange(match(tau, valori_tau), match(coda, c("lower", "upper"))) %>%
    mutate(
      tau = formatta_decimali(tau, 3),
      across(c(rid_var_med, rid_mse_med), ~ paste0(formatta_decimali(.x, 1), "%"))
    )
}

riepilogo_confronto_alpha_console <- prepara_confronto_console()

# Per c_hat riassumo soltanto la riduzione della varianza
riepilogo_confronto_c_console <- confronto_basi %>%
  group_by(tau, coda) %>%
  summarise(
    var_vinte = paste0(sum(is.finite(riduzione_varianza_c) & riduzione_varianza_c > 0),
                       "/", sum(is.finite(riduzione_varianza_c))),
    rid_var_med = median(riduzione_varianza_c[is.finite(riduzione_varianza_c)]),
    .groups = "drop"
  ) %>%
  arrange(match(tau, valori_tau), match(coda, c("lower", "upper"))) %>%
  mutate(
    tau = formatta_decimali(tau, 3),
    rid_var_med = paste0(formatta_decimali(rid_var_med, 1), "%")
  )


# Stampa dei risultati
cat("\n\nRISULTATI PRINCIPALI\n")

cat("\nBias della media Monte Carlo di alpha_hat rispetto ad alpha vero\n\n")
print(as.data.frame(riepilogo_alpha_console), row.names = FALSE)

cat("\nMedia Monte Carlo di c_hat nella finestra di stima\n\n")
print(as.data.frame(riepilogo_media_c_console), row.names = FALSE)

cat(
  "\nCONFRONTI APPAIATI RANGHI–PIT\n",
  "Vittorie = numerosità con bias assoluto, varianza o MSE inferiori nei Ranghi / numerosità confrontabili.\n",
  "Riduzioni percentuali mediane sulle numerosità; valori positivi favoriscono i Ranghi.\n",
  sep = ""
)

cat("\nConfronto appaiato: alpha\n\n")
print(as.data.frame(riepilogo_confronto_alpha_console), row.names = FALSE)

cat("\nConfronto appaiato: varianza di c_hat (non accuratezza)\n\n")
print(as.data.frame(riepilogo_confronto_c_console), row.names = FALSE)

cat("\nEsecuzione completata\n")
