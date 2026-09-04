# Monte Carlo - Copula Gaussiana


# 1. Ambiente, parametri e benchmark -------------------------------------------

# 1.1 Ambiente
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(44)

# Per tibble in console
#
# options(
#   pillar.sigfig = 6,
#   pillar.neg = TRUE
# )


# 1.2 Disegno Monte Carlo
# Stesse Kendall tau utilizzate per Clayton e Gumbel
tabella_parametri <- tibble::tibble(
  tau_kendall = c(0.20, 1 / 3, 0.50)
) |>
  mutate(
    rho = sin(pi * tau_kendall / 2)
  )

numerosita_campioni <- c(5000L, 10000L, 20000L, 50000L, 100000L)

numero_repliche <- 1000L


# 1.3 Finestra fissa di stima
# La copula gaussiana è simmetrica: lower e upper utilizzano
# la stessa finestra già definita nell'analisi singola.
x_min_stima <- 0.01
x_max_stima <- 0.10
passo_x <- 0.005

soglie_stima <- seq(
  from = x_min_stima,
  to = x_max_stima,
  by = passo_x
)


# 1.4 Benchmark teorici Gaussiani
# Ottenuti applicando la stessa regressione log-log utilizzata sulle probabilità empiriche.
#
# Lower e upper condividono gli stessi benchmark per la simmetria della copula.
#
# Per ogni Rho vengono conservati:
# - livello GTCF teorico della finestra
# - beta teorico della finestra;
# - alpha teorico della finestra;
# - benchmark asintotici di alpha e del livello GTCF.
lista_benchmark <- lapply(

  seq_len(nrow(tabella_parametri)),

  function(indice_parametro) {

    tau_corrente <- tabella_parametri$tau_kendall[[indice_parametro]]
    rho_corrente <- tabella_parametri$rho[[indice_parametro]]


    # 1.4.1 ---

    # Probabilità teoriche nella finestra [0.01, 0.10]
    # Lower e upper coincidono teoricamente per simmetria
    probabilita_teoriche <-
      copBasic::NORMcop(
        u = soglie_stima,
        v = soglie_stima,
        para = rho_corrente
      )


    # 1.4.2 ---

    # Regressione teorica log-log sulla finestra utilizzata poi anche nelle repliche Monte Carlo
    modello_teorico <- lm(log(probabilita_teoriche) ~ log(soglie_stima))

    intercetta_teorica_finestra <- coef(modello_teorico)[[1]]          # a
    livello_gtcf_teorico_finestra <- exp(intercetta_teorica_finestra) # c
    beta_teorico_finestra <- coef(modello_teorico)[[2]]                # beta
    alpha_teorico_finestra <- 1 / beta_teorico_finestra                # alfa


    # 1.4.3 ---

    # Benchmark asintotici Gaussiani
    beta_asintotico <- 2 / (1 + rho_corrente)
    alpha_asintotico <- (1 + rho_corrente) / 2

    # Alla normalizzazione asintotica la GTCF tende a zero
    livello_gtcf_asintotico <- 0


    # 1.4.4 ---

    # Riga di informazioni per quel Rho
    tibble::tibble(
      rho = rho_corrente,
      tau_kendall = tau_corrente,
      livello_gtcf_teorico_finestra = livello_gtcf_teorico_finestra,
      beta_teorico_finestra = beta_teorico_finestra,
      alpha_teorico_finestra = alpha_teorico_finestra,
      beta_asintotico = beta_asintotico,
      alpha_asintotico = alpha_asintotico,
      livello_gtcf_asintotico = livello_gtcf_asintotico
    )
  }
)

# Unisci tutto
tabella_benchmark <- bind_rows(lista_benchmark)

# Mostra in console
cat("\nBenchmark teorici\n")
print(tabella_benchmark, n = Inf, width = Inf)



# 2. Funzioni essenziali -------------------------------------------------------

# 2.1 Trasformazione in pseudo-osservazioni
# Converte le osservazioni nei ranghi uniformi compresi tra 0 e 1.
crea_pseudo_osservazioni <- function(x) {
  rank(x, ties.method = "average") / (length(x) + 1)
}


# 2.2 Probabilità empiriche di coda
# Calcola conteggi e probabilità congiunte lower e upper per ogni soglia.
calcola_probabilita_coda <- function(U, V, soglie) {

  numerosita <- length(U)

  conteggi_lower <- vapply(
    soglie,
    function(x) sum(U <= x & V <= x),
    integer(1)
  )

  conteggi_upper <- vapply(
    soglie,
    function(x) sum(U > 1 - x & V > 1 - x),
    integer(1)
  )

  data.frame(
    soglia_x = soglie,
    conteggi_lower = conteggi_lower,
    probabilita_lower = conteggi_lower / numerosita,
    conteggi_upper = conteggi_upper,
    probabilita_upper = conteggi_upper / numerosita
  )
}


# 2.3 Stima dei parametri GTCF
# Regressione log-log con controlli numerici
stima_parametri_gtcf <- function(soglie, probabilita) {

  # ! CONTROLLO
  # La regressione non è stimabile se una probabilità è nulla o non finita.
  if (any(!is.finite(probabilita)) || any(probabilita <= 0)) {
    return(
      list(
        intercetta = NA_real_,
        livello_gtcf = NA_real_,
        beta = NA_real_,
        alpha = NA_real_
      )
    )
  }

  # Stima
  modello <- lm(log(probabilita) ~ log(soglie))
  intercetta <- coef(modello)[[1]]
  beta <- coef(modello)[[2]]

  # ! CONTROLLO
  # Beta deve esistere e rappresentare una velocità di decadimento positiva
  if (!is.finite(beta) || beta <= 0) {
    return(
      list(
        intercetta = NA_real_,
        livello_gtcf = NA_real_,
        beta = NA_real_,
        alpha = NA_real_
      )
    )
  }

  # Output
  list(
    intercetta = intercetta,
    livello_gtcf = exp(intercetta),
    beta = beta,
    alpha = 1 / beta
  )
}


# 2.4 Singola replica Monte Carlo
# Simula un campione Gaussiano e applica separatamente alle due code
# la stessa procedura utilizzata nell'Analisi singola.
esegui_replica_gaussiana <- function(rho,
                                     numerosita,
                                     soglie,
                                     replica) {


  # 2.4.1 ---

  # INIZIO BLOCCO GENERAZIONE GAUSSIANA

  # Prima normale standard
  X <- rnorm(numerosita)

  # Componente indipendente
  errore_indipendente <- rnorm(numerosita)

  # Seconda normale standard correlata con X
  Y <- rho * X + sqrt(1 - rho^2) * errore_indipendente

  # Trasformazione in pseudo-osservazioni
  U <- crea_pseudo_osservazioni(X)
  V <- crea_pseudo_osservazioni(Y)

  # FINE BLOCCO GENERAZIONE GAUSSIANA


  # 2.4.2 ---

  # Calcolo conteggi e probabilità empiriche lower e upper
  probabilita_coda <- calcola_probabilita_coda(
    U = U,
    V = V,
    soglie = soglie
  )

  # Stime lower e upper sulla stessa finestra fissa
  stima_lower <- stima_parametri_gtcf(
    soglie = soglie,
    probabilita = probabilita_coda$probabilita_lower
  )

  stima_upper <- stima_parametri_gtcf(
    soglie = soglie,
    probabilita = probabilita_coda$probabilita_upper
  )


  # 2.4.3 ---

  # Una riga per ciascuna coda
  tibble::tibble(

    rho = rho,
    numerosita = as.integer(numerosita),
    replica = as.integer(replica),

    coda = c(
      "Lower",
      "Upper"
    ),

    intercetta_loglog = c(
      stima_lower$intercetta,
      stima_upper$intercetta
    ),

    livello_gtcf = c(
      stima_lower$livello_gtcf,
      stima_upper$livello_gtcf
    ),

    beta = c(
      stima_lower$beta,
      stima_upper$beta
    ),

    alpha = c(
      stima_lower$alpha,
      stima_upper$alpha
    ),

    eventi_xmin = c(
      probabilita_coda$conteggi_lower[1],
      probabilita_coda$conteggi_upper[1]
    ),

    valida = c(
      all(
        is.finite(
          c(
            stima_lower$intercetta,
            stima_lower$livello_gtcf,
            stima_lower$beta,
            stima_lower$alpha
          )
        )
      ),
      all(
        is.finite(
          c(
            stima_upper$intercetta,
            stima_upper$livello_gtcf,
            stima_upper$beta,
            stima_upper$alpha
          )
        )
      )
    )
  )
}



# 3. Esecuzione del Monte Carlo ------------------------------------------------

# Inizializzo
numero_risultati <- nrow(tabella_parametri) * length(numerosita_campioni) * numero_repliche
lista_risultati_montecarlo <- vector(mode = "list", length = numero_risultati)
indice_risultato <- 1L # Prima posizione

# Avvio
cat("\nEsecuzione Monte Carlo\n")

for (indice_parametro in seq_len(nrow(tabella_parametri))) {

  tau_corrente <- tabella_parametri$tau_kendall[[indice_parametro]]
  rho_corrente <- tabella_parametri$rho[[indice_parametro]]

  for (numerosita_corrente in numerosita_campioni) {

    for (replica_corrente in seq_len(numero_repliche)) {

      # Salvo i risultati della singola esecuzione
      lista_risultati_montecarlo[[indice_risultato]] <-
        esegui_replica_gaussiana(
          rho = rho_corrente,
          numerosita = numerosita_corrente,
          soglie = soglie_stima,
          replica = replica_corrente
        )

      # Incremento
      indice_risultato <- indice_risultato + 1L

      # Stampa lo status
      if (replica_corrente %% 100L == 0L || replica_corrente == numero_repliche) {
        cat(
          "tau =", round(tau_corrente, 6),
          "| rho =", round(rho_corrente, 6),
          "| n =", numerosita_corrente,
          "| repliche completate:",
          replica_corrente, "/", numero_repliche,
          "\n"
        )
      }
    }
  }
}

# Unisco in tabella
risultati_grezzi <- bind_rows(lista_risultati_montecarlo)



# 4. Preparazione dei risultati ------------------------------------------------

# 4.1 Benchmark in formato long (Verticalizzo)
# Duplico i benchmark comuni in una riga lower e una riga upper.
benchmark_long <- bind_rows(

  # Tabella Lower
  tabella_benchmark |>
    transmute(
      rho,
      tau_kendall,
      coda = "Lower",
      livello_gtcf_teorico_finestra,
      beta_teorico_finestra,
      alpha_teorico_finestra,
      alpha_asintotico,
      livello_gtcf_asintotico
    ),

  # Tabella Upper
  tabella_benchmark |>
    transmute(
      rho,
      tau_kendall,
      coda = "Upper",
      livello_gtcf_teorico_finestra,
      beta_teorico_finestra,
      alpha_teorico_finestra,
      alpha_asintotico,
      livello_gtcf_asintotico
    )
)


# 4.2 Collegamento tra stime empiriche e benchmark teorici
# Per ogni riga di risultati_grezzi, cerco in benchmark_long lo stesso rho e la stessa coda,
# e aggiungo i benchmark corrispondenti
risultati_completi <- risultati_grezzi |>
  left_join(
    benchmark_long,
    by = c(
      "rho",
      "coda"
    )
  )



# 5. Sintesi dei risultati -----------------------------------------------------

# 5.1 Indicatori sintetici
# Statistiche Monte Carlo per tau, rho, numerosità e coda

sintesi_risultati <- risultati_completi |>
  group_by(
    tau_kendall,
    rho,
    numerosita,
    coda
  ) |>
  summarise(
    percentuale_stime_numericamente_valide = 100 * mean(valida),

    alpha_teorico_finestra = first(alpha_teorico_finestra),
    alpha_asintotico = first(alpha_asintotico),

    media_alpha = mean(alpha, na.rm = TRUE),
    mediana_alpha = median(alpha, na.rm = TRUE),

    # Bias alfa
    bias_alpha =
      mean(
        alpha - alpha_teorico_finestra,
        na.rm = TRUE
      ),

    deviazione_standard_alpha = sd(alpha, na.rm = TRUE),
    iqr_alpha = IQR(alpha, na.rm = TRUE),

    # RMSE alfa
    rmse_alpha =
      sqrt(
        mean(
          (alpha - alpha_teorico_finestra)^2,
          na.rm = TRUE
        )
      ),

    livello_gtcf_teorico_finestra = first(livello_gtcf_teorico_finestra),
    livello_gtcf_asintotico = first(livello_gtcf_asintotico),

    media_livello_gtcf = mean(livello_gtcf, na.rm = TRUE),
    mediana_livello_gtcf = median(livello_gtcf, na.rm = TRUE),

    # Bias GTCF
    bias_livello_gtcf =
      mean(
        livello_gtcf -
          livello_gtcf_teorico_finestra,
        na.rm = TRUE
      ),

    deviazione_standard_livello_gtcf = sd(livello_gtcf, na.rm = TRUE),
    iqr_livello_gtcf = IQR(livello_gtcf, na.rm = TRUE),

    # RMSE GTCF
    rmse_livello_gtcf =
      sqrt(
        mean(
          (
            livello_gtcf -
              livello_gtcf_teorico_finestra
          )^2,
          na.rm = TRUE
        )
      ),

    .groups = "drop"
  ) |>
  mutate(

    bias_relativo_alpha_percentuale = 100 * bias_alpha / alpha_teorico_finestra,
    rmse_relativo_alpha_percentuale = 100 * rmse_alpha / alpha_teorico_finestra,

    bias_relativo_livello_gtcf_percentuale = 100 * bias_livello_gtcf / livello_gtcf_teorico_finestra,
    rmse_relativo_livello_gtcf_percentuale = 100 * rmse_livello_gtcf / livello_gtcf_teorico_finestra

  ) |>
  arrange(
    tau_kendall,
    numerosita,
    coda
  )


# 5.2 Risultati principali
tabella_risultati_principali <- sintesi_risultati |>
  select(
    tau_kendall,
    rho,
    numerosita,
    coda,
    percentuale_stime_numericamente_valide,

    alpha_teorico_finestra,
    media_alpha,
    mediana_alpha,
    bias_alpha,
    deviazione_standard_alpha,
    rmse_alpha,

    livello_gtcf_teorico_finestra,
    media_livello_gtcf,
    mediana_livello_gtcf,
    bias_livello_gtcf,
    deviazione_standard_livello_gtcf,
    rmse_livello_gtcf
  )

# Console
cat("\nRisultati principali\n")
print(tabella_risultati_principali, n = Inf, width = Inf)



# 6. Preparazione dei dati grafici ---------------------------------------------

# 6.1 Livelli ed etichette comuni

livelli_coda <- c(
  "Coda inferiore",
  "Coda superiore"
)

livelli_parametro <- sprintf(
  "τ = %.3f | ρ = %.3f",
  tabella_parametri$tau_kendall,
  tabella_parametri$rho
)

livelli_numerosita <- format(
  numerosita_campioni,
  big.mark = ".",
  decimal.mark = ",",
  scientific = FALSE,
  trim = TRUE
)


# 6.2 Dati per le distribuzioni di alpha e livello GTCF

dati_distribuzioni <- risultati_completi |>
  mutate(
    coda_grafico = recode(
      coda,
      Lower = "Coda inferiore",
      Upper = "Coda superiore"
    ),

    coda_grafico = factor(
      coda_grafico,
      levels = livelli_coda
    ),

    parametro_pannello = factor(
      sprintf(
        "τ = %.3f | ρ = %.3f",
        tau_kendall,
        rho
      ),
      levels = livelli_parametro
    ),

    numerosita_pannello = factor(
      format(
        numerosita,
        big.mark = ".",
        decimal.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      levels = livelli_numerosita
    )
  )


# 6.3 Benchmark per i grafici delle distribuzioni

benchmark_distribuzioni <- benchmark_long |>
  mutate(
    coda_grafico = recode(
      coda,
      Lower = "Coda inferiore",
      Upper = "Coda superiore"
    ),

    coda_grafico = factor(
      coda_grafico,
      levels = livelli_coda
    ),

    parametro_pannello = factor(
      sprintf(
        "τ = %.3f | ρ = %.3f",
        tau_kendall,
        rho
      ),
      levels = livelli_parametro
    )
  )


# 6.4 Dati per bias e RMSE relativi di alpha

dati_performance_alpha <- sintesi_risultati |>
  transmute(
    parametro_pannello = factor(
      sprintf(
        "τ = %.3f | ρ = %.3f",
        tau_kendall,
        rho
      ),
      levels = livelli_parametro
    ),

    numerosita_pannello = factor(
      format(
        numerosita,
        big.mark = ".",
        decimal.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      levels = livelli_numerosita
    ),

    coda_grafico = factor(
      recode(
        coda,
        Lower = "Coda inferiore",
        Upper = "Coda superiore"
      ),
      levels = livelli_coda
    ),

    Bias = bias_relativo_alpha_percentuale,
    RMSE = rmse_relativo_alpha_percentuale
  ) |>
  pivot_longer(
    cols = c(
      Bias,
      RMSE
    ),
    names_to = "statistica",
    values_to = "valore_percentuale"
  ) |>
  mutate(
    statistica = factor(
      statistica,
      levels = c(
        "Bias",
        "RMSE"
      )
    )
  )


# 6.5 Dati per bias e RMSE relativi del livello GTCF

dati_performance_livello_gtcf <- sintesi_risultati |>
  transmute(
    parametro_pannello = factor(
      sprintf(
        "τ = %.3f | ρ = %.3f",
        tau_kendall,
        rho
      ),
      levels = livelli_parametro
    ),

    numerosita_pannello = factor(
      format(
        numerosita,
        big.mark = ".",
        decimal.mark = ",",
        scientific = FALSE,
        trim = TRUE
      ),
      levels = livelli_numerosita
    ),

    coda_grafico = factor(
      recode(
        coda,
        Lower = "Coda inferiore",
        Upper = "Coda superiore"
      ),
      levels = livelli_coda
    ),

    Bias = bias_relativo_livello_gtcf_percentuale,
    RMSE = rmse_relativo_livello_gtcf_percentuale
  ) |>
  pivot_longer(
    cols = c(
      Bias,
      RMSE
    ),
    names_to = "statistica",
    values_to = "valore_percentuale"
  ) |>
  mutate(
    statistica = factor(
      statistica,
      levels = c(
        "Bias",
        "RMSE"
      )
    )
  )



# 7. Sistema grafico comune -----------------------------------------------------

# 7.1 Palette
# Mantiene lo stesso linguaggio visivo dell'Analisi Gaussiana
# e del Monte Carlo Clayton.

col_lower <- "#67A9CF"
col_upper <- "#EF8A62"
col_riferimento <- "#9AA1A8"
col_griglia <- "#E6E6E6"
col_testo <- "#222222"


# 7.2 Tema comune

tema_tesi <- function(base_size = 12, griglia_x = FALSE) {
  theme_minimal(
    base_size = base_size,
    base_family = "sans"
  ) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        face = "bold",
        size = base_size * 1.28,
        color = col_testo,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = base_size * 0.90,
        color = "#5A5A5A",
        margin = margin(b = 10)
      ),
      plot.caption.position = "plot",
      plot.caption = element_text(
        size = base_size * 0.76,
        color = "#6A6A6A",
        hjust = 0,
        margin = margin(t = 8)
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(
        color = col_griglia,
        linewidth = 0.35
      ),
      panel.grid.major.x = if (griglia_x) {
        element_line(
          color = col_griglia,
          linewidth = 0.35
        )
      } else {
        element_blank()
      },
      axis.line = element_line(
        color = "#9A9A9A",
        linewidth = 0.35
      ),
      axis.title = element_text(
        size = base_size * 0.95,
        color = col_testo
      ),
      axis.text = element_text(
        color = col_testo
      ),
      axis.ticks = element_line(
        color = "#9A9A9A",
        linewidth = 0.35
      ),
      axis.ticks.length = grid::unit(2.5, "pt"),
      strip.background = element_rect(
        fill = "#F4F4F4",
        color = NA
      ),
      strip.text = element_text(
        face = "bold",
        color = col_testo,
        margin = margin(6, 6, 6, 6)
      ),
      panel.spacing = grid::unit(1.4, "lines"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = grid::unit(1.8, "lines"),
      legend.spacing.x = grid::unit(0.6, "lines"),
      plot.margin = margin(10, 14, 10, 10)
    )
}



# 8. Grafici Monte Carlo --------------------------------------------------------

# 8.1 Distribuzione delle stime di alpha

benchmark_alpha_grafico <- benchmark_distribuzioni |>
  select(
    parametro_pannello,
    alpha_teorico_finestra,
    alpha_asintotico
  ) |>
  distinct() |>
  pivot_longer(
    cols = c(
      alpha_teorico_finestra,
      alpha_asintotico
    ),
    names_to = "benchmark",
    values_to = "alpha_benchmark"
  ) |>
  mutate(
    benchmark = recode(
      benchmark,
      alpha_teorico_finestra = "Teorico della finestra",
      alpha_asintotico = "Asintotico"
    ),

    benchmark = factor(
      benchmark,
      levels = c(
        "Teorico della finestra",
        "Asintotico"
      )
    )
  )


grafico_distribuzione_alpha <- ggplot(
  dati_distribuzioni,
  aes(
    x = numerosita_pannello,
    y = alpha,
    fill = coda_grafico
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.60,
    outlier.alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_hline(
    data = benchmark_alpha_grafico,
    aes(
      yintercept = alpha_benchmark,
      linetype = benchmark
    ),
    colour = col_testo,
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ parametro_pannello,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Coda inferiore" = col_lower,
      "Coda superiore" = col_upper
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Teorico della finestra" = "solid",
      "Asintotico" = "22"
    )
  ) +
  labs(
    title = "Distribuzione delle stime di α",
    subtitle = paste0(
      "Linea continua = benchmark di finestra; ",
      "linea tratteggiata = riferimento asintotico"
    ),
    x = "Numerosità campionaria",
    y = expression(hat(alpha)),
    caption = paste0(
      numero_repliche,
      " repliche per configurazione. ",
      "Finestra comune lower e upper [0,01; 0,10]. ",
      "Le stime non sono vincolate a [0,5; 1]."
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  tema_tesi()


# 8.2 Distribuzione delle stime del livello GTCF

benchmark_livello_gtcf_grafico <- benchmark_distribuzioni |>
  select(
    parametro_pannello,
    livello_gtcf_teorico_finestra
  ) |>
  distinct()


grafico_distribuzione_livello_gtcf <- ggplot(
  dati_distribuzioni,
  aes(
    x = numerosita_pannello,
    y = livello_gtcf,
    fill = coda_grafico
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.60,
    outlier.alpha = 0.45,
    na.rm = TRUE
  ) +
  geom_hline(
    data = benchmark_livello_gtcf_grafico,
    aes(
      yintercept = livello_gtcf_teorico_finestra
    ),
    colour = col_testo,
    linewidth = 0.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ parametro_pannello,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Coda inferiore" = col_lower,
      "Coda superiore" = col_upper
    )
  ) +
  scale_y_log10() +
  labs(
    title = "Distribuzione delle stime del livello GTCF",
    subtitle = paste0(
      "Scala logaritmica; linea continua = ",
      "benchmark teorico della finestra"
    ),
    x = "Numerosità campionaria",
    y = expression(hat(c)),
    caption = paste0(
      numero_repliche,
      " repliche per configurazione. ",
      "Tutte le stime, inclusi gli outlier, sono conservate. ",
      "Per la Gaussiana non esiste un livello asintotico positivo."
    )
  ) +
  guides(
    fill = guide_legend(order = 1)
  ) +
  tema_tesi()


# 8.3 Funzione comune per i grafici di performance relativa

crea_grafico_performance <- function(dati,
                                     titolo,
                                     sottotitolo) {

  ggplot(
    dati,
    aes(
      x = numerosita_pannello,
      y = valore_percentuale,
      colour = coda_grafico,
      linetype = coda_grafico,
      shape = coda_grafico,
      group = coda_grafico
    )
  ) +
    geom_hline(
      data = dati |>
        distinct(
          parametro_pannello,
          statistica
        ) |>
        filter(statistica == "Bias"),
      aes(yintercept = 0),
      inherit.aes = FALSE,
      linetype = "22",
      colour = col_riferimento,
      linewidth = 0.7
    ) +
    geom_line(
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    geom_point(
      size = 2.4,
      na.rm = TRUE
    ) +
    facet_grid(
      statistica ~ parametro_pannello,
      scales = "free_y"
    ) +
    scale_colour_manual(
      values = c(
        "Coda inferiore" = col_lower,
        "Coda superiore" = col_upper
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Coda inferiore" = "solid",
        "Coda superiore" = "22"
      )
    ) +
    scale_shape_manual(
      values = c(
        "Coda inferiore" = 16,
        "Coda superiore" = 17
      )
    ) +
    scale_y_continuous(
      labels = scales::label_number(
        suffix = "%",
        decimal.mark = ",",
        accuracy = 0.1
      )
    ) +
    labs(
      title = titolo,
      subtitle = sottotitolo,
      x = "Numerosità campionaria",
      y = "Valore relativo al benchmark di finestra",
      caption = paste0(
        "Bias e RMSE sono espressi in percentuale rispetto al benchmark ",
        "teorico della stessa finestra. I valori assoluti restano nel CSV."
      )
    ) +
    tema_tesi(
      griglia_x = TRUE
    )
}


# 8.4 Bias e RMSE relativi di alpha

grafico_performance_alpha <- crea_grafico_performance(
  dati = dati_performance_alpha,
  titolo = "Performance relativa delle stime di α",
  sottotitolo = paste0(
    "Bias e RMSE percentuali al crescere ",
    "della numerosità campionaria"
  )
)


# 8.5 Bias e RMSE relativi del livello GTCF

grafico_performance_livello_gtcf <- crea_grafico_performance(
  dati = dati_performance_livello_gtcf,
  titolo = "Performance relativa delle stime del livello GTCF",
  sottotitolo = paste0(
    "Bias e RMSE percentuali al crescere ",
    "della numerosità campionaria"
  )
)



# 9. Salvataggio degli output ---------------------------------------------------

# 9.1 CSV

tabella_output_repliche <- risultati_completi |>
  select(
    tau_kendall,
    rho,
    numerosita,
    replica,
    coda,
    intercetta_loglog,
    livello_gtcf,
    beta,
    alpha,
    eventi_xmin,
    valida
  ) |>
  arrange(
    tau_kendall,
    numerosita,
    replica,
    coda
  )

write.csv(
  tabella_output_repliche,
  "output/Gaussiana/Monte Carlo Gaussiana/MonteCarlo_Repliche_Gaussiana.csv",
  row.names = FALSE
)

write.csv(
  sintesi_risultati,
  "output/Gaussiana/Monte Carlo Gaussiana/MonteCarlo_Risultati_Gaussiana.csv",
  row.names = FALSE
)


# 9.2 Raccolta degli oggetti grafici

lista_grafici <- list(
  distribuzione_alpha = grafico_distribuzione_alpha,
  distribuzione_livello_gtcf = grafico_distribuzione_livello_gtcf,
  performance_alpha = grafico_performance_alpha,
  performance_livello_gtcf = grafico_performance_livello_gtcf
)


# 9.3 Funzione di salvataggio dei grafici

salva_grafico <- function(grafico,
                          nome_file,
                          larghezza,
                          altezza) {

  ggsave(
    filename = nome_file,
    plot = grafico,
    width = larghezza,
    height = altezza,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}


# 9.4 Specifiche dei file PNG

specifiche_grafici <- tibble::tribble(
  ~nome_grafico, ~nome_file, ~larghezza, ~altezza,

  "distribuzione_alpha",
  "output/Gaussiana/Monte Carlo Gaussiana/MonteCarlo_Distribuzione_Alpha_Gaussiana.png",
  14, 8,

  "distribuzione_livello_gtcf",
  "output/Gaussiana/Monte Carlo Gaussiana/MonteCarlo_Distribuzione_Livello_GTCF_Gaussiana.png",
  14, 8,

  "performance_alpha",
  "output/Gaussiana/Monte Carlo Gaussiana/MonteCarlo_Performance_Alpha_Gaussiana.png",
  14, 8,

  "performance_livello_gtcf",
  "output/Gaussiana/Monte Carlo Gaussiana/MonteCarlo_Performance_Livello_GTCF_Gaussiana.png",
  14, 8
)


for (i in seq_len(nrow(specifiche_grafici))) {
  salva_grafico(
    grafico = lista_grafici[[specifiche_grafici$nome_grafico[[i]]]],
    nome_file = specifiche_grafici$nome_file[[i]],
    larghezza = specifiche_grafici$larghezza[[i]],
    altezza = specifiche_grafici$altezza[[i]]
  )
}


# 9.5 PDF multipagina

grDevices::cairo_pdf(
  filename =
    "output/Gaussiana/Monte Carlo Gaussiana/Grafici_MonteCarlo_Gaussiana.pdf",
  width = 14,
  height = 9,
  onefile = TRUE
)

for (grafico in lista_grafici) {
  print(grafico)
}

dev.off()



# 10. Riepilogo finale ----------------------------------------------------------

cat(
  "\nEsecuzione completata: risultati, tabelle e grafici salvati.\n",
  "Output: output/Gaussiana/Monte Carlo Gaussiana/\n"
)
