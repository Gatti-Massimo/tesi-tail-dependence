# Analisi - Copula Gaussiana


# 1. Ambiente e parametri -------------------------------------------------------

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


# 1.2 Parametri dell'analisi
n_osservazioni <- 10000L

tau_kendall <- 1 / 3
rho <- sin(pi * tau_kendall / 2)

x_min_stima <- 0.01
x_max_stima <- 0.10
x_max_gtcf <- 0.50
passo_x <- 0.005


# 1.3 Soglie di analisi
soglie_x <- seq(
  from = x_min_stima,
  to = x_max_gtcf,
  by = passo_x
)



# 2. Funzioni essenziali --------------------------------------------------------

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

  tibble::tibble(
    soglia_x = soglie,
    conteggi_lower = conteggi_lower,
    probabilita_lower = conteggi_lower / numerosita,
    conteggi_upper = conteggi_upper,
    probabilita_upper = conteggi_upper / numerosita
  )
}


# 2.3 Stima log-log di beta, alpha e livello della GTCF
# Stima la regressione log-log, ricava alpha = 1 / beta
# e il livello c = exp(intercetta) associato alla potenza stimata.
stima_beta_alpha <- function(soglie, probabilita) {

  modello <- lm(
    log(probabilita) ~ log(soglie)
  )

  intercetta <- coef(modello)[[1]]
  beta <- coef(modello)[[2]]

  list(
    intercetta = intercetta,
    coefficiente_potenza = exp(intercetta),
    beta = beta,
    alpha = 1 / beta,
    r_quadro_descrittivo = summary(modello)$r.squared
  )
}


# 2.4 Costruzione della GTCF
# Applica la parametrizzazione G_s,alpha(x) = P_s(x) / x^(1 / alpha).
calcola_gtcf <- function(probabilita_congiunta, soglia_x, alpha) {
  probabilita_congiunta / soglia_x^(1 / alpha)
}



# 3. Simulazione e pseudo-osservazioni -----------------------------------------

# 3.1 Simulazione delle osservazioni originali
# Simula dalla copula gaussiana e applica marginali normali standard.
dati_copula <- copBasic::simCOP(
  n = n_osservazioni,
  cop = copBasic::NORMcop,
  para = rho,
  graphics = FALSE,
  resamv01 = TRUE
)

dati_originali <- data.frame(
  X = qnorm(dati_copula$U),
  Y = qnorm(dati_copula$V)
)


# 3.2 Trasformazione in pseudo-osservazioni

U <- crea_pseudo_osservazioni(dati_originali[, "X"])
V <- crea_pseudo_osservazioni(dati_originali[, "Y"])



# 4. Probabilità di coda e stime ------------------------------------------------

# Perchè le teoriche?
# L'empirico è stimato su una finestra finita, quindi non va confrontato direttamente solo con l'alpha asintotico.
# Questo elimina l'effetto dovuto semplicemente al fatto che stiamo lavorando
# su una regione finita e non nel limite x -> 0.


# 4.1 Probabilità di coda

# 4.1A Probabilità empiriche
# Calcola conteggi e probabilità congiunte osservate
tabella_probabilita <- calcola_probabilita_coda(U = U, V = V, soglie = soglie_x)

# 4.1B Probabilità teoriche
# Ci dà, per ogni soglia (x), la probabilità congiunta teorica esatta prevista dalla copula Gaussiana
tabella_probabilita <- tabella_probabilita |>
  mutate(

    # Lower tail: P(U <= x, V <= x) = C(x, x)
    probabilita_teorica_lower = copBasic::NORMcop(u = soglia_x, v = soglia_x, para = rho),

    # Upper tail: P(U > 1-x, V > 1-x)
    probabilita_teorica_upper = 2 * soglia_x - 1 + copBasic::NORMcop(u = 1 - soglia_x, v = 1 - soglia_x, para = rho)
  )


# 4.2 Finestra finale di stima
# La Gaussiana è simmetrica, quindi lower e upper usano la stessa finestra.
tabella_probabilita_stima <- tabella_probabilita |>
  filter(
    soglia_x >= x_min_stima,
    soglia_x <= x_max_stima
  )


# 4.3 Stima log-log nella finestra

# 4.3A Stime empiriche
# Prendo dalla tabella filtrata Lower
#' @soglie(x)
#' @probabilità_empiriche(x)
stima_empirica_lower <- stima_beta_alpha(
  soglie = tabella_probabilita_stima$soglia_x,
  probabilita = tabella_probabilita_stima$probabilita_lower
)

# Uguale per la Upper
stima_empirica_upper <- stima_beta_alpha(
  soglie = tabella_probabilita_stima$soglia_x,
  probabilita = tabella_probabilita_stima$probabilita_upper
)

# Prendo alfa e beta
alpha_empirico_lower <- stima_empirica_lower$alpha
alpha_empirico_upper <- stima_empirica_upper$alpha
beta_empirico_lower <- stima_empirica_lower$beta
beta_empirico_upper <- stima_empirica_upper$beta

# 4.3B Benchmark teorici di finestra
# Stessa identica cosa, ma invece delle probabilità empiriche uso quelle teoriche.
stima_teorica_lower <- stima_beta_alpha(
  soglie = tabella_probabilita_stima$soglia_x,
  probabilita = tabella_probabilita_stima$probabilita_teorica_lower
)

# Uguale per la Upper
stima_teorica_upper <- stima_beta_alpha(
  soglie = tabella_probabilita_stima$soglia_x,
  probabilita = tabella_probabilita_stima$probabilita_teorica_upper
)

# Prendo alfa e beta
alpha_teorico_finestra_lower <- stima_teorica_lower$alpha
alpha_teorico_finestra_upper <- stima_teorica_upper$alpha
beta_teorico_finestra_lower <- stima_teorica_lower$beta
beta_teorico_finestra_upper <- stima_teorica_upper$beta


# 4.4 Benchmark asintotici
# Valori teorici di riferimento della copula Gaussiana nel limite x -> 0.

# Lower e Upper sono uguali per simmetria:
# beta = 2/(1+rho), alpha = (1+rho)/2, limite GTCF = 0.
beta_asintotico_lower <- 2 / (1 + rho)
alpha_asintotico_lower <- (1 + rho) / 2
limite_gtcf_asintotico_lower <- 0

beta_asintotico_upper <- 2 / (1 + rho)
alpha_asintotico_upper <- (1 + rho) / 2
limite_gtcf_asintotico_upper <- 0


# 4.5 Livelli della GTCF nella finestra
# Ricava c = exp(intercetta) che rappresenta il livello di riferimento
# della GTCF normalizzata con l'alpha stimato nella rispettiva finestra.

# 4.5A Livelli empirici
livello_gtcf_empirico_lower <- stima_empirica_lower$coefficiente_potenza
livello_gtcf_empirico_upper <- stima_empirica_upper$coefficiente_potenza

# 4.5B Livelli teorici
livello_gtcf_teorico_finestra_lower <- stima_teorica_lower$coefficiente_potenza
livello_gtcf_teorico_finestra_upper <- stima_teorica_upper$coefficiente_potenza



# 5. GTCF -----------------------------------------------------

# 5.1 Costruzione della GTCF

# Crea una nuova tabella contenente solo alcune colonne specificate
tabella_gtcf <- tabella_probabilita |>
  transmute(

    soglia_x,

    # 5.1A GTCF empiriche
    # Normalizzazione empirica di finestra
    gtcf_empirica_lower = calcola_gtcf(
      probabilita_congiunta = probabilita_lower,
      soglia_x = soglia_x,
      alpha = alpha_empirico_lower
    ),

    gtcf_empirica_upper = calcola_gtcf(
      probabilita_congiunta = probabilita_upper,
      soglia_x = soglia_x,
      alpha = alpha_empirico_upper
    ),

    # 5.1B GTCF teoriche
    # Normalizzazione teorica di finestra
    gtcf_teorica_finestra_lower = calcola_gtcf(
      probabilita_congiunta = probabilita_teorica_lower,
      soglia_x = soglia_x,
      alpha = alpha_teorico_finestra_lower
    ),

    gtcf_teorica_finestra_upper = calcola_gtcf(
      probabilita_congiunta = probabilita_teorica_upper,
      soglia_x = soglia_x,
      alpha = alpha_teorico_finestra_upper
    ),

    # 5.1C GTCF asintotiche
    # Normalizzazione asintotica
    gtcf_teorica_asintotica_lower = calcola_gtcf(
      probabilita_congiunta = probabilita_teorica_lower,
      soglia_x = soglia_x,
      alpha = alpha_asintotico_lower
    ),

    gtcf_teorica_asintotica_upper = calcola_gtcf(
      probabilita_congiunta = probabilita_teorica_upper,
      soglia_x = soglia_x,
      alpha = alpha_asintotico_upper
    )
  )



# 6. Tabella finale dei risultati ----------------------------------------------

# Raccoglie in un'unica tabella la configurazione dell'analisi,
# le stime empiriche e i relativi benchmark teorici.

tabella_risultati_finali <- tibble::tibble(

  # 6.1 Configurazione

  # Righe
  coda = c("Lower", "Upper"),

  # Colonne
  numerosita = n_osservazioni,
  rho = rho,
  tau_kendall = tau_kendall,

  soglia_minima = x_min_stima,
  soglia_massima = x_max_stima,

  numero_soglie = nrow(tabella_probabilita_stima),

  eventi_soglia_minima = c(
    tabella_probabilita_stima$conteggi_lower[1],
    tabella_probabilita_stima$conteggi_upper[1]
  ),


  # 6.2 Risultati empirici

  beta_empirico = c(
    beta_empirico_lower,
    beta_empirico_upper
  ),

  alpha_empirico = c(
    alpha_empirico_lower,
    alpha_empirico_upper
  ),

  livello_gtcf_empirico = c(
    livello_gtcf_empirico_lower,
    livello_gtcf_empirico_upper
  ),

  r_quadro_descrittivo = c(
    stima_empirica_lower$r_quadro_descrittivo,
    stima_empirica_upper$r_quadro_descrittivo
  ),


  # 6.3 Benchmark teorici di finestra

  beta_teorico_finestra = c(
    beta_teorico_finestra_lower,
    beta_teorico_finestra_upper
  ),

  alpha_teorico_finestra = c(
    alpha_teorico_finestra_lower,
    alpha_teorico_finestra_upper
  ),

  livello_gtcf_teorico_finestra = c(
    livello_gtcf_teorico_finestra_lower,
    livello_gtcf_teorico_finestra_upper
  ),

  # Aggiungo confronto tramite differenze numeriche
  errore_alpha = c(
    alpha_empirico_lower - alpha_teorico_finestra_lower,
    alpha_empirico_upper - alpha_teorico_finestra_upper
  ),

  errore_livello_gtcf = c(
    livello_gtcf_empirico_lower - livello_gtcf_teorico_finestra_lower,
    livello_gtcf_empirico_upper - livello_gtcf_teorico_finestra_upper
  ),


  # 6.4 Benchmark teorici asintotici

  beta_asintotico = c(
    beta_asintotico_lower,
    beta_asintotico_upper
  ),

  alpha_asintotico = c(
    alpha_asintotico_lower,
    alpha_asintotico_upper
  ),

  livello_gtcf_asintotico = c(
    limite_gtcf_asintotico_lower,
    limite_gtcf_asintotico_upper
  )
)



# 7. Preparazione dei dati per i grafici ---------------------------------------

# 7.1 Etichette comuni

livelli_coda <- c(
  "Coda inferiore",
  "Coda superiore"
)


# 7.2 Probabilità nella finestra di stima
# Prepara probabilità empiriche, teoriche e potenze stimate
# utilizzate nel grafico log-log.

dati_probabilita_lower <- tabella_probabilita_stima |>
  transmute(
    soglia_x,
    empirica = probabilita_lower,
    teorica = probabilita_teorica_lower,
    potenza =
      stima_empirica_lower$coefficiente_potenza *
      soglia_x^beta_empirico_lower
  ) |>
  pivot_longer(
    cols = c(empirica, teorica, potenza),
    names_to = "serie",
    values_to = "probabilita"
  ) |>
  mutate(
    coda = "Coda inferiore"
  )


dati_probabilita_upper <- tabella_probabilita_stima |>
  transmute(
    soglia_x,
    empirica = probabilita_upper,
    teorica = probabilita_teorica_upper,
    potenza =
      stima_empirica_upper$coefficiente_potenza *
      soglia_x^beta_empirico_upper
  ) |>
  pivot_longer(
    cols = c(empirica, teorica, potenza),
    names_to = "serie",
    values_to = "probabilita"
  ) |>
  mutate(
    coda = "Coda superiore"
  )


dati_probabilita_grafico <- bind_rows(
  dati_probabilita_lower,
  dati_probabilita_upper
) |>
  mutate(
    coda = factor(
      coda,
      levels = livelli_coda
    ),
    serie = recode(
      serie,
      empirica = "Empirica",
      teorica = "Teorica",
      potenza = "Potenza stimata"
    )
  )


# 7.3 GTCF teorica: alpha di finestra e alpha asintotico
# La copula gaussiana è simmetrica: le GTCF teoriche lower e upper
# coincidono e vengono quindi rappresentate una sola volta.

dati_gtcf_alpha_grafico <- tabella_gtcf |>
  transmute(
    soglia_x,
    finestra = gtcf_teorica_finestra_lower,
    asintotico = gtcf_teorica_asintotica_lower
  ) |>
  pivot_longer(
    cols = c(finestra, asintotico),
    names_to = "normalizzazione",
    values_to = "gtcf"
  ) |>
  mutate(
    normalizzazione = recode(
      normalizzazione,
      finestra = "Alpha di finestra",
      asintotico = "Alpha asintotico"
    ),
    normalizzazione = factor(
      normalizzazione,
      levels = c(
        "Alpha di finestra",
        "Alpha asintotico"
      )
    )
  )


# 7.4 GTCF empirica finale
# Seleziona nella finestra finale la GTCF normalizzata
# con l'alpha empirico stimato in ciascuna coda.

dati_gtcf_finale_empirica <- bind_rows(

  tabella_gtcf |>
    filter(
      soglia_x >= x_min_stima,
      soglia_x <= x_max_stima
    ) |>
    transmute(
      soglia_x,
      gtcf = gtcf_empirica_lower,
      coda = "Coda inferiore"
    ),

  tabella_gtcf |>
    filter(
      soglia_x >= x_min_stima,
      soglia_x <= x_max_stima
    ) |>
    transmute(
      soglia_x,
      gtcf = gtcf_empirica_upper,
      coda = "Coda superiore"
    )

) |>
  mutate(
    coda = factor(
      coda,
      levels = livelli_coda
    )
  )


# Livelli locali mostrati nel grafico della GTCF empirica.
# Non viene inserito un limite asintotico positivo.
dati_riferimenti_gtcf <- tibble::tibble(
  coda = factor(
    c(
      "Coda inferiore",
      "Coda inferiore",
      "Coda superiore",
      "Coda superiore"
    ),
    levels = livelli_coda
  ),

  tipo = c(
    "Livello stimato",
    "Benchmark di finestra",
    "Livello stimato",
    "Benchmark di finestra"
  ),

  livello = c(
    livello_gtcf_empirico_lower,
    livello_gtcf_teorico_finestra_lower,
    livello_gtcf_empirico_upper,
    livello_gtcf_teorico_finestra_upper
  )
)


# 7.5 Visione completa delle due code
# La lower viene rappresentata verso 0 e la upper verso 1.

dati_curva_intera <- bind_rows(

  tabella_probabilita |>
    transmute(
      posizione = soglia_x,
      probabilita = probabilita_lower,
      coda = "Coda inferiore"
    ),

  tabella_probabilita |>
    transmute(
      posizione = 1 - soglia_x,
      probabilita = probabilita_upper,
      coda = "Coda superiore"
    )

) |>
  mutate(
    coda = factor(
      coda,
      levels = livelli_coda
    )
  )


dati_finestre_curva_intera <- tibble::tibble(
  coda = factor(
    livelli_coda,
    levels = livelli_coda
  ),

  xmin = c(
    x_min_stima,
    1 - x_max_stima
  ),

  xmax = c(
    x_max_stima,
    1 - x_min_stima
  )
)


# 7.6 Confronto tra probabilità grezza e GTCF normalizzata
# Entrambe le curve vengono indicizzate rispetto al valore osservato
# nella soglia massima della finestra, così da confrontarne soltanto la forma.

gtcf_overlay_lower <- tabella_gtcf |>
  filter(
    soglia_x >= x_min_stima,
    soglia_x <= x_max_stima
  ) |>
  pull(gtcf_empirica_lower)


gtcf_overlay_upper <- tabella_gtcf |>
  filter(
    soglia_x >= x_min_stima,
    soglia_x <= x_max_stima
  ) |>
  pull(gtcf_empirica_upper)


dati_overlay_lower <- tibble::tibble(
  soglia_x = tabella_probabilita_stima$soglia_x,

  probabilita_grezza =
    tabella_probabilita_stima$probabilita_lower /
    tail(
      tabella_probabilita_stima$probabilita_lower,
      1
    ),

  gtcf_normalizzata =
    gtcf_overlay_lower /
    tail(
      gtcf_overlay_lower,
      1
    ),

  coda = "Coda inferiore"
)


dati_overlay_upper <- tibble::tibble(
  soglia_x = tabella_probabilita_stima$soglia_x,

  probabilita_grezza =
    tabella_probabilita_stima$probabilita_upper /
    tail(
      tabella_probabilita_stima$probabilita_upper,
      1
    ),

  gtcf_normalizzata =
    gtcf_overlay_upper /
    tail(
      gtcf_overlay_upper,
      1
    ),

  coda = "Coda superiore"
)


dati_overlay <- bind_rows(
  dati_overlay_lower,
  dati_overlay_upper
) |>
  pivot_longer(
    cols = c(
      probabilita_grezza,
      gtcf_normalizzata
    ),
    names_to = "serie",
    values_to = "valore"
  ) |>
  mutate(
    coda = factor(
      coda,
      levels = livelli_coda
    ),

    serie = recode(
      serie,
      probabilita_grezza =
        "Probabilità grezza indicizzata",
      gtcf_normalizzata =
        "GTCF normalizzata indicizzata"
    )
  )



# 8. Sistema grafico comune -----------------------------------------------------

# 8.1 Palette
# Colori comuni utilizzati nei grafici dell'analisi.

col_lower <- "#67A9CF"
col_upper <- "#EF8A62"
col_teorico <- "#333333"
col_asintotico <- "#999999"
col_riferimento <- "#9AA1A8"
col_finestra <- "#F7F7F7"
col_griglia <- "#E6E6E6"
col_testo <- "#222222"

palette_code <- c(
  "Coda inferiore" = col_lower,
  "Coda superiore" = col_upper
)


# 8.2 Tema grafico comune
# Definisce lo stile uniforme utilizzato in tutti i grafici.

tema_tesi <- function(base_size = 12) {

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
      panel.grid.major.x = element_blank(),

      panel.grid.major.y = element_line(
        color = col_griglia,
        linewidth = 0.35
      ),

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


# 8.3 Finestre di stima
# La stessa finestra finale viene utilizzata nelle due code.

dati_finestre_stima <- tibble::tibble(
  coda = factor(
    c(
      "Coda inferiore",
      "Coda superiore"
    ),
    levels = livelli_coda
  ),

  soglia_minima = c(
    x_min_stima,
    x_min_stima
  ),

  soglia_massima = x_max_stima
)


# 8.4 Evidenziazione delle finestre di stima
# Aggiunge area ombreggiata e limite superiore della finestra ai grafici.

evidenzia_finestre_stima <- function() {

  list(
    geom_rect(
      data = dati_finestre_stima,
      aes(
        xmin = soglia_minima,
        xmax = soglia_massima,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = col_finestra,
      alpha = 0.72
    ),

    geom_vline(
      data = dati_finestre_stima,
      aes(
        xintercept = soglia_massima
      ),
      inherit.aes = FALSE,
      color = col_riferimento,
      linetype = "22",
      linewidth = 0.55
    )
  )
}



# 9. Grafici -------------------------------------------------------------------

# 9.1 Probabilità congiunta lungo entrambe le code

grafico_curva_intera <- function() {

  minimo_positivo <- min(
    dati_curva_intera$probabilita[
      dati_curva_intera$probabilita > 0
    ]
  ) * 0.80

  ggplot() +
    geom_rect(
      data = dati_finestre_curva_intera,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = minimo_positivo,
        ymax = Inf,
        fill = coda
      ),
      inherit.aes = FALSE,
      alpha = 0.08,
      show.legend = FALSE
    ) +
    geom_vline(
      xintercept = 0.50,
      color = col_riferimento,
      linetype = "22",
      linewidth = 0.70
    ) +
    geom_line(
      data = dati_curva_intera,
      aes(
        x = posizione,
        y = probabilita,
        color = coda,
        group = coda
      ),
      linewidth = 1.00
    ) +
    geom_point(
      data = dati_curva_intera,
      aes(
        x = posizione,
        y = probabilita,
        color = coda
      ),
      size = 1.75,
      alpha = 0.85
    ) +
    scale_color_manual(
      values = palette_code
    ) +
    scale_fill_manual(
      values = palette_code
    ) +
    scale_x_continuous(
      limits = c(
        x_min_stima,
        1 - x_min_stima
      ),
      breaks = c(
        0.01,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
        0.99
      ),
      labels = scales::label_number(
        accuracy = 0.01,
        trim = TRUE
      ),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_log10(
      breaks = scales::breaks_log(n = 6),
      labels = scales::label_number(
        accuracy = 0.001,
        trim = TRUE
      ),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    labs(
      title = "Probabilità congiunta lungo entrambe le code",
      subtitle = paste0(
        "Dal centro x = 0.50, la lower tail si sviluppa verso 0 ",
        "e la upper tail verso 1."
      ),
      x = "Soglia sul quadrato unitario",
      y = "Probabilità congiunta empirica",
      caption = paste0(
        "Le aree evidenziate identificano la finestra [0.01; 0.10] ",
        "utilizzata per stimare α. ",
        "La linea verticale indica il centro x = 0.50."
      )
    ) +
    tema_tesi()
}


# 9.2 Probabilità congiunte nella finestra di stima

grafico_probabilita_loglog <- function() {

  etichette_facet <- c(
    "Coda inferiore" = sprintf(
      "Coda inferiore  |  α̂ = %.4f",
      alpha_empirico_lower
    ),
    "Coda superiore" = sprintf(
      "Coda superiore  |  α̂ = %.4f",
      alpha_empirico_upper
    )
  )

  ggplot() +
    geom_line(
      data = dati_probabilita_grafico |>
        filter(serie == "Teorica"),
      aes(
        x = soglia_x,
        y = probabilita,
        group = coda
      ),
      color = col_teorico,
      linewidth = 0.90
    ) +
    geom_line(
      data = dati_probabilita_grafico |>
        filter(serie == "Potenza stimata"),
      aes(
        x = soglia_x,
        y = probabilita,
        color = coda,
        group = coda
      ),
      linewidth = 0.85,
      linetype = "22"
    ) +
    geom_point(
      data = dati_probabilita_grafico |>
        filter(serie == "Empirica"),
      aes(
        x = soglia_x,
        y = probabilita,
        color = coda
      ),
      size = 2.05,
      alpha = 0.88
    ) +
    facet_wrap(
      vars(coda),
      nrow = 1,
      labeller = as_labeller(etichette_facet)
    ) +
    scale_color_manual(
      values = palette_code,
      guide = "none"
    ) +
    scale_x_log10(
      limits = c(x_min_stima, x_max_stima),
      breaks = c(
        0.01,
        0.02,
        0.03,
        0.05,
        0.10
      ),
      labels = scales::label_number(
        accuracy = 0.001,
        trim = TRUE
      ),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    scale_y_log10(
      breaks = scales::breaks_log(n = 6),
      labels = scales::label_number(
        accuracy = 0.001,
        trim = TRUE
      ),
      expand = expansion(mult = c(0.05, 0.09))
    ) +
    labs(
      title = "Probabilità congiunte nella finestra di stima",
      subtitle = sprintf(
        paste0(
          "Punti = probabilità empiriche; linea scura = teoriche; ",
          "linea tratteggiata = potenza stimata. ",
          "n = %s, τ = %.3f, ρ = %.2f."
        ),
        format(
          n_osservazioni,
          big.mark = ".",
          decimal.mark = ",",
          scientific = FALSE
        ),
        tau_kendall,
        rho
      ),
      x = "Soglia x",
      y = "Probabilità congiunta",
      caption = sprintf(
        "Finestra comune di regressione lower e upper: [%.2f; %.2f].",
        x_min_stima,
        x_max_stima
      )
    ) +
    tema_tesi()
}


# 9.3 GTCF empirica dopo la normalizzazione

grafico_gtcf_normalizzazione <- function() {

  etichette_facet <- c(
    "Coda inferiore" = sprintf(
      "Coda inferiore  |  α̂ = %.4f  |  ĉ = %.4f",
      alpha_empirico_lower,
      livello_gtcf_empirico_lower
    ),

    "Coda superiore" = sprintf(
      "Coda superiore  |  α̂ = %.4f  |  ĉ = %.4f",
      alpha_empirico_upper,
      livello_gtcf_empirico_upper
    )
  )

  ggplot(
    data = dati_gtcf_finale_empirica,
    aes(
      x = soglia_x,
      y = gtcf,
      color = coda,
      group = coda
    )
  ) +
    geom_hline(
      data = dati_riferimenti_gtcf |>
        filter(tipo == "Livello stimato"),
      aes(yintercept = livello),
      inherit.aes = FALSE,
      color = col_teorico,
      linewidth = 0.95
    ) +
    geom_hline(
      data = dati_riferimenti_gtcf |>
        filter(tipo == "Benchmark di finestra"),
      aes(yintercept = livello),
      inherit.aes = FALSE,
      color = col_riferimento,
      linetype = "22",
      linewidth = 0.85
    ) +
    geom_line(linewidth = 0.95) +
    geom_point(
      size = 2.0,
      alpha = 0.88
    ) +
    facet_wrap(
      vars(coda),
      nrow = 1,
      labeller = as_labeller(etichette_facet)
    ) +
    scale_color_manual(
      values = palette_code,
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(
        x_min_stima,
        x_max_stima
      ),
      breaks = c(
        0.01,
        0.03,
        0.05,
        0.075,
        0.10
      ),
      labels = scales::label_number(
        accuracy = 0.01,
        trim = TRUE
      ),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.01,
        trim = TRUE
      ),
      expand = expansion(mult = c(0.06, 0.08))
    ) +
    labs(
      title = "GTCF empirica dopo la normalizzazione",
      subtitle = paste0(
        "La normalizzazione con α̂ rimuove la componente di potenza ",
        "e isola il livello locale stimato nella finestra."
      ),
      x = "Soglia x",
      y = "GTCF empirica",
      caption = paste0(
        "Linea continua scura = livello stimato; ",
        "linea tratteggiata = benchmark teorico di finestra. ",
        "Con la sola normalizzazione di potenza non esiste ",
        "un livello asintotico positivo."
      )
    ) +
    tema_tesi()
}


# 9.4 Confronto tra curva grezza e GTCF normalizzata

grafico_overlay_normalizzazione <- function() {

  colori_overlay <- c(
    "Probabilità grezza indicizzata" = col_riferimento,
    "GTCF normalizzata indicizzata" = col_teorico
  )

  linee_overlay <- c(
    "Probabilità grezza indicizzata" = "solid",
    "GTCF normalizzata indicizzata" = "22"
  )

  ggplot(
    data = dati_overlay,
    aes(
      x = soglia_x,
      y = valore,
      color = serie,
      linetype = serie,
      group = serie
    )
  ) +
    evidenzia_finestre_stima() +
    geom_hline(
      yintercept = 1,
      color = col_asintotico,
      linetype = "13",
      linewidth = 0.80
    ) +
    geom_line(
      linewidth = 1.00
    ) +
    geom_point(
      size = 1.9,
      alpha = 0.88
    ) +
    facet_wrap(
      vars(coda),
      nrow = 1
    ) +
    scale_color_manual(
      values = colori_overlay,
      breaks = names(colori_overlay)
    ) +
    scale_linetype_manual(
      values = linee_overlay,
      breaks = names(linee_overlay)
    ) +
    scale_x_continuous(
      limits = c(
        x_min_stima,
        x_max_stima
      ),
      breaks = c(
        0.01,
        0.02,
        0.03,
        0.05,
        0.08,
        0.10
      ),
      labels = scales::label_number(
        accuracy = 0.01,
        trim = TRUE
      ),
      expand = expansion(
        mult = c(
          0.01,
          0.01
        )
      )
    ) +
    scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.1,
        trim = TRUE
      ),
      expand = expansion(
        mult = c(
          0.06,
          0.08
        )
      )
    ) +
    labs(
      title =
        "Confronto tra probabilità grezza e GTCF normalizzata",

      subtitle = paste0(
        "Entrambe le serie sono indicizzate a 1 in x = 0.10 ",
        "per confrontarne la forma nella finestra di stima."
      ),

      x = "Soglia x",
      y = "Livello relativo",

      caption = paste0(
        "Spostandosi verso sinistra, x diminuisce e la soglia diventa ",
        "più estrema. Il confronto riguarda la forma, non i livelli assoluti."
      )
    ) +
    tema_tesi()
}


# 9.5 GTCF teorica: alpha di finestra e alpha asintotico

grafico_gtcf_teorica <- function() {

  colori_normalizzazione <- c(
    "Alpha di finestra" = col_teorico,
    "Alpha asintotico" = col_asintotico
  )

  linee_normalizzazione <- c(
    "Alpha di finestra" = "solid",
    "Alpha asintotico" = "22"
  )

  dati_finestra_teorica <- dati_finestre_stima |>
    distinct(
      soglia_minima,
      soglia_massima
    )

  ggplot(
    data = dati_gtcf_alpha_grafico,
    aes(
      x = soglia_x,
      y = gtcf,
      color = normalizzazione,
      linetype = normalizzazione
    )
  ) +
    geom_rect(
      data = dati_finestra_teorica,
      aes(
        xmin = soglia_minima,
        xmax = soglia_massima,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = col_finestra,
      alpha = 0.72
    ) +
    geom_vline(
      data = dati_finestra_teorica,
      aes(
        xintercept = soglia_massima
      ),
      inherit.aes = FALSE,
      color = col_riferimento,
      linetype = "22",
      linewidth = 0.55
    ) +
    geom_line(linewidth = 0.95) +
    scale_color_manual(
      values = colori_normalizzazione
    ) +
    scale_linetype_manual(
      values = linee_normalizzazione
    ) +
    scale_x_continuous(
      limits = c(x_min_stima, x_max_gtcf),
      breaks = c(
        0.01,
        0.10,
        0.20,
        0.30,
        0.40,
        0.50
      ),
      labels = scales::label_number(
        accuracy = 0.01,
        trim = TRUE
      ),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.01,
        trim = TRUE
      ),
      expand = expansion(mult = c(0.06, 0.08))
    ) +
    labs(
      title = "Effetto della normalizzazione sulla GTCF teorica",
      subtitle = paste0(
        "Confronto tra il benchmark efficace nella finestra osservabile ",
        "e il benchmark asintotico comune alle due code."
      ),
      x = "Soglia x",
      y = "GTCF teorica",
      caption = sprintf(
        paste0(
          "Gaussiana con τ = %.3f e ρ = %.2f: ",
          "α di finestra = %.4f; α asintotico = %.2f."
        ),
        tau_kendall,
        rho,
        alpha_teorico_finestra_lower,
        alpha_asintotico_lower
      )
    ) +
    tema_tesi()
}


# 10. Creazione degli oggetti grafici ------------------------------------------

# Crea i cinque grafici finali e li raccoglie nella lista
# utilizzata per il salvataggio in PNG e PDF.
lista_grafici <- list(
  curva_intera = grafico_curva_intera(),
  probabilita_loglog = grafico_probabilita_loglog(),
  gtcf_normalizzazione = grafico_gtcf_normalizzazione(),
  overlay_normalizzazione = grafico_overlay_normalizzazione(),
  gtcf_teorica = grafico_gtcf_teorica()
)



# 11. Salvataggio degli output -------------------------------------------------

# 11.1 Tabelle CSV

# Unisce probabilità di coda e GTCF calcolate sulle stesse soglie.
tabella_output_soglie <- tabella_probabilita |>
  left_join(
    tabella_gtcf,
    by = "soglia_x"
  )


write.csv(
  tabella_output_soglie,
  file = "output/Gaussiana/Analisi Gaussiana/Probabilita_GTCF_Gaussiana.csv",
  row.names = FALSE
)

write.csv(
  tabella_risultati_finali,
  file = "output/Gaussiana/Analisi Gaussiana/Risultati_Finali_Gaussiana.csv",
  row.names = FALSE
)


# 11.2 Funzione di salvataggio dei grafici

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


# 11.3 Specifiche dei file PNG

specifiche_grafici <- tibble::tribble(
  ~nome_grafico, ~nome_file, ~larghezza, ~altezza,

  "curva_intera",
  "output/Gaussiana/Analisi Gaussiana/Curva_Intera_Gaussiana.png",
  12, 6.5,

  "probabilita_loglog",
  "output/Gaussiana/Analisi Gaussiana/Probabilita_LogLog_Gaussiana.png",
  12, 6.5,

  "gtcf_normalizzazione",
  "output/Gaussiana/Analisi Gaussiana/GTCF_Normalizzazione_Gaussiana.png",
  12, 6.5,

  "overlay_normalizzazione",
  "output/Gaussiana/Analisi Gaussiana/Overlay_Normalizzazione_Gaussiana.png",
  12, 6.5,

  "gtcf_teorica",
  "output/Gaussiana/Analisi Gaussiana/GTCF_Teorica_Gaussiana.png",
  12, 6.5
)


for (i in seq_len(nrow(specifiche_grafici))) {
  salva_grafico(
    grafico = lista_grafici[[specifiche_grafici$nome_grafico[[i]]]],
    nome_file = specifiche_grafici$nome_file[[i]],
    larghezza = specifiche_grafici$larghezza[[i]],
    altezza = specifiche_grafici$altezza[[i]]
  )
}


# 11.4 PDF multipagina

grDevices::cairo_pdf(
  filename = "output/Gaussiana/Analisi Gaussiana/Grafici_Analisi_Gaussiana.pdf",
  width = 12,
  height = 7,
  onefile = TRUE
)

for (grafico in lista_grafici) {
  print(grafico)
}

dev.off()



# 12. Riepilogo finale ----------------------------------------------------------

tabella_riepilogo_finale <- tabella_risultati_finali |>
  select(
    coda,
    numerosita,
    soglia_minima,
    soglia_massima,
    beta_empirico,
    alpha_empirico,
    alpha_teorico_finestra,
    alpha_asintotico,
    livello_gtcf_empirico,
    livello_gtcf_teorico_finestra,
    livello_gtcf_asintotico
  )

cat("\nRisultati finali dell'analisi Gaussiana\n")

print(
  tabella_riepilogo_finale,
  n = Inf,
  width = Inf
)

cat(
  "\nRisultati empirici principali:",
  "\nLower: alpha =", round(alpha_empirico_lower, 6),
  "| livello GTCF stimato nella finestra =", round(livello_gtcf_empirico_lower, 6),
  "\nUpper: alpha =", round(alpha_empirico_upper, 6),
  "| livello GTCF stimato nella finestra =", round(livello_gtcf_empirico_upper, 6),
  "\n"
)

cat(
  "\nFinestra finale:",
  "\nLower = [",
  x_min_stima,
  ", ",
  x_max_stima,
  "]",
  "\nUpper = [",
  x_min_stima,
  ", ",
  x_max_stima,
  "]\n",
  sep = ""
)

cat(
  "\nParametri della copula:",
  "\nKendall tau =", round(tau_kendall, 6),
  "\nRho gaussiano =", round(rho, 6),
  "\n"
)

cat(
  "\nNota metodologica:",
  "\nIl livello GTCF stimato è un livello locale della finestra.",
  "\nPer la copula gaussiana non è definito un limite asintotico positivo.",
  "\n"
)

cat("\nEsecuzione completata: tabelle e grafici salvati.\n")
