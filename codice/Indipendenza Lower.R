# ============================================================
# GTCF - INDIPENDENZA, LOWER TAIL
# Stima della TCF mediante curva di potenza e calcolo di alpha*
# ============================================================


# ------------------------------------------------------------
# 1. IMPOSTAZIONI INIZIALI
# ------------------------------------------------------------

rm(list = ls())

library(copBasic)

source("codice/funzioni/curva_potenza.R")
source("codice/funzioni/stima_curva_potenza.R")
source("codice/funzioni/calcola_gtcf.R")
source("codice/funzioni/controlla_monotonia.R")
source("codice/funzioni/trova_alpha_star.R")

set.seed(44)

n_osservazioni <- 10000

soglie <- seq(
  from = 0.01,
  to = 0.50,
  by = 0.01
)

min_eventi_congiunti <- 20

griglia_alpha <- seq(
  from = 0.5,
  to = 1,
  by = 0.0001
)

tolleranza_monotonia <- 1e-10

tipo_copula_empirica <- "weibull"


# ------------------------------------------------------------
# 2. SIMULAZIONE DALLA COPULA DI INDIPENDENZA
# ------------------------------------------------------------

osservazioni <- simCOP(
  n = n_osservazioni,
  cop = P,
  graphics = FALSE
)

correlazione_campione <- cor(
  osservazioni$U,
  osservazioni$V
)


# ------------------------------------------------------------
# 3. COPULA EMPIRICA E TCF - LOWER TAIL
# ------------------------------------------------------------

# P_L(x) = P(U <= x, V <= x) = C(x, x)

probabilita_lower_empirica <- EMPIRcop(
  u = soglie,
  v = soglie,
  para = osservazioni,
  ctype = tipo_copula_empirica
)

conteggi_lower <- round(
  probabilita_lower_empirica * n_osservazioni
)

tcf_lower_empirica <-
  probabilita_lower_empirica / soglie


# Valori teorici sotto indipendenza

probabilita_lower_teorica <- P(
  u = soglie,
  v = soglie
)

tcf_lower_teorica <-
  probabilita_lower_teorica / soglie


# Tabella della lower tail

tabella_lower <- data.frame(
  x = soglie,
  conteggi_lower = conteggi_lower,
  probabilita_lower_empirica = probabilita_lower_empirica,
  probabilita_lower_teorica = probabilita_lower_teorica,
  tcf_lower_empirica = tcf_lower_empirica,
  tcf_lower_teorica = tcf_lower_teorica
)


# ------------------------------------------------------------
# 4. REGIONE INFORMATIVA E NORMALIZZAZIONE
# ------------------------------------------------------------

dati_stima_lower <- subset(
  tabella_lower,
  conteggi_lower >= min_eventi_congiunti
)

indice_riferimento <- which.max(
  tabella_lower$x
)

tcf_riferimento <-
  tabella_lower$tcf_lower_empirica[indice_riferimento]

dati_stima_lower$tcf_normalizzata <-
  dati_stima_lower$tcf_lower_empirica / tcf_riferimento


# ------------------------------------------------------------
# 5. STIMA DELLA CURVA LOWER
# ------------------------------------------------------------

stima_lower <- stima_curva_potenza(
  x = dati_stima_lower$x,
  tcf_normalizzata = dati_stima_lower$tcf_normalizzata
)

stopifnot(
  stima_lower$convergence == 0
)

a_lower <- unname(
  stima_lower$par["a"]
)

gamma_lower <- unname(
  stima_lower$par["gamma"]
)

# Curva riportata sulla scala originale della TCF.

tabella_lower$tcf_lower_stimata <-
  tcf_riferimento *
  curva_potenza(
    x = tabella_lower$x,
    a = a_lower,
    gamma = gamma_lower
  )


# ------------------------------------------------------------
# 6. CALCOLO DI ALPHA* PER LA LOWER TAIL
# ------------------------------------------------------------

x_lower <- dati_stima_lower$x

tcf_lower_regolarizzata <-
  tabella_lower$tcf_lower_stimata[
    match(x_lower, tabella_lower$x)
  ]

# q_L(x) = P_L(x) / x,
# quindi P_L(x) = x * q_L(x).

probabilita_lower_regolarizzata <-
  x_lower * tcf_lower_regolarizzata

risultato_alpha_lower <- trova_alpha_star(
  probabilita_congiunta = probabilita_lower_regolarizzata,
  soglia_x = x_lower,
  griglia_alpha = griglia_alpha,
  tolleranza = tolleranza_monotonia
)

alpha_star_lower <-
  risultato_alpha_lower$alpha_star

gtcf_alpha_star_lower <-
  risultato_alpha_lower$gtcf_alpha_star


# Tabella finale della lower tail

tabella_alpha_lower <- data.frame(
  x = x_lower,
  probabilita_regolarizzata = probabilita_lower_regolarizzata,
  gtcf_alpha_star = gtcf_alpha_star_lower
)


# ------------------------------------------------------------
# 7. GRAFICO DELLA TCF LOWER
# ------------------------------------------------------------

# Cartella e percorso di salvataggio

cartella_grafici <- "grafici"

dir.create(
  path = cartella_grafici,
  recursive = TRUE,
  showWarnings = FALSE
)

percorso_grafico <- file.path(
  cartella_grafici,
  "TCF_indipendenza_lower.png"
)


# Funzione grafica

disegna_tcf_lower <- function() {

  par(
    pty = "s",
    mar = c(5, 5, 4, 2)
  )

  plot(
    x = tabella_lower$x,
    y = tabella_lower$tcf_lower_empirica,
    type = "n",
    xlim = c(0, 0.5),
    ylim = c(0, 0.5),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "Soglia x",
    ylab = "Tail Concentration Function",
    main = "TCF lower nel caso di indipendenza"
  )

  axis(
    side = 1,
    at = seq(0, 0.5, by = 0.1)
  )

  axis(
    side = 2,
    at = seq(0, 0.5, by = 0.1),
    las = 1
  )

  box()

  lines(
    x = tabella_lower$x,
    y = tabella_lower$tcf_lower_teorica,
    col = "dodgerblue3",
    lwd = 2.5,
    lty = 2
  )

  lines(
    x = tabella_lower$x,
    y = tabella_lower$tcf_lower_stimata,
    col = "red3",
    lwd = 3
  )

  points(
    x = tabella_lower$x,
    y = tabella_lower$tcf_lower_empirica,
    col = "black",
    pch = 16,
    cex = 0.8
  )

  legend(
    "topleft",
    legend = c(
      "TCF empirica",
      "Curva di potenza stimata",
      "TCF teorica: q(x) = x"
    ),
    col = c(
      "black",
      "red3",
      "dodgerblue3"
    ),
    pch = c(16, NA, NA),
    lty = c(NA, 1, 2),
    lwd = c(NA, 3, 2.5),
    bty = "n",
    cex = 0.85
  )
}


# Apertura in una finestra grafica separata

if (.Platform$OS.type == "windows") {

  windows(width = 7, height = 7)

} else if (Sys.info()[["sysname"]] == "Darwin") {

  quartz(width = 7, height = 7)

} else {

  x11(width = 7, height = 7)
}

disegna_tcf_lower()


# Salvataggio del grafico

png(
  filename = percorso_grafico,
  width = 1800,
  height = 1800,
  res = 250
)

disegna_tcf_lower()

dev.off()


# ------------------------------------------------------------
# 8. RISULTATI FINALI
# ------------------------------------------------------------

gtcf_lower_non_decrescente <- if (!is.na(alpha_star_lower)) {

  controlla_monotonia(
    curva = gtcf_alpha_star_lower,
    tolleranza = tolleranza_monotonia
  )

} else {

  NA
}


cat(
  "\nRISULTATI DELLA STIMA LOWER\n",
  "===========================\n",

  "\nSIMULAZIONE\n",
  "-----------\n",
  "Correlazione campionaria tra U e V: ",
  round(correlazione_campione, 6),
  "\n",

  "\nLOWER TAIL\n",
  "----------\n",
  "Parametro a stimato: ",
  round(a_lower, 6),
  "\n",
  "Parametro gamma stimato: ",
  round(gamma_lower, 6),
  "\n",
  "Convergenza di optim: ",
  stima_lower$convergence,
  "\n",
  "Numero di soglie informative: ",
  nrow(dati_stima_lower),
  "\n",
  "Alpha* stimato: ",
  if (is.na(alpha_star_lower)) {
    NA
  } else {
    round(alpha_star_lower, 6)
  },
  "\n",
  "Alpha* teorico sotto indipendenza: 0.5\n",
  "GTCF non decrescente: ",
  gtcf_lower_non_decrescente,
  "\n",

  "\nVALORE IN x = 0.5\n",
  "------------------\n",
  "TCF empirica di riferimento: ",
  round(tcf_riferimento, 6),
  "\n",

  "\nGRAFICO SALVATO IN\n",
  "------------------\n",
  percorso_grafico,
  "\n",
  sep = ""
)
