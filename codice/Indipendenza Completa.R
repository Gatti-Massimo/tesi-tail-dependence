# ============================================================
# GTCF - INDIPENDENZA, LOWER E UPPER TAIL
# Stima della TCF mediante curve di potenza e calcolo di alpha*
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
# 3. COPULA EMPIRICA E TCF - LOWER E UPPER TAIL
# ------------------------------------------------------------

# LOWER TAIL
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


# UPPER TAIL
# P_U(x) = P(U > 1-x, V > 1-x)
#        = 2x - 1 + C(1-x, 1-x)

copula_upper_empirica <- EMPIRcop(
  u = 1 - soglie,
  v = 1 - soglie,
  para = osservazioni,
  ctype = tipo_copula_empirica
)

probabilita_upper_empirica <-
  2 * soglie - 1 + copula_upper_empirica

conteggi_upper <- round(
  probabilita_upper_empirica * n_osservazioni
)

tcf_upper_empirica <-
  probabilita_upper_empirica / soglie


# Valori teorici sotto indipendenza

probabilita_lower_teorica <- P(
  u = soglie,
  v = soglie
)

copula_upper_teorica <- P(
  u = 1 - soglie,
  v = 1 - soglie
)

probabilita_upper_teorica <-
  2 * soglie - 1 + copula_upper_teorica

tcf_lower_teorica <-
  probabilita_lower_teorica / soglie

tcf_upper_teorica <-
  probabilita_upper_teorica / soglie


# Tabella delle due code

tabella_code <- data.frame(
  x = soglie,
  soglia_lower = soglie,
  soglia_upper = 1 - soglie,
  conteggi_lower = conteggi_lower,
  conteggi_upper = conteggi_upper,
  probabilita_lower_empirica = probabilita_lower_empirica,
  probabilita_upper_empirica = probabilita_upper_empirica,
  tcf_lower_empirica = tcf_lower_empirica,
  tcf_upper_empirica = tcf_upper_empirica,
  tcf_lower_teorica = tcf_lower_teorica,
  tcf_upper_teorica = tcf_upper_teorica
)


# ------------------------------------------------------------
# 4. REGIONI INFORMATIVE E NORMALIZZAZIONE
# ------------------------------------------------------------

dati_stima_lower <- subset(
  tabella_code,
  conteggi_lower >= min_eventi_congiunti
)

dati_stima_upper <- subset(
  tabella_code,
  conteggi_upper >= min_eventi_congiunti
)

# In x = 0.5 lower e upper coincidono.
indice_riferimento <- which.max(tabella_code$x)

tcf_riferimento <-
  tabella_code$tcf_lower_empirica[indice_riferimento]

dati_stima_lower$tcf_normalizzata <-
  dati_stima_lower$tcf_lower_empirica / tcf_riferimento

dati_stima_upper$tcf_normalizzata <-
  dati_stima_upper$tcf_upper_empirica / tcf_riferimento


# ------------------------------------------------------------
# 5. STIMA DELLE CURVE LOWER E UPPER
# ------------------------------------------------------------

stima_lower <- stima_curva_potenza(
  x = dati_stima_lower$x,
  tcf_normalizzata = dati_stima_lower$tcf_normalizzata
)

stima_upper <- stima_curva_potenza(
  x = dati_stima_upper$x,
  tcf_normalizzata = dati_stima_upper$tcf_normalizzata
)

stopifnot(
  stima_lower$convergence == 0,
  stima_upper$convergence == 0
)

a_lower <- unname(stima_lower$par["a"])
gamma_lower <- unname(stima_lower$par["gamma"])

a_upper <- unname(stima_upper$par["a"])
gamma_upper <- unname(stima_upper$par["gamma"])

# Curve riportate sulla scala originale della TCF.

tabella_code$tcf_lower_stimata <-
  tcf_riferimento *
  curva_potenza(
    x = tabella_code$x,
    a = a_lower,
    gamma = gamma_lower
  )

tabella_code$tcf_upper_stimata <-
  tcf_riferimento *
  curva_potenza(
    x = tabella_code$x,
    a = a_upper,
    gamma = gamma_upper
  )


# ------------------------------------------------------------
# 6. CALCOLO DI ALPHA* PER LOWER E UPPER TAIL
# ------------------------------------------------------------

# LOWER TAIL

x_lower <- dati_stima_lower$x

tcf_lower_regolarizzata <-
  tabella_code$tcf_lower_stimata[
    match(x_lower, tabella_code$x)
  ]

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


# UPPER TAIL

x_upper <- dati_stima_upper$x

tcf_upper_regolarizzata <-
  tabella_code$tcf_upper_stimata[
    match(x_upper, tabella_code$x)
  ]

probabilita_upper_regolarizzata <-
  x_upper * tcf_upper_regolarizzata

risultato_alpha_upper <- trova_alpha_star(
  probabilita_congiunta = probabilita_upper_regolarizzata,
  soglia_x = x_upper,
  griglia_alpha = griglia_alpha,
  tolleranza = tolleranza_monotonia
)

alpha_star_upper <-
  risultato_alpha_upper$alpha_star

gtcf_alpha_star_upper <-
  risultato_alpha_upper$gtcf_alpha_star


# Tabelle finali delle due code

tabella_alpha_lower <- data.frame(
  x = x_lower,
  probabilita_regolarizzata = probabilita_lower_regolarizzata,
  gtcf_alpha_star = gtcf_alpha_star_lower
)

tabella_alpha_upper <- data.frame(
  x = x_upper,
  probabilita_regolarizzata = probabilita_upper_regolarizzata,
  gtcf_alpha_star = gtcf_alpha_star_upper
)


# ------------------------------------------------------------
# 7. GRAFICO DELLA TCF COMPLETA
# ------------------------------------------------------------

# La soglia upper viene invertita per rappresentarla
# da t = 0.5 a t = 0.99.
# Il valore t = 0.5 non viene duplicato.

soglia_completa <- c(
  tabella_code$soglia_lower,
  rev(tabella_code$soglia_upper)[-1]
)

tcf_empirica_completa <- c(
  tabella_code$tcf_lower_empirica,
  rev(tabella_code$tcf_upper_empirica)[-1]
)

tcf_stimata_completa <- c(
  tabella_code$tcf_lower_stimata,
  rev(tabella_code$tcf_upper_stimata)[-1]
)

tcf_teorica_completa <- c(
  tabella_code$tcf_lower_teorica,
  rev(tabella_code$tcf_upper_teorica)[-1]
)


# Cartella e percorso di salvataggio

cartella_grafici <- "grafici"

dir.create(
  path = cartella_grafici,
  recursive = TRUE,
  showWarnings = FALSE
)

percorso_grafico <- file.path(
  cartella_grafici,
  "TCF_indipendenza_completa.png"
)


# Funzione grafica

disegna_tcf_completa <- function() {

  par(
    pty = "s",
    mar = c(5, 5, 4, 2)
  )

  plot(
    x = soglia_completa,
    y = tcf_empirica_completa,
    type = "n",
    xlim = c(0, 1),
    ylim = c(0, 0.5),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "Soglia t",
    ylab = "Tail Concentration Function",
    main = "TCF completa nel caso di indipendenza"
  )

  axis(
    side = 1,
    at = seq(0, 1, by = 0.1)
  )

  axis(
    side = 2,
    at = seq(0, 0.5, by = 0.1),
    las = 1
  )

  box()

  abline(
    v = 0.5,
    col = "grey75",
    lty = 3
  )

  lines(
    x = soglia_completa,
    y = tcf_teorica_completa,
    col = "dodgerblue3",
    lwd = 2.5,
    lty = 2
  )

  lines(
    x = soglia_completa,
    y = tcf_stimata_completa,
    col = "red3",
    lwd = 3
  )

  points(
    x = soglia_completa,
    y = tcf_empirica_completa,
    col = "black",
    pch = 16,
    cex = 0.8
  )

  legend(
    "topleft",
    legend = c(
      "TCF empirica",
      "Curve di potenza stimate",
      "TCF teorica"
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

disegna_tcf_completa()


# Salvataggio del grafico

png(
  filename = percorso_grafico,
  width = 1800,
  height = 1800,
  res = 250
)

disegna_tcf_completa()

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

gtcf_upper_non_decrescente <- if (!is.na(alpha_star_upper)) {

  controlla_monotonia(
    curva = gtcf_alpha_star_upper,
    tolleranza = tolleranza_monotonia
  )

} else {

  NA
}


cat(
  "\nRISULTATI DELLA STIMA COMPLETA\n",
  "==============================\n",

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

  "\nUPPER TAIL\n",
  "----------\n",
  "Parametro a stimato: ",
  round(a_upper, 6),
  "\n",
  "Parametro gamma stimato: ",
  round(gamma_upper, 6),
  "\n",
  "Convergenza di optim: ",
  stima_upper$convergence,
  "\n",
  "Numero di soglie informative: ",
  nrow(dati_stima_upper),
  "\n",
  "Alpha* stimato: ",
  if (is.na(alpha_star_upper)) {
    NA
  } else {
    round(alpha_star_upper, 6)
  },
  "\n",
  "Alpha* teorico sotto indipendenza: 0.5\n",
  "GTCF non decrescente: ",
  gtcf_upper_non_decrescente,
  "\n",

  "\nVALORE COMUNE IN x = 0.5\n",
  "-------------------------\n",
  "TCF empirica di riferimento: ",
  round(tcf_riferimento, 6),
  "\n",

  "\nGRAFICO SALVATO IN\n",
  "------------------\n",
  percorso_grafico,
  "\n",
  sep = ""
)
