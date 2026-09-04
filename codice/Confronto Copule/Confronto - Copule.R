# Confronto - Copule Gaussiana, Clayton e Gumbel


# 1. Impostazioni --------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

ordine_copule <- c("Gaussiana", "Clayton", "Gumbel")
ordine_code <- c("Lower", "Upper")
livelli_coda <- c("Coda inferiore", "Coda superiore")
tau_attesi <- round(c(0.20, 1 / 3, 0.50), 12)
cartella_output <- file.path("output", "Confronto Copule")

parametri_copule <- c(
  "Gaussiana" = "rho",
  "Clayton" = "theta",
  "Gumbel" = "theta"
)



# 2. Importazione automatica ---------------------------------------------------

# Per ogni copula importa solo i tre CSV necessari al confronto.
# I nomi dei file e delle cartelle sono costruiti automaticamente.

importa_copula <- function(copula) {

  cartella_analisi <- file.path(
    "output", copula, paste("Analisi", copula)
  )

  cartella_montecarlo <- file.path(
    "output", copula, paste("Monte Carlo", copula)
  )

  list(
    probabilita = read.csv(
      file.path(
        cartella_analisi,
        paste0("Probabilita_GTCF_", copula, ".csv")
      ),
      check.names = FALSE
    ),
    risultati = read.csv(
      file.path(
        cartella_analisi,
        paste0("Risultati_Finali_", copula, ".csv")
      ),
      check.names = FALSE
    ),
    montecarlo = read.csv(
      file.path(
        cartella_montecarlo,
        paste0("MonteCarlo_Risultati_", copula, ".csv")
      ),
      check.names = FALSE
    )
  )
}

input_copule <- setNames(
  lapply(ordine_copule, importa_copula),
  ordine_copule
)

probabilita <- lapply(input_copule, `[[`, "probabilita")
risultati <- lapply(input_copule, `[[`, "risultati")
montecarlo <- lapply(input_copule, `[[`, "montecarlo")



# 3. Aggregazione dei risultati ------------------------------------------------

# 3.1 Analisi singole a Kendall tau = 1/3
# Aggiunge il nome della copula e rende omogenea la colonna
# che contiene il rispettivo parametro.

prepara_risultati <- function(dati, copula) {

  nome_parametro <- parametri_copule[[copula]]
  valore_parametro <- dati[[nome_parametro]]

  dati |>
    mutate(
      copula = .env$copula,
      nome_parametro = .env$nome_parametro,
      valore_parametro = .env$valore_parametro
    ) |>
    select(
      copula, nome_parametro, valore_parametro, coda, numerosita,
      tau_kendall, soglia_minima, soglia_massima, numero_soglie,
      eventi_soglia_minima, beta_empirico, alpha_empirico,
      livello_gtcf_empirico, r_quadro_descrittivo,
      beta_teorico_finestra, alpha_teorico_finestra,
      livello_gtcf_teorico_finestra, errore_alpha,
      errore_livello_gtcf, beta_asintotico, alpha_asintotico,
      livello_gtcf_asintotico
    )
}

confronto_analisi_singole <- bind_rows(Map(
  prepara_risultati,
  risultati,
  names(risultati)
)) |>
  mutate(
    copula = factor(copula, levels = ordine_copule),
    coda = factor(coda, levels = ordine_code),
    tau_kendall = round(tau_kendall, 12)
  ) |>
  arrange(copula, coda)


# 3.2 Risultati Monte Carlo standardizzati
# Ricava Kendall tau dal parametro della copula quando necessario
# e seleziona le misure comuni alle tre famiglie.

prepara_montecarlo <- function(dati, copula) {

  tau_kendall <- switch(
    copula,
    "Gaussiana" = dati$tau_kendall,
    "Clayton" = dati$theta / (dati$theta + 2),
    "Gumbel" = 1 - 1 / dati$theta
  )

  dati |>
    transmute(
      copula = .env$copula,
      tau_kendall = round(.env$tau_kendall, 12),
      numerosita, coda, alpha_teorico_finestra, alpha_asintotico,
      livello_gtcf_teorico_finestra, livello_gtcf_asintotico,
      media_alpha, bias_alpha, rmse_alpha, media_livello_gtcf,
      bias_livello_gtcf, rmse_livello_gtcf,
      rmse_relativo_alpha_percentuale,
      rmse_relativo_livello_gtcf_percentuale,
      percentuale_stime_valide =
        percentuale_stime_numericamente_valide
    )
}

confronto_montecarlo_completo <- bind_rows(Map(
  prepara_montecarlo,
  montecarlo,
  names(montecarlo)
)) |>
  mutate(
    copula = factor(copula, levels = ordine_copule),
    coda = factor(coda, levels = ordine_code)
  ) |>
  arrange(tau_kendall, copula, coda, numerosita)

validazione_montecarlo_comparata <- confronto_montecarlo_completo |>
  select(
    copula, tau_kendall, numerosita, coda,
    alpha_teorico_finestra, media_alpha, bias_alpha, rmse_alpha,
    rmse_relativo_alpha_percentuale,
    livello_gtcf_teorico_finestra, media_livello_gtcf,
    bias_livello_gtcf, rmse_livello_gtcf,
    rmse_relativo_livello_gtcf_percentuale,
    percentuale_stime_valide
  )


# 3.3 Firma teorica delle code
# Raccoglie alpha e livello della GTCF per entrambe le code
# a parità di famiglia e Kendall tau.

benchmark_code <- confronto_montecarlo_completo |>
  select(
    copula, tau_kendall, coda, alpha_teorico_finestra,
    alpha_asintotico, livello_gtcf_teorico_finestra,
    livello_gtcf_asintotico
  ) |>
  distinct()

firma_lower <- benchmark_code |>
  filter(coda == "Lower") |>
  select(-coda) |>
  rename(
    alpha_lower_finestra = alpha_teorico_finestra,
    alpha_lower_asintotico = alpha_asintotico,
    c_lower_finestra = livello_gtcf_teorico_finestra,
    c_lower_asintotico = livello_gtcf_asintotico
  )

firma_upper <- benchmark_code |>
  filter(coda == "Upper") |>
  select(-coda) |>
  rename(
    alpha_upper_finestra = alpha_teorico_finestra,
    alpha_upper_asintotico = alpha_asintotico,
    c_upper_finestra = livello_gtcf_teorico_finestra,
    c_upper_asintotico = livello_gtcf_asintotico
  )

firma_coda_gtcf <- firma_lower |>
  inner_join(firma_upper, by = c("copula", "tau_kendall")) |>
  arrange(tau_kendall, copula)


# 3.4 Curve GTCF comparate a Kendall tau = 1/3
# Trasforma le colonne lower e upper in una struttura comune.

prepara_curve <- function(dati, copula) {

  bind_rows(
    dati |>
      transmute(
        copula = .env$copula, tau_kendall = round(1 / 3, 12),
        soglia_x, coda = "Lower", conteggi = conteggi_lower,
        probabilita_empirica = probabilita_lower,
        probabilita_teorica = probabilita_teorica_lower,
        gtcf_empirica = gtcf_empirica_lower,
        gtcf_teorica_finestra = gtcf_teorica_finestra_lower,
        gtcf_teorica_asintotica = gtcf_teorica_asintotica_lower
      ),
    dati |>
      transmute(
        copula = .env$copula, tau_kendall = round(1 / 3, 12),
        soglia_x, coda = "Upper", conteggi = conteggi_upper,
        probabilita_empirica = probabilita_upper,
        probabilita_teorica = probabilita_teorica_upper,
        gtcf_empirica = gtcf_empirica_upper,
        gtcf_teorica_finestra = gtcf_teorica_finestra_upper,
        gtcf_teorica_asintotica = gtcf_teorica_asintotica_upper
      )
  )
}

curve_gtcf_comparate <- bind_rows(Map(
  prepara_curve,
  probabilita,
  names(probabilita)
)) |>
  mutate(
    copula = factor(copula, levels = ordine_copule),
    coda = factor(coda, levels = ordine_code)
  ) |>
  left_join(
    confronto_analisi_singole |>
      select(
        copula, coda, soglia_minima_stima = soglia_minima,
        soglia_massima_stima = soglia_massima, alpha_empirico,
        alpha_teorico_finestra, alpha_asintotico,
        livello_gtcf_empirico, livello_gtcf_teorico_finestra,
        livello_gtcf_asintotico
      ),
    by = c("copula", "coda")
  ) |>
  arrange(copula, coda, soglia_x)



# 4. Salvataggio delle quattro tabelle -----------------------------------------

# Il numero di righe Monte Carlo si adatta automaticamente alle
# numerosità comuni importate dai tre CSV.

righe_montecarlo_attese <-
  length(ordine_copule) *
  length(tau_attesi) *
  length(ordine_code) *
  length(numerosita_confronto)

stopifnot(
  nrow(validazione_montecarlo_comparata) == righe_montecarlo_attese,
  nrow(distinct(
    validazione_montecarlo_comparata,
    copula, tau_kendall, numerosita, coda
  )) == righe_montecarlo_attese,
  identical(sort(unique(firma_coda_gtcf$tau_kendall)), tau_attesi)
)

write.csv(
  confronto_analisi_singole,
  file = file.path(
    cartella_output,
    "Confronto_Analisi_Singole_tau_1_3.csv"
  ),
  row.names = FALSE
)

write.csv(
  firma_coda_gtcf,
  file = file.path(cartella_output, "Firma_Coda_GTCF.csv"),
  row.names = FALSE
)

write.csv(
  validazione_montecarlo_comparata,
  file = file.path(
    cartella_output,
    "Validazione_MonteCarlo_Comparata.csv"
  ),
  row.names = FALSE
)

write.csv(
  curve_gtcf_comparate,
  file = file.path(
    cartella_output,
    "Curve_GTCF_Comparate_tau_1_3.csv"
  ),
  row.names = FALSE
)



# 5. Sistema grafico comune ----------------------------------------------------

# 5.1 Palette
# Riprende i toni usati negli script delle singole copule.

col_riferimento <- "#9AA1A8"
col_finestra <- "#ECEFF1"
col_griglia <- "#E6E6E6"
col_testo <- "#222222"

palette_copule <- c(
  "Gaussiana" = "#67A9CF",
  "Clayton" = "#EF8A62",
  "Gumbel" = "#66C2A5"
)

forme_copule <- c(
  "Gaussiana" = 16,
  "Clayton" = 17,
  "Gumbel" = 15
)


# 5.2 Tema grafico comune
# Mantiene titolo allineato, griglia orizzontale, facet neutri
# e legenda in basso come nelle analisi singole.

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
      axis.text = element_text(color = col_testo),
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



# 6. Grafici ------------------------------------------------------------------

# 6.1 Mappa della firma di coda

dati_mappa <- firma_coda_gtcf |>
  transmute(
    copula, tau_kendall, alpha_lower_asintotico,
    alpha_upper_asintotico,
    tau_label = case_when(
      tau_kendall == 0.20 ~ "τ = 0.20",
      tau_kendall == round(1 / 3, 12) ~ "τ = 1/3",
      tau_kendall == 0.50 ~ "τ = 0.50"
    ),
    tau_label = factor(
      tau_label,
      levels = c("τ = 0.20", "τ = 1/3", "τ = 0.50")
    )
  )

etichette_mappa <- dati_mappa |>
  filter(copula != "Clayton") |>
  mutate(
    label_x = alpha_lower_asintotico +
      if_else(copula == "Gaussiana", 0.008, 0),
    label_y = alpha_upper_asintotico - 0.020
  )

etichetta_clayton <- dati_mappa |>
  filter(copula == "Clayton") |>
  slice(1) |>
  mutate(
    tau_label = "τ = 0.20, 1/3, 0.50",
    label_x = 0.985,
    label_y = 0.525
  )

grafico_mappa_firma_coda <-
  ggplot(
    data = dati_mappa,
    aes(
      x = alpha_lower_asintotico,
      y = alpha_upper_asintotico,
      color = copula,
      group = copula
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = col_riferimento,
      linetype = "22",
      linewidth = 0.70
    ) +
    geom_path(linewidth = 0.95) +
    geom_point(size = 3.3) +
    geom_text(
      data = etichette_mappa,
      aes(x = label_x, y = label_y, label = tau_label, color = copula),
      inherit.aes = FALSE,
      size = 3.25,
      show.legend = FALSE
    ) +
    geom_text(
      data = etichetta_clayton,
      aes(x = label_x, y = label_y, label = tau_label, color = copula),
      inherit.aes = FALSE,
      hjust = 1,
      size = 3.25,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = 0.59,
      y = 0.90,
      label = "Coda superiore più persistente",
      color = "#5A5A5A",
      size = 3.25
    ) +
    annotate(
      "text",
      x = 0.90,
      y = 0.59,
      label = "Coda inferiore più persistente",
      color = "#5A5A5A",
      size = 3.25
    ) +
    scale_color_manual(values = palette_copule) +
    scale_x_continuous(
      breaks = seq(0.5, 1, 0.1),
      limits = c(0.48, 1.02)
    ) +
    scale_y_continuous(
      breaks = seq(0.5, 1, 0.1),
      limits = c(0.48, 1.02)
    ) +
    coord_fixed(clip = "off") +
    labs(
      title = "Mappa della firma di coda",
      subtitle = "Esponenti α asintotici a parità di Kendall τ.",
      x = expression(alpha[L * "," * infinity]),
      y = expression(alpha[U * "," * infinity]),
      caption = paste0(
        "La diagonale indica simmetria tra le code. ",
        "Per Clayton i tre livelli di τ coincidono in (1; 0,5)."
      )
    ) +
    tema_tesi()


# 6.2 Alpha di finestra e alpha asintotico

dati_alpha <- bind_rows(
  firma_coda_gtcf |>
    transmute(
      copula,
      tau_kendall,
      coda = livelli_coda[[1]],
      finestra = alpha_lower_finestra,
      asintotico = alpha_lower_asintotico
    ),
  firma_coda_gtcf |>
    transmute(
      copula,
      tau_kendall,
      coda = livelli_coda[[2]],
      finestra = alpha_upper_finestra,
      asintotico = alpha_upper_asintotico
    )
) |>
  pivot_longer(
    c(finestra, asintotico),
    names_to = "riferimento_alpha",
    values_to = "alpha"
  ) |>
  mutate(
    copula = factor(copula, levels = ordine_copule),
    coda = factor(coda, levels = livelli_coda),
    riferimento_alpha = recode(
      riferimento_alpha,
      finestra = "Alpha di finestra",
      asintotico = "Alpha asintotico"
    ),
    riferimento_alpha = factor(
      riferimento_alpha,
      levels = c("Alpha di finestra", "Alpha asintotico")
    )
  )

linee_alpha <- c(
  "Alpha di finestra" = "solid",
  "Alpha asintotico" = "22"
)

punti_alpha <- c(
  "Alpha di finestra" = 16,
  "Alpha asintotico" = 1
)

grafico_confronto_alpha <-
  ggplot(
    data = dati_alpha,
    aes(
      x = tau_kendall,
      y = alpha,
      color = copula,
      linetype = riferimento_alpha,
      shape = riferimento_alpha,
      group = interaction(copula, riferimento_alpha)
    )
  ) +
    geom_line(linewidth = 0.90) +
    geom_point(size = 2.75) +
    facet_grid(coda ~ copula) +
    scale_color_manual(
      values = palette_copule,
      guide = "none"
    ) +
    scale_linetype_manual(
      values = linee_alpha,
      breaks = names(linee_alpha)
    ) +
    scale_shape_manual(
      values = punti_alpha,
      breaks = names(punti_alpha)
    ) +
    scale_x_continuous(
      breaks = tau_attesi,
      labels = c("0,20", "1/3", "0,50")
    ) +
    scale_y_continuous(
      breaks = seq(0.5, 1, 0.1),
      limits = c(0.48, 1.02),
      labels = scales::label_number(
        accuracy = 0.1,
        decimal.mark = ",",
        trim = TRUE
      )
    ) +
    labs(
      title = "Alpha di finestra e alpha asintotico",
      subtitle = "Confronto per copula, coda e Kendall τ.",
      x = expression(tau),
      y = expression(alpha),
      caption = paste0(
        "La distanza tra le linee misura l'effetto della finestra finita ",
        "rispetto al comportamento limite."
      )
    ) +
    tema_tesi()


# 6.3 Curve GTCF teoriche comparate

dati_grafico_gtcf <- curve_gtcf_comparate |>
  select(
    copula, coda, soglia_x,
    gtcf_teorica_finestra, gtcf_teorica_asintotica
  ) |>
  pivot_longer(
    c(gtcf_teorica_finestra, gtcf_teorica_asintotica),
    names_to = "normalizzazione",
    values_to = "gtcf_teorica"
  ) |>
  mutate(
    coda = recode(
      as.character(coda),
      "Lower" = livelli_coda[[1]],
      "Upper" = livelli_coda[[2]]
    ),
    coda = factor(coda, levels = livelli_coda),
    normalizzazione = recode(
      normalizzazione,
      "gtcf_teorica_finestra" = "Alpha di finestra",
      "gtcf_teorica_asintotica" = "Alpha asintotico"
    ),
    normalizzazione = factor(
      normalizzazione,
      levels = c("Alpha di finestra", "Alpha asintotico")
    )
  )

finestre_gtcf <- curve_gtcf_comparate |>
  distinct(
    copula, coda, soglia_minima_stima, soglia_massima_stima
  ) |>
  mutate(
    coda = recode(
      as.character(coda),
      "Lower" = livelli_coda[[1]],
      "Upper" = livelli_coda[[2]]
    ),
    coda = factor(coda, levels = livelli_coda)
  )

linee_normalizzazione <- c(
  "Alpha di finestra" = "solid",
  "Alpha asintotico" = "22"
)

grafico_curve_gtcf <-
  ggplot(
    data = dati_grafico_gtcf,
    aes(
      x = soglia_x,
      y = gtcf_teorica,
      color = copula,
      linetype = normalizzazione,
      group = interaction(copula, normalizzazione)
    )
  ) +
    geom_rect(
      data = finestre_gtcf,
      aes(
        xmin = soglia_minima_stima,
        xmax = soglia_massima_stima,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = col_finestra,
      alpha = 0.88
    ) +
    geom_vline(
      data = finestre_gtcf,
      aes(xintercept = soglia_minima_stima),
      inherit.aes = FALSE,
      color = col_riferimento,
      linetype = "22",
      linewidth = 0.55
    ) +
    geom_vline(
      data = finestre_gtcf,
      aes(xintercept = soglia_massima_stima),
      inherit.aes = FALSE,
      color = col_riferimento,
      linetype = "22",
      linewidth = 0.55
    ) +
    geom_line(linewidth = 0.95) +
    facet_grid(coda ~ copula) +
    scale_color_manual(
      values = palette_copule,
      guide = "none"
    ) +
    scale_linetype_manual(
      values = linee_normalizzazione,
      breaks = names(linee_normalizzazione)
    ) +
    scale_x_continuous(
      limits = c(0, 0.50),
      breaks = c(0.01, 0.10, 0.20, 0.30, 0.40, 0.50),
      labels = scales::label_number(
        accuracy = 0.01,
        decimal.mark = ",",
        trim = TRUE
      ),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.1,
        decimal.mark = ",",
        trim = TRUE
      ),
      expand = expansion(mult = c(0.06, 0.08))
    ) +
    labs(
      title = "Curve GTCF teoriche comparate",
      subtitle = "Normalizzazione di finestra e asintotica, Kendall τ = 1/3.",
      x = "Soglia x",
      y = "GTCF teorica",
      caption = paste0(
        "Le scale verticali sono comuni a tutti i pannelli. ",
        "L'area grigia e le linee tratteggiate delimitano la finestra ",
        "effettiva di stima di α."
      )
    ) +
    tema_tesi() +
    theme(
      panel.spacing.x = grid::unit(2.2, "lines")
    )


# 6.4 RMSE relativo Monte Carlo

riduzioni_rmse <- validazione_montecarlo_comparata |>
  filter(numerosita %in% c(20000, 100000)) |>
  select(
    copula, tau_kendall, coda, numerosita,
    alpha = rmse_relativo_alpha_percentuale,
    livello_gtcf = rmse_relativo_livello_gtcf_percentuale
  ) |>
  pivot_longer(
    c(alpha, livello_gtcf),
    names_to = "parametro",
    values_to = "rmse_percentuale"
  ) |>
  pivot_wider(
    names_from = numerosita,
    values_from = rmse_percentuale,
    names_prefix = "n_"
  ) |>
  mutate(
    riduzione_percentuale =
      100 * (1 - n_100000 / n_20000)
  ) |>
  group_by(parametro) |>
  summarise(
    riduzione_mediana = median(riduzione_percentuale),
    .groups = "drop"
  )

riduzione_mediana_alpha <- riduzioni_rmse |>
  filter(parametro == "alpha") |>
  pull(riduzione_mediana)

riduzione_mediana_livello <- riduzioni_rmse |>
  filter(parametro == "livello_gtcf") |>
  pull(riduzione_mediana)

didascalia_riduzione_rmse <- paste0(
  "Da n=20.000 a n=100.000, l'RMSE relativo mediano diminuisce del ",
  formatC(
    riduzione_mediana_alpha,
    format = "f",
    digits = 1,
    decimal.mark = ","
  ),
  "% per α e del ",
  formatC(
    riduzione_mediana_livello,
    format = "f",
    digits = 1,
    decimal.mark = ","
  ),
  "% per il livello GTCF."
)

dati_rmse <- validazione_montecarlo_comparata |>
  select(
    copula, tau_kendall, numerosita, coda,
    rmse_relativo_alpha_percentuale,
    rmse_relativo_livello_gtcf_percentuale
  ) |>
  pivot_longer(
    c(
      rmse_relativo_alpha_percentuale,
      rmse_relativo_livello_gtcf_percentuale
    ),
    names_to = "parametro",
    values_to = "rmse_percentuale"
  ) |>
  mutate(
    numerosita = factor(
      numerosita,
      levels = numerosita_confronto,
      labels = format(
        numerosita_confronto,
        big.mark = ".",
        decimal.mark = ",",
        scientific = FALSE,
        trim = TRUE
      )
    ),
    coda = recode(
      as.character(coda),
      "Lower" = livelli_coda[[1]],
      "Upper" = livelli_coda[[2]]
    ),
    coda = factor(coda, levels = livelli_coda),
    parametro = recode(
      parametro,
      "rmse_relativo_alpha_percentuale" = "Esponente α",
      "rmse_relativo_livello_gtcf_percentuale" = "Livello GTCF"
    ),
    parametro = factor(
      parametro,
      levels = c("Esponente α", "Livello GTCF")
    ),
    pannello_rmse = paste(
      recode(
        as.character(parametro),
        "Esponente α" = "α",
        "Livello GTCF" = "GTCF"
      ),
      recode(
        as.character(coda),
        "Coda inferiore" = "inferiore",
        "Coda superiore" = "superiore"
      ),
      sep = " · "
    ),
    pannello_rmse = factor(
      pannello_rmse,
      levels = c(
        "α · inferiore",
        "α · superiore",
        "GTCF · inferiore",
        "GTCF · superiore"
      )
    ),
    tau_label = case_when(
      tau_kendall == 0.20 ~ "τ = 0,20",
      tau_kendall == round(1 / 3, 12) ~ "τ = 1/3",
      tau_kendall == 0.50 ~ "τ = 0,50"
    ),
    tau_label = factor(
      tau_label,
      levels = c("τ = 0,20", "τ = 1/3", "τ = 0,50")
    )
  )

grafico_rmse_montecarlo <-
  ggplot(
    data = dati_rmse,
    aes(
      x = numerosita,
      y = rmse_percentuale,
      color = copula,
      shape = copula,
      group = copula
    )
  ) +
    geom_blank(aes(y = 0)) +
    geom_line(linewidth = 0.75, alpha = 0.88) +
    geom_point(size = 2.65) +
    facet_grid(
      pannello_rmse ~ tau_label,
      scales = "free_y"
    ) +
    scale_color_manual(
      values = palette_copule,
      breaks = names(palette_copule)
    ) +
    scale_shape_manual(
      values = forme_copule,
      breaks = names(forme_copule)
    ) +
    scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.1,
        decimal.mark = ",",
        suffix = "%",
        trim = TRUE
      ),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(
      title = "RMSE relativo di α e del livello GTCF",
      subtitle = paste0(
        "Validazione Monte Carlo per copula, coda, Kendall τ ",
        "e numerosità campionaria."
      ),
      x = "Numerosità campionaria",
      y = "RMSE relativo",
      caption = paste0(
        didascalia_riduzione_rmse,
        "\n",
        "Le scale verticali sono separate per combinazione di parametro ",
        "e coda per rendere leggibili ordini di grandezza differenti."
      )
    ) +
    tema_tesi()



# 7. Salvataggio dei grafici ---------------------------------------------------

salva_grafico <- function(grafico,
                          nome_file,
                          larghezza,
                          altezza) {

  ggsave(
    filename = file.path(cartella_output, nome_file),
    plot = grafico,
    width = larghezza,
    height = altezza,
    units = "in",
    dpi = 320,
    bg = "white"
  )
}

salva_grafico(
  grafico_mappa_firma_coda,
  "Mappa_Firma_Coda_Alpha.png",
  8.5, 7.5
)

salva_grafico(
  grafico_confronto_alpha,
  "Confronto_Alpha_Finestra_Asintotico.png",
  12, 6.5
)

salva_grafico(
  grafico_curve_gtcf,
  "Curve_GTCF_Comparate_tau_1_3.png",
  12, 6.5
)

salva_grafico(
  grafico_rmse_montecarlo,
  "RMSE_MonteCarlo_Comparato.png",
  12, 9
)



# 8. Fine ----------------------------------------------------------

cat(
  "\nEsecuzione completata: quattro CSV e quattro grafici salvati.\n"
)
