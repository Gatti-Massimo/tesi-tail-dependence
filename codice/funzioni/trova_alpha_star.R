# Cerca il minimo valore di alpha per cui la GTCF
# risulta non decrescente.

trova_alpha_star <- function(
    probabilita_congiunta,
    soglia_x,
    griglia_alpha,
    tolleranza = 1e-10
) {

  alpha_ammissibile <- vapply(
    griglia_alpha,
    function(alpha) {

      gtcf_corrente <- calcola_gtcf(
        probabilita_congiunta = probabilita_congiunta,
        soglia_x = soglia_x,
        alpha = alpha
      )

      controlla_monotonia(
        curva = gtcf_corrente,
        tolleranza = tolleranza
      )
    },
    logical(1)
  )

  valori_ammissibili <- griglia_alpha[alpha_ammissibile]

  if (length(valori_ammissibili) == 0) {
    return(
      list(
        alpha_star = NA_real_,
        gtcf_alpha_star = rep(NA_real_, length(soglia_x)),
        alpha_ammissibili = numeric(0)
      )
    )
  }

  alpha_star <- min(valori_ammissibili)

  list(
    alpha_star = alpha_star,
    gtcf_alpha_star = calcola_gtcf(
      probabilita_congiunta = probabilita_congiunta,
      soglia_x = soglia_x,
      alpha = alpha_star
    ),
    alpha_ammissibili = valori_ammissibili
  )
}
