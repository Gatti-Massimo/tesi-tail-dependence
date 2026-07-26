# Calcola la GTCF:
# G_alpha(x) = P(x) / x^(1 / alpha)

calcola_gtcf <- function(
    probabilita_congiunta,
    soglia_x,
    alpha
) {

  if (length(probabilita_congiunta) != length(soglia_x)) {
    stop("Probabilità e soglie devono avere la stessa lunghezza.")
  }

  if (any(soglia_x <= 0)) {
    stop("Le soglie devono essere strettamente positive.")
  }

  if (length(alpha) != 1L || alpha < 0.5 || alpha > 1) {
    stop("Alpha deve essere compreso tra 0.5 e 1.")
  }

  probabilita_congiunta / soglia_x^(1 / alpha)
}
