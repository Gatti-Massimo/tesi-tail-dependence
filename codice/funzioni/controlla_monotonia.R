# Verifica se una curva è non decrescente.

controlla_monotonia <- function(
    curva,
    tolleranza = 1e-10
) {

  if (length(curva) < 2L) {
    stop("La curva deve contenere almeno due valori.")
  }

  if (tolleranza < 0) {
    stop("La tolleranza deve essere non negativa.")
  }

  all(diff(curva) >= -tolleranza)
}
