# Stima i parametri della curva di potenza
# minimizzando la somma degli errori quadratici.

stima_curva_potenza <- function(
    x,
    tcf_normalizzata
) {

  optim(
    par = c(
      a = 0.05,
      gamma = 1
    ),

    fn = function(parametri) {

      residui <- tcf_normalizzata -
        curva_potenza(
          x = x,
          a = parametri["a"],
          gamma = parametri["gamma"]
        )

      sum(residui^2)
    },

    method = "L-BFGS-B",

    lower = c(
      a = 0,
      gamma = 0.001
    ),

    upper = c(
      a = 0.999,
      gamma = 20
    )
  )
}
