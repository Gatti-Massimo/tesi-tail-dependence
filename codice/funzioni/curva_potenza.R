# Famiglia crescente di curve di potenza:
# f(x) = a + (1 - a)(2x)^gamma
#
# con:
# 0 <= a <= 1
# gamma > 0
# f(0.5) = 1

curva_potenza <- function(x, a, gamma) {
  a + (1 - a) * (2 * x)^gamma
}
