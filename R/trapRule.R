trapRule <- function(y, x = 1:length(y), maxT = max(x)) {
  y <- y[x <= maxT]
  x <- x[x <= maxT]
  
  y <- y - y[1]
  dx <- diff(x)
  dy <- diff(y)
  miny <- pmin(head(y, -1), tail(y, -1))

  sum(dx*(miny + dy/2)/max(x))
}
