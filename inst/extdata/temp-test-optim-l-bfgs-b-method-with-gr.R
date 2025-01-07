start = c(-18,18)
true_params <- c(-5,5)

foo <- function(x, y, a, b) {8*2*(x - a)^2 + (y - b)^2}
bar <- function(x, y, a, b) {c(8*2*(x - a), 2*(y - b))}

foo_wrapped <- function(args,a,b){foo(args[1], args[2], a, b)}
bar_wrapped <- function(args,a,b){bar(args[1], args[2], a, b)}

xseq = seq(-20,20,by=0.1)
yseq = seq(-20,20,by=0.1)
grid <- expand.grid(x = xseq, y = yseq)


image(x = xseq, y = yseq,
      z = matrix(
        data = foo(x = grid$x,
                   y = grid$y,
                   a = true_params[1],
                   b = true_params[2]),
        nrow = length(xseq),
        byrow = TRUE))
points(x = -true_params[1], y = -true_params[2], col = 'black')

t0 <- Sys.time()
optim_no_gr <- optim(par = start,
      fn = foo_wrapped,
      a = true_params[1],
      b = true_params[2],
      method = "L-BFGS-B",
      lower = c(-20,-20),
      upper = c(20,20))
print(Sys.time() - t0)

t0 <- Sys.time()
optim_gr <- optim(par = start,
                  fn = foo_wrapped,
                  gr = bar_wrapped,
                  a = true_params[1],
                  b = true_params[2],
                  method = "L-BFGS-B",
                  lower = c(-10,-10),
                  upper = c(10,10))
print(Sys.time() - t0)
