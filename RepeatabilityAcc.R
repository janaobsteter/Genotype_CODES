vg = 1
vpe = 0.5
vy = 4
ve = 1.5

t = (vg + vpe)/vy
t

n = 30
h2 = vg / vy
b = (n * h2) / (1 + (n - 1)*t)



vg = 1
vpe = 0.0

ve = 1.5



t = (vg + vpe)/vy
t

n = 3
h2 = 0.25
t = h2
b = (n * h2) / (1 + (n - 1)*t)


r2 <- function(h2, N, M) { (h2 * N) / (h2*N + M)}

r2(0.63, 10000, 1000)
