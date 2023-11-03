a <- matrix(1:3, 1:3)
b <- matrix(1:3, 1:3)

c <- cbind(a, b)

write.table(c, "C.txt")
