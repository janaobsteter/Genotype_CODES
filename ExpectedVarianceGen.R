
V(t) = V(t-1)*(1-1/(2Ne)) 
V(t) = V(0) (1- 1/(2Ne))^t

V0 = 1
t = 20
Ne=172

#Pricakovana SD v t=200
sd <- sqrt(V0 * (1 - 1/(2*Ne)) ^ t)

#Ne v v alternativnem scenariju
NEa = (1-0.5)*Ne

# Pricakovana SD v t=20 v alternativnem scenariju
sdA <- sqrt(V0 * (1 - 1/(2*NEa)) ^ t)

#Relativna sprememba
diff <- (1 - sd/sdA) * 100


print("Pricakovana SD v t=200")
print(sd)
print("Pricakovana SD v t=200 v alternativnem scenariju")
print(sdA)
print("Relativna razlika")
print(diff)