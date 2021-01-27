members = ["Ana", "Joze", "miha", "Sonja", "Boris", "Eva"]

{member: len(member) for member in members}

double = {(member + str(number)): number for member in members for number in range(4)}
print(double)


[x * x2 for x in [1,2,3] for x2 in [3,4,5]]
