library(reticulate)

# Load python shell
repl_python()

# Check python version
reticulate::py_config()

# Use another python version
use_python("/home/jana/bin/Selection1/bin/python2.7", required=TRUE)
repl_python()

# Import python libraries into R as an object
plt <- import('matplotlib.pyplot')
np <- import('numpy')

# Generate random number
np$random$seed(19680801L)

N <- 50L
#python: x = np.random.rand(50)
x <- np$random$rand(N) 
y <- np$random$rand(N) 
colours <- np$random$rand(N) 
area <- (30L * as.numeric(np$random$rand(N))) ** 2L

#python: plt.scatter(x, y, s=area, c=colours, alpha=0.5)
plt$scatter(x, y, s=area, c=colours, alpha=0.5)

plt$show()
