import pandas as pd
import matplotlib.pyplot as plt


x = pd.read_csv("output.txt", header=None)

plt.plot(x[0], x[2], 'o', color='black')
plt.show()


        