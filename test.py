import math
import numpy as np

nr_splits = 10 #If sufficient data is generated, increase nr_splits as it provides more accurate cutoff points
pval_array = [0.1532, 0.4532, 0.4532, 2.44532, 0.16532, 0.718532]
xaxis = []
yaxis = []
sum = 0
max_val = math.ceil(max(pval_array))
min_val = math.floor(min(pval_array))
step = (max_val - min_val) / nr_splits
for i in np.arange(min_val, max_val, step):
    xaxis.append(i)
    y = len([j for j in pval_array if j < i])
    diff = y - sum
    sum = y
    yaxis.append(diff / len(pval_array))

print(xaxis)
print(yaxis)


file = open("Matlab/xaxis.txt", "w")
for i in xaxis:
    file.write(str(i) + " ")
file.close()
variable_name = "correct"
file = open("Matlab/yaxis.txt", "w")
for i in xaxis:
    file.write(str(i) + " ")
file.close()
