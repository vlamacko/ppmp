#!/usr/bin/env python3

""" prediction_analysis.py
Further analysis of the predictions. Currently separate from the package.
It may be integrate into ppmp later.
"""

import glob
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


CSV_PATH = './out/prediction/'
q = 0
single_count = 0
triplet_count = 0
single = 0
triplet = 0
total_single_std = 0
total_triplet_std = 0
single_error = 0
triplet_error = 0

data = {"x": [], "yre": [], "y": [], "ID": [], "std": [], "zeros": []}

for path in glob.glob(CSV_PATH + '*.csv'):
    file_name = os.path.basename(path)
    df = pd.read_csv(path)

    for i in range(len(df['rmsd mean'])):
        # Calculate absolute error
        if df['prediction mean'][i] == df['single mean'][i]:
            single += 1
            single_error += abs((df['prediction mean'][i] - df['rmsd mean'][i]))
        else:
            triplet += 1
            triplet_error += abs((df['prediction mean'][i] - df['rmsd mean'][i]))
        # Reject correctly predicted
        if df['prediction mean'][i] - df['prediction std'][i] \
                <= df['rmsd mean'][i] \
                <= df['prediction mean'][i] + df['prediction std'][i]:
            q += 1
        # Populate dictionary for the graph
        else:
            if df['prediction mean'][i] == df['single mean'][i]:
                single_count += 1
                total_single_std += df['prediction std'][i]
                data["ID"].append(1)
            else:
                triplet_count += 1
                total_triplet_std += df['prediction std'][i]
                data["ID"].append(2)
            data["x"].append(df['module'][i])
            data["y"].append(df['prediction mean'][i])
            data["yre"].append(df['rmsd mean'][i])
            data["std"].append(df["prediction std"][i])

x = np.linspace(0, 12, 13)
myx_ticks = data["x"]

fig = plt.figure(figsize=(17, 14))
plt.xticks(x, myx_ticks)
plt.title('Plot Demonstrating the Error of the Incorrectly Predicted Modules', fontsize=20)
plt.xlabel('Module', fontsize=20)
plt.ylabel('RMSD', fontsize=20)

for i in range(len(data["x"])):
    if data["ID"][i] == 1:
        plt.scatter(x[i], data["y"][i], c='r', marker='o', s=75)
        plt.errorbar(x[i], data["y"][i], yerr=data["std"][i], fmt='o', c='r', elinewidth=2)
    if data["ID"][i] == 2:
        plt.scatter(x[i], data["y"][i], c='b', marker='o', s=75)
        plt.errorbar(x[i], data["y"][i], yerr=data["std"][i], fmt='o', c='b', elinewidth=2)

plot = plt.scatter(x, data["yre"], marker='o', c='black', label='The average RMSD from simulation')
red_patch = mpatches.Patch(color='r', label='Predictions made from a module distribution')
blue_patch = mpatches.Patch(color='b', label='Predictions made from the triplet distribution')
data_patch = mpatches.Patch(color='black', label='The average RMSD from simulation')
plt.legend(handles=[red_patch, blue_patch, plot], prop={'size': 20})

plt.savefig("./out/prediction/analysis.pdf")
