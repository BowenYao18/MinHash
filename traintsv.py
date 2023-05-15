import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
import pandas as pd
import random

# data_path = "/scratch1/zx22/bio/promoter-1k/train.tsv"
data_path = "/scratch1/zx22/bio/dnabert_ft/train.tsv"
df = pd.read_table(data_path)
df.dropna()

length = 100
# Populate the list
lst = []
arr = np.random.permutation(len(df))[:length]
for i in range(length):
    sequence = str(df.iloc[arr[i]].sequence).replace(" ", "")
    lst.append(sequence)

def jaccard_similarity(str1, str2, k):
    set1 = set()
    set2 = set()
    for i in range(len(str1) - k + 1):
        set1.add(str1[i : i + k])
    for j in range(len(str2) - k + 1):
        set2.add(str2[j : j + k])
    intersection = set1.intersection(set2)
    union = set1.union(set2)
    return len(intersection) / len(union)

def jaccard_matrix(lst, k):
    '''Return the jaccard distance matrix for all sequence. '''
    length = len(lst)
    jaccard_mat = np.zeros((length, length))
    for i in range(length):
        for j in range(i, length):
            sim = jaccard_similarity(lst[i], lst[j], k)
            jaccard_mat[i, j] = sim
            jaccard_mat[j, i] = sim
    return jaccard_mat

jaccard_mat = jaccard_matrix(lst, k=5)

# # Save the jaccard_mat to a file
# np.save('jaccard_mat2.npy', jaccard_mat)

# # Load the saved matrix from file
# jaccard_mat = np.load('jaccard_mat2.npy')

# # Save the jaccard_mat to a file
# np.save('jaccard_mat_dnabert_ft.npy', jaccard_mat)

# # Load the saved matrix from file
# jaccard_mat = np.load('jaccard_mat_dnabert_ft.npy')

# Save the jaccard_mat to a file
np.save('jaccard_mat_promoter.npy', jaccard_mat)

# Load the saved matrix from file
jaccard_mat = np.load('jaccard_mat_promoter.npy')

fig, ax = plt.subplots()
heatmap = ax.pcolor(jaccard_mat[:length, :length], cmap="bwr")
plt.colorbar(heatmap)

# Show values in each cell
for i in range(length):
    for j in range(length):
        value = "{:.2f}".format(jaccard_mat[i, j])

plt.tight_layout()
# plt.title("promoter-1k")
# plt.savefig("promoter_train_data.png")
plt.title("dnabert")
plt.savefig("dnabert_train_data.png")
