import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

record = list(SeqIO.parse("HiSeq_accuracy.fa", "fasta"))

## Populate the list
lst = []
length = len(record)
for i in range(length):
    lst.append(str(record[i].seq))

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

# Save the jaccard_mat to a file
np.save('jaccard_mat.npy', jaccard_mat)

# Load the saved matrix from file
jaccard_mat = np.load('jaccard_mat.npy')

fig, ax = plt.subplots()
heatmap = ax.pcolor(jaccard_mat, cmap=plt.cm.Blues)
plt.colorbar(heatmap)

# Show values in each cell
for i in range(length):
    for j in range(length):
        value = "{:.2f}".format(jaccard_mat[i, j])

plt.tight_layout()
plt.savefig("heatmap_HiSeq.png")
