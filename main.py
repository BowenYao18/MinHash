from Bio import SeqIO
from sklearn.utils.murmurhash import murmurhash3_32
import sys
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

record = list(SeqIO.parse("HiSeq_accuracy.fa", "fasta"))

## Populate the list
lst = []
id = []
for i in range(len(record)):
        lst.append(str(record[i].seq))
        id.append(str(record[i].id))
# print(lst)
# print(id)

def hash_min(string, k, seed, hash_range):
        '''Calculate the min hash number. '''
        min = sys.maxsize
        for i in range(len(string) - k + 1):
                value = murmurhash3_32(string[i : i + k], seed, positive=True)
                if value < min:
                        min = value
        return min % hash_range

def build_matrix(list_seq, hash_range=20, time=20, k=5):
        '''Return the martix for all hashing process. '''
        matrix = [[0 for j in range(hash_range)] for i in range(time)]
        ## populate the matrix
        for seed in range(time):
                for seq in list_seq:
                        min = hash_min(seq, k, seed, hash_range)
                        matrix[seed][min] += 1
        return matrix

def find(matrix, string, hash_range=20, time=20, k=5):
        '''Find the number of collide with other string.'''
        count = 0
        for seed in range(time):
                min = hash_min(string, k, seed, hash_range)
                count += matrix[seed][min]
        count /= time
        return count

def prob_dist(matrix, list_seq, hash_range=20, time=20, k=5):
        '''Find the probability distribution of each string. '''
        output = []
        for seq in list_seq:
                output.append(1 / find(matrix, seq, hash_range, time, k))
        s = sum(output)
        return [x / s for x in output]

global_matrix = build_matrix(lst, hash_range=20, time=20, k=5)
#print(global_matrix)
probability = prob_dist(global_matrix, lst, hash_range=20, time=20, k=5)
#print(probability)


def random_unique(id, pick_num):
        random_elements = random.sample(id, pick_num)
        #print(random_elements)
        lst = []
        for i in random_elements:
                name = i.split('.')[0]
                if name not in lst:
                        lst.append(name)
        #print(lst)
        return len(lst)

def check_unique_species(id, pick_num, probability):
        choice = np.random.choice(id, pick_num, probability)
        #print(choice)
        unique = []
        for item in choice:
                name = item.split('.')[0]
                if name not in unique:
                        unique.append(name)
        #print(unique)
        return len(unique)

# print("Random selction result: ", random_unique(id, 10))
# print("Algorithm selction result: ", check_unique_species(id, 10, probability))

# def test(time, pick_num):
#         count = 0
#         for i in range(time):
#                 if (random_unique(id, pick_num) < check_unique_species(id, pick_num, probability)):
#                         count += 1
#         return count/time

# print(test(100, 5))


def compare(str1, str2, k, hash_range, time):
        same = 0
        for seed in range(time):
                min1 = hash_min(str1, k, seed, hash_range)
                min2 = hash_min(str2, k, seed, hash_range)
                # print("min1",min1,"min2",min2)
                if (min1 == min2):
                        same += 1
        return same / time

def jaccard_similarity(str1, str2, k):
        set1 = set()
        set2 = set()
        for i in range(len(str1) - k + 1):
                set1.add(str1[i : i + k])
        for j in range(len(str2) - k + 1):
                set2.add(str2[j : j + k])
        intersection = set1.intersection(set2)
        union = set1.union(set2)
        # print(union)
        # print(intersection)
        return len(intersection) / len(union)

# print(jaccard_similarity(dna1, dna2, 3))
# print(compare(dna1, dna2, 3, 1000000000, 10000))

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

fig, ax = plt.subplots()
heatmap = ax.pcolor(jaccard_mat, cmap=plt.cm.Blues)
plt.colorbar(heatmap)

# Show values in each cell
for i in range(len(id)):
    for j in range(len(id)):
        value = "{:.2f}".format(jaccard_mat[i, j])

plt.tight_layout()
plt.savefig("heatmap.png")
