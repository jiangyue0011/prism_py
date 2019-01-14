import sys

def load_clusters(cluster_file):
    with open(cluster_file, 'r') as f:
        clusters = []
        for line in f:
            tmp = line.split()
            ref = tmp[0]
            sv_start = int(tmp[1])
            sv_end = int(tmp[2])
            clusters.append((ref, sv_start, sv_end, "DEL"))
    return clusters
