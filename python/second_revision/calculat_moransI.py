"""
Calculate Moran's I for cline and patchy maps in 10x10 grid
tianlin.duan42@gmail.com
2026.09.21
"""
############################# modules #########################################
#import scanpy
import numpy as np
from libpysal.weights import lat2W, DistanceBand
from esda.moran import Moran

############################# program #########################################
patchy_map = np.loadtxt("/home/anadem/github/arg-for-gea/map/random_map1.csv", delimiter=",")
cline_map = np.tile(np.arange(0.05, 1.0, 0.1), (10, 1))

# Create a matrix of weights, using Contiguity-Based Weights: Queen
w = lat2W(nrows=10, ncols=10, rook=False)
#w.transform="r"

#Crate the pysal Moran object and calculate Moran's I
# Patchy: No row normalization: -0.03755
patchy_stats = Moran(patchy_map, w)
patchy_I = patchy_stats.I

#Cline: No row normalization: 0.9331
cline_stats = Moran(cline_map, w)
cline_I = cline_stats.I


# # Create a matrix of weights, using Distance-Based Weights: Not used
# Method 1:
# coords = np.array([(x, y)
#                    for y in range(10)
#                    for x in range(10)])
# w = DistanceBand(coords, threshold=5, binary=False, alpha=-1)
# w.transform = "r"
#
# # Method 2: 8 nearest cells
# w_knn = KNN.from_array(
#     [(x, y_) for y_ in range(10) for x in range(10)],
#     k=8
# )
# w_knn.transform = "r"