import numpy as np
from scipy import stats
from sys import argv

with open(argv[1], 'r') as fin:
    depth = [int(line) for line in fin.readlines()]
    fin.close()

print(f"5th percentile: {np.percentile(depth, 5)}")
print(f"95th percentile: {np.percentile(depth, 95)}\n")

median_dp = np.median(depth)
mad_dp = stats.median_abs_deviation(depth)

print(f"Median: {round(median_dp, 2)}")
print(f"Median absolute deviation: {round(mad_dp, 2)}")

min_cutoff = median_dp - (5 * mad_dp)
max_cutoff = median_dp + (5 * mad_dp)

print(f"Minimum depth cutoff (MEDIAN - 5 * MAD): {round(min_cutoff)}")
print(f"Maximum depth cutoff (MEDIAN + 5 * MAD): {round(max_cutoff)}\n")
