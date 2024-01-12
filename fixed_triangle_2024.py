import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
from matplotlib.colors import LinearSegmentedColormap


heatmap_data_g = pd.read_csv("./fixjan2024_all_regions_sum_nPix_perk_green_channel_PaperData_thres50.csv",index_col=0)
heatmap_data_r = pd.read_csv("./fixjan2024_all_regions_sum_nPix_perk_red_channel_PaperData_thres50.csv",index_col=0)

heatmap_data_g = heatmap_data_g.T
heatmap_data_r = heatmap_data_r.T

zs_g = heatmap_data_g.values
zs_r = heatmap_data_r.values
M,N = zs_g.shape

print(M,N)

x = np.arange(M + 1)
y = np.arange(N + 1)
#gap_size = 0.1
#y = np.arange(N * (1 + gap_size) + 1)
xs, ys = np.meshgrid(x, y)

triangles1 = [(i + j*(M+1), i+1 + j*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
triangles2 = [(i+1 + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
triang1 = Triangulation(xs.ravel()-0.5, ys.ravel()-0.5, triangles1)
triang2 = Triangulation(xs.ravel()-0.5, ys.ravel()-0.5, triangles2)


# Flatten the 2D arrays
#zs_g_1d = zs_g.ravel()
zs_g_1d = zs_g.flatten(order='F')
zs_r_1d = zs_r.flatten(order='F')

cmap_pinks = LinearSegmentedColormap.from_list('more_magenta', ['#ffffff', '#F702F7'], N=256)
cmap_greens = LinearSegmentedColormap.from_list('more_magenta', ['#ffffff', '#026a2c'], N=256)

fig, ax = plt.subplots(figsize=(9,9))

img1 = plt.tripcolor(triang1, zs_g_1d, cmap=plt.get_cmap(cmap_greens), vmax=0.15)
#img1 = plt.tripcolor(triang1, zs_g_1d, cmap=plt.get_cmap('Greens', 1000), vmax=10000)
img2 = plt.tripcolor(triang2, zs_r_1d, cmap=plt.get_cmap(cmap_pinks), vmax=0.15)
#img2 = plt.tripcolor(triang2, zs_r_1d, cmap=plt.get_cmap(cmap_pinks, 1000), vmax=10000)

print(heatmap_data_g.index)
#ax.set_xticks(x)
ax.set_xticklabels(heatmap_data_g.index, rotation=90)
#ax.set_yticks(y)
ax.set_yticklabels(heatmap_data_g.columns, style='italic')

fig.colorbar(img2)
fig.colorbar(img1)
plt.xlim(x[0]-0.5, x[-1]-0.5)
plt.ylim(y[0]-0.5, y[-1]-0.5)
plt.xticks(x[:-1], rotation=90)
plt.yticks(y[:-1]) 
plt.tight_layout()
plt.savefig("autism_pErk_Jan2024_trianglepErk_reordered_percents.pdf")
plt.show()
