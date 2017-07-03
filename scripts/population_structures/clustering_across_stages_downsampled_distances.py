import numpy as np
import os
import argparse
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("--algorithm", "-a", default="MDS")
parser.add_argument("--chromosome", "-c", default=None, type=int)
args = parser.parse_args()

datasets = ["rings", "schizonts", "trophozoites"]
algorithm = args.algorithm

labels = None
X = None
for i, dataset in enumerate(datasets):
    if args.chromosome is None:
        filename = (
            "results/ay2013/%s_10000_raw_%s_distances_downsampled.npy" % (
                dataset, algorithm))
    else:
        filename = (
            "results/ay2013/%s_10000_raw_%s_distances_chr%02d.npy" % (
                dataset, algorithm, args.chromosome))
    x = np.load(filename)
    x = x[:1000]
    if labels is not None:
        X[i, :len(x)] = x
        labels = np.concatenate([labels, i * np.ones(len(x))])
    else:
        labels = i * np.ones(len(x))
        X = np.nan * np.zeros((len(datasets), 1000, x.shape[1]))
        X[i, :len(x)] = x

X = X.reshape(-1, x.shape[1])
# Remove missing labels
X = X[np.invert(np.all(np.isnan(X), axis=1))]
# Remove any bead with missing value
X = X[:, np.invert(np.isnan(X.sum(axis=0)))]

X /= (X).sum(axis=1)[:, np.newaxis] / 1000

colors = ['#FFB773', '#FF7C00', "r",  # '#A65100',
          '#04859D', 'skyblue',
          '#015666', 'green', 'darkgreen']
markers = [".", ".", ".", "D", "D", "D"]
stages = ["Rin", "Sch", "Trop"]

pca = PCA(n_components=3, copy=False, random_state=1, svd_solver="randomized")
X = pca.fit_transform(X)
fig, (ax1, ax2) = plt.subplots(
    ncols=2,
    figsize=(16, 6),
    subplot_kw={"aspect": 1})

for i in range(len(datasets)):
    ax1.plot(
        X[labels == i, 0], X[labels == i, 1],
        marker=".",
        markerfacecolor=colors[i],
        markeredgecolor="white", linewidth=0, label=stages[i],
        markeredgewidth=0,
        markersize=10,
        zorder=10,
        alpha=0.9)
    ax2.plot(
        X[labels == i, 2], X[labels == i, 1],
        marker=".",
        markerfacecolor=colors[i], markeredgecolor="white",
        linewidth=0, label=stages[i],
        markeredgewidth=0,
        markersize=10,
        zorder=10,
        alpha=0.9)

ax1.spines['right'].set_color('none')
ax1.spines['top'].set_color('none')
ax1.xaxis.set_ticks_position('bottom')
ax1.spines['bottom'].set_position(('data', 0))
ax1.yaxis.set_ticks_position('left')
ax1.spines['left'].set_position(('data', 0))

ax2.spines['right'].set_color('none')
ax2.spines['top'].set_color('none')
ax2.xaxis.set_ticks_position('bottom')
ax2.spines['bottom'].set_position(('data', 0))
ax2.yaxis.set_ticks_position('left')
ax2.spines['left'].set_position(('data', 0))

xlabels = ax1.get_xticklabels()
#xlabels = [label.get_text() if (i in [0, len(xlabels)-1]) or
#           int(label.get_text()) == 0
#           else "" for i, label in enumerate(xlabels)]

#ax1.set_xticklabels(xlabels)

for label in ax1.get_xticklabels() + ax1.get_yticklabels():
    label.set_bbox(dict(facecolor='white', edgecolor='None', alpha=0.65))
    label.set_zorder(0)

for label in ax2.get_xticklabels() + ax2.get_yticklabels():
    label.set_bbox(dict(facecolor='white', edgecolor='None', alpha=0.65))
    label.set_zorder(0)

lgnd = ax2.legend(prop={'size': 14}, bbox_to_anchor=(1.05, 1),
                  borderaxespad=0, numpoints=3)
for legend_handle in lgnd.legendHandles:
    legend_handle._legmarker.set_markersize(10)
    legend_handle._legmarker.set_alpha(1)

xmin, xmax = ax1.get_xlim()
ymin, ymax = ax1.get_ylim()

ax1.text(xmin - (xmax - xmin) * 0.05,
         0, "First component",
         fontweight="bold", rotation=90, horizontalalignment="center",
         verticalalignment="center")
ax1.text(0,
         ymin - (ymax - ymin) * 0.05,
         "Second component",
         fontweight="bold", horizontalalignment="center",
         verticalalignment="center")
ymin, ymax = ax2.get_ylim()
ax2.text(0,
         ymin - (ymax - ymin) * 0.05,
         "Third component",
         fontweight="bold", horizontalalignment="center",
         verticalalignment="center")

try:
    os.makedirs("images")
except OSError:
    pass
if args.chromosome is None:
    outfile = "images/pca_downsampled_distances_%s.png" % algorithm
else:
    outfile = "images/pca_downsampled_distances_%s_chr%02d.png" % (
        algorithm, args.chromosome)

fig.savefig(outfile)
