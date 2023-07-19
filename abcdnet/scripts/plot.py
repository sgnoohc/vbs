import argparse
import json
import sys
import os
srcdir = f"{os.getcwd()}/python"
if srcdir not in sys.path:
    sys.path.insert(0, srcdir)
    
import yahist
import uproot
import numpy as np
import awkward as ak
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({"figure.facecolor":  (1,1,1,0)})
from sklearn.metrics import roc_curve

from utils import VBSConfig
import ingress
import train

parser = argparse.ArgumentParser(description="Run inference")
parser.add_argument("config_json", type=str, help="config JSON")
parser.add_argument(
    "--epoch", type=int, default=50, metavar="N",
    help="training epoch of model to use for inference (default: 50)"
)
parser.add_argument(
    "--corr2D_rebin", type=int, default=0, help="rebin 2D correlation plots"
)
parser.add_argument(
    "--corr2D_ymax", type=float, default=8., 
    help="max value for y-axis of nominal 2D correlation plot"
)
parser.add_argument(
    "--corr2D_ymin", type=float, default=0., 
    help="min value for y-axis of nominal 2D correlation plot"
)
parser.add_argument(
    "--corr2D_ybins", type=int, default=20, 
    help="number of bins for the y-axis of the nominal 2D correlation plot"
)
parser.add_argument(
    "--loss_logy", action="store_true",
    help="make y-axis of loss curve log-scale"
)
args = parser.parse_args()

config = VBSConfig.from_json(args.config_json)
plots_dir = f"{config.basedir}/{config.name}/plots"

# Get history JSON
history_json = {}
with open(train.get_outfile(config, tag="history", ext="json")) as f:
    history_json = json.load(f)

# Get testing inferences
test_csv = train.get_outfile(
    config, 
    epoch=args.epoch, 
    tag="test_inferences",
    ext="csv",
    subdir="inferences"
)
test_df = pd.read_csv(test_csv)

# Get training inferences
train_csv = train.get_outfile(
    config, 
    epoch=args.epoch, 
    tag="train_inferences",
    ext="csv",
    subdir="inferences"
)
train_df = pd.read_csv(train_csv)

# Merge testing and training dfs
total_df = pd.concat([test_df, train_df])

# --- Plot loss curve ---
fig, axes = plt.subplots(figsize=(16, 12))

epochs = range(1, len(history_json["train_loss"])+1)

N = 10 # granularity of epochs shown (every Nth epoch is plotted)
axes.plot(epochs[::N], history_json["test_loss"][::N], label=r"$\mathcal{L}_{test} = \mathcal{L}_{BCE}$ + $\lambda$dCorr", color="C0");
axes.plot(epochs[::N], history_json["train_loss"][::N], label=r"$\mathcal{L}_{train} = \mathcal{L}_{BCE}$ + $\lambda$dCorr", color="C1");
axes.plot(epochs[::N], history_json["test_bce"][::N], label=r"$\mathcal{L}_{BCE}$", color="C0", linestyle="dashed")
axes.plot(epochs[::N], history_json["train_bce"][::N], label=r"$\mathcal{L}_{BCE}$", color="C1", linestyle="dashed")
axes.plot(epochs[::N], history_json["test_disco"][::N], label=r"$\lambda$dCorr$^2$", color="C0", linestyle="dotted")
axes.plot(epochs[::N], history_json["train_disco"][::N], label=r"$\lambda$dCorr$^2$", color="C1", linestyle="dotted")

axes.axvline(args.epoch, color="k", alpha=0.25)

axes.tick_params(axis="both", which="both", direction="in", labelsize=32, top=True, right=True)
axes.set_xlabel("Epoch", size=32);
axes.set_ylabel("Avg. Loss", size=32);
if args.loss_logy:
    axes.set_yscale("log")
axes.autoscale()
axes.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, fontsize=24)

plt.savefig(f"{plots_dir}/loss.png", bbox_inches="tight")
plt.savefig(f"{plots_dir}/loss.pdf", bbox_inches="tight")
print(f"Wrote loss curve to {plots_dir}/loss.png")
plt.close()


# --- Plot ROC curve ---
fig, axes = plt.subplots(figsize=(12,12))

# Plot training ROC curve
fpr, tpr, thresh = roc_curve(train_df.truth, train_df.score, sample_weight=train_df.weight)
axes.plot(fpr, tpr, label=f"DisCo train (AUC = {np.trapz(tpr, fpr):.2f})");

# Plot testing ROC curve
fpr, tpr, thresh = roc_curve(test_df.truth, test_df.score, sample_weight=test_df.weight)
axes.plot(fpr, tpr, label=f"DisCo test (AUC = {np.trapz(tpr, fpr):.2f})");

total_sig = total_df[total_df.truth == 1].weight.sum()
total_bkg = total_df[total_df.truth == 0].weight.sum()
axes.scatter(fpr[tpr >= 5/total_sig][0], tpr[tpr >= 5/total_sig][0], marker="*", s=100, zorder=104, color="c", label=f"5 sig ({fpr[tpr >= 5/total_sig][0]*total_bkg:0.2f} bkg)");
axes.scatter(fpr[tpr >= 4/total_sig][0], tpr[tpr >= 4/total_sig][0], marker="v", s=100, zorder=103, color="b", label=f"4 sig ({fpr[tpr >= 4/total_sig][0]*total_bkg:0.2f} bkg)");
axes.scatter(fpr[tpr >= 3/total_sig][0], tpr[tpr >= 3/total_sig][0], marker="^", s=100, zorder=102, color="m", label=f"3 sig ({fpr[tpr >= 3/total_sig][0]*total_bkg:0.2f} bkg)");
axes.scatter(fpr[tpr >= 2/total_sig][0], tpr[tpr >= 2/total_sig][0], marker="o", s=100, zorder=101, color="r", label=f"2 sig ({fpr[tpr >= 2/total_sig][0]*total_bkg:0.2f} bkg)");
axes.scatter(fpr[tpr >= 1/total_sig][0], tpr[tpr >= 1/total_sig][0], marker="s", s=100, zorder=100, color="y", label=f"1 sig ({fpr[tpr >= 1/total_sig][0]*total_bkg:0.2f} bkg)");

# Format axes
axes.tick_params(axis="both", which="both", direction="in", labelsize=32, top=True, right=True)
axes.set_xlabel("Background efficiency", size=32);
axes.set_ylabel("Signal efficiency", size=32);
axes.legend(fontsize=24);

plt.savefig(f"{plots_dir}/roc_epoch{args.epoch}.png", bbox_inches="tight")
plt.savefig(f"{plots_dir}/roc_epoch{args.epoch}.pdf", bbox_inches="tight")
print(f"Wrote ROC curve to {plots_dir}/roc_epoch{args.epoch}.png")
plt.close()


# --- Plot score histogram ---
fig, axes = plt.subplots(figsize=(12, 12))

bins = np.linspace(0, 1, 101)

test_n_true = test_df[test_df.truth == 1].weight.sum()
test_n_false = test_df[test_df.truth == 0].weight.sum()
axes.hist(
    test_df[test_df.truth == 1].score, 
    weights=test_df[test_df.truth == 1].weight/test_n_true, 
    bins=bins,
    histtype="step",
    linewidth=2,
    label="test (sig)"
);
axes.hist(
    test_df[test_df.truth == 0].score, 
    weights=test_df[test_df.truth == 0].weight/test_n_false, 
    bins=bins,
    histtype="step",
    linewidth=2,
    label="test (bkg)"
);
axes.set_yscale("log");
axes.legend(fontsize=16)

axes.tick_params(axis="both", which="both", direction="in", labelsize=20, top=True, right=True)
axes.set_xlabel("score", size=20);
axes.set_ylabel("a.u.", size=20);

plt.savefig(f"{plots_dir}/scores_epoch{args.epoch}.png", bbox_inches="tight")
plt.savefig(f"{plots_dir}/scores_epoch{args.epoch}.pdf", bbox_inches="tight")
print(f"Wrote scores histogram to {plots_dir}/scores_epoch{args.epoch}.png")
plt.close()

# --- Plot correlation histogram ---
fig, axes = plt.subplots(figsize=(12, 12))

xbins = [0, 0.05, 0.1, 0.2, 0.4, 0.7, 1]
ybins = np.linspace(args.corr2D_ymin, args.corr2D_ymax, args.corr2D_ybins+1)

hist = yahist.Hist2D(
    (test_df[test_df.truth == 0].score, test_df[test_df.truth == 0].disco_target), 
    weights=test_df[test_df.truth == 0].weight/test_n_false, 
    bins=[xbins, ybins]
);

if args.corr2D_rebin > 0:
    hist = hist.rebin(args.corr2D_rebin)

hist.plot(ax=axes, hide_empty=False)

hist.profile().plot(errors=True, color="w", fmt=".", label=r"DisCo target profile")

axes.set_xlabel("ABCDNet score", size=32)
axes.set_ylabel("DisCo target", size=32)
axes.set_xlim([0, 1])
axes.set_ylim([args.corr2D_ymin, args.corr2D_ymax])

plt.savefig(f"{plots_dir}/correlation2D_epoch{args.epoch}.png", bbox_inches="tight")
plt.savefig(f"{plots_dir}/correlation2D_epoch{args.epoch}.pdf", bbox_inches="tight")
print(f"Wrote correlation histogram to {plots_dir}/correlation2D_epoch{args.epoch}.png")
plt.close()


# --- Plot other correlation histogram ---
fig, axes = plt.subplots(figsize=(12, 12))

hist.transpose().plot(ax=axes, hide_empty=False)

hist.transpose().profile().plot(errors=True, color="w", fmt=".", label="ABCDNet score profile")

axes.set_ylabel("ABCDNet score", size=32)
axes.set_xlabel("DisCo target", size=32)
axes.set_ylim([0, 1]);
axes.set_xlim([args.corr2D_ymin, args.corr2D_ymax])

plt.savefig(f"{plots_dir}/correlation2D_flipped_epoch{args.epoch}.png", bbox_inches="tight")
plt.savefig(f"{plots_dir}/correlation2D_flipped_epoch{args.epoch}.pdf", bbox_inches="tight")
print(f"Wrote other correlation histogram to {plots_dir}/correlation2D_flipped_epoch{args.epoch}.png")
plt.close()

# --- Wrap up ---
print(f"Done. All plots can be found here: {plots_dir}")
