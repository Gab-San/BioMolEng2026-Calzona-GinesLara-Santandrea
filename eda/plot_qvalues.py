from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
eda = Path(__file__).resolve().parent
source = eda.parent / "pipeline/output/postprocessed"
variants = sorted((path for path in source.iterdir() if path.is_dir()),
                  key=lambda path: (path.name != "WT", path.name.count("_"), path.name))
labels, samples = [path.name for path in variants], []
for variant in variants:
    time, q = np.loadtxt(variant / "analysis/qvalue.dat", skiprows=1, usecols=(0, 1), unpack=True)
    samples.append(q[time >= time.max() - 2500])
fig, ax = plt.subplots(figsize=(max(8, 0.53 * len(labels)), 5.5), layout="constrained")
ax.boxplot(samples, tick_labels=labels, showmeans=True, meanline=True,
           meanprops={"linestyle": "-", "color": "black"}, medianprops={"visible": False})
ax.axhline(samples[labels.index("WT")].mean(), color="tab:blue", linestyle="--", label="WT mean")
ax.set(xlabel="Mutation", ylabel="Q value", title="Q values over the final 2500 ps")
ax.tick_params(axis="x", labelrotation=90)
ax.legend()
(eda / "output").mkdir(exist_ok=True)
fig.savefig(eda / "output/qvalue_boxplots.png", dpi=180)
