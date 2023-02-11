import numpy as np
import matplotlib.pyplot as plt

data = {
    "Barton LLC": 109438.50,
    "Frami, Hills and Schmidt": 103569.59,
    "Fritsch, Russel and Anderson": 112214.71,
    "Jerde-Hilpert": 112591.43,
    "Keeling LLC": 100934.30,
    "Koepp Ltd": 103660.54,
    "Kulas Inc": 137351.96,
    "Trantow-Barrows": 123381.38,
    "White-Trantow": 135841.99,
    "Will LLC": 104437.60
}

group_data = list(data.values())
group_names = list(data.keys())
group_mean = np.mean(group_data)

fig, ax = plt.subplots()
ax.barh(group_names, group_data)
plt.savefig("IMG_01.png", dpi=150)

# Some customization
plt.style.use("ggplot")
fig, ax = plt.subplots()
ax.barh(group_names, group_data)
labels = ax.get_xticklabels()
plt.setp(labels, rotation=45, horizontalalignment="right")
plt.savefig("IMG_02.png", dpi=150)

# Update style
plt.rcParams.update({"figure.autolayout": True})
fig, ax = plt.subplots()
ax.barh(group_names, group_data)
labels = ax.get_xticklabels()
plt.setp(labels, rotation=45, horizontalalignment="right")
plt.savefig("IMG_03.png", dpi=150)

# Add labels
fig, ax = plt.subplots()
ax.barh(group_names, group_data)
labels = ax.get_xticklabels()
plt.setp(labels, rotation=45, horizontalalignment="right")
ax.set(xlim=[-10000, 140000],
       xlabel="Total revenue",
       ylabel="Company",
       title="Company Revenue")
plt.savefig("IMG_04.png", dpi=150)


# Change figure size
fig, ax = plt.subplots(figsize=(8,4))
ax.barh(group_names, group_data)
labels = ax.get_xticklabels()
plt.setp(labels, rotation=45, horizontalalignment="right")
ax.set(xlim=[-10000, 140000],
       xlabel="Total revenue",
       ylabel="Company",
       title="Company Revenue")
plt.savefig("IMG_05.png", dpi=150)

# For label formatting
def currency(x, pos):
    """The two arguments are the value and tick position"""
    if x >= 1e6:
        s = "${:1.1f}M".format(x*1e-6)
    else:
        s = "${:1.0f}K".format(x*1e-3)
    return s


fig, ax = plt.subplots(figsize=(6, 8))
ax.barh(group_names, group_data)
labels = ax.get_xticklabels()
plt.setp(labels, rotation=45, horizontalalignment="right")
ax.set(xlim=[-10000, 140000],
       xlabel="Total Revenue",
       ylabel="Company",
       title="Company Revenue")
ax.xaxis.set_major_formatter(currency)
plt.savefig("IMG_06.png", dpi=150)

