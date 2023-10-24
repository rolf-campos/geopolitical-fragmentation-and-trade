# %% Imports and constants

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['figure.dpi'] = 72*6

ROOTDIR = os.path.realpath('.')
FIGDIR = os.path.join(ROOTDIR, 'results')
SOURCE_DATA_DIR = os.path.join(ROOTDIR, 'results')


# %% Functions

def load_data():
    file_path = os.path.join(SOURCE_DATA_DIR, 'fragmentation_3_prod.csv')
    df = pd.read_csv(file_path)
    return df

def figA_trade(df):
    fig, ax = plt.subplots()

    s_1 = df['Trade']
    s_2 = df['Trade_prod']

    ax.scatter(s_1, s_2, color="blue", alpha=0.2)

    for i, txt in enumerate(df['iso']):
        if s_1[i] < -8 or s_2[i] < -8:
            ax.annotate(txt, (s_1[i], s_2[i]))

    v_min = -45
    v_max = 5
    ax.set_xlim(v_min, v_max)
    ax.set_ylim(v_min, v_max)
    ax.plot(range(v_min, v_max+1), range(v_min, v_max+1), color="black", alpha=0.5, lw=1)
    ax.axhline(0, v_min, v_max, color="k", lw=0.5, ls="--")
    ax.axvline(0, v_min, v_max, color="k", lw=0.5, ls="--")
    ax.set_xlabel("Baseline results")
    ax.set_ylabel("Replication with data from TiVA")
    rotn = np.degrees(np.arctan2(1, 1))
    ax.text(v_min+2, v_min+3, f"$45^o$", fontsize=8, va='bottom',
               rotation=rotn, rotation_mode='anchor', transform_rotates_text=True)


    # Save
    fig.savefig(
        os.path.join(FIGDIR, 'figA_trade.png'),
        bbox_inches='tight',
        pad_inches=0
    )

    return fig

def figA_welfare(df):
    fig, ax = plt.subplots()

    s_1 = df['Welfare']
    s_2 = df['Welfare_prod']

    ax.scatter(s_1, s_2, color="blue", alpha=0.2)

    for i, txt in enumerate(df['iso']):
        if s_1[i] < -2.5 or s_2[i] < -2.5:
            ax.annotate(txt, (s_1[i], s_2[i]))

    v_min = -8
    v_max = 1
    ax.set_xlim(v_min, v_max)
    ax.set_ylim(v_min, v_max)
    ax.plot(range(v_min, v_max+1), range(v_min, v_max+1), color="black", alpha=0.5, lw=1)
    ax.axhline(0, v_min, v_max, color="k", lw=0.5, ls="--")
    ax.axvline(0, v_min, v_max, color="k", lw=0.5, ls="--")
    ax.set_xlabel("Baseline results")
    ax.set_ylabel("Replication with data from TiVA")
    rotn = np.degrees(np.arctan2(1, 1))
    ax.text(v_min+0.3, v_min+0.4, f"$45^o$", fontsize=8, va='bottom',
               rotation=rotn, rotation_mode='anchor', transform_rotates_text=True)

    # Save
    fig.savefig(
        os.path.join(FIGDIR, 'figA_welfare.png'),
        bbox_inches='tight',
        pad_inches=0
    )

    return fig

# %% Main

df = load_data()

# Plot the trade figure
figA_trade(df)
figA_welfare(df)


# %%
