# %% Imports and constants

import os
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams['figure.dpi'] = 72*6

ROOTDIR = os.path.realpath('.')
FIGDIR = os.path.join(ROOTDIR, 'results')
SOURCE_DATA_DIR = os.path.join(ROOTDIR, 'results')
MAP_DIR = os.path.join(ROOTDIR, 'map')

# %% Functions

def load_data():
    file_path = os.path.join(SOURCE_DATA_DIR, 'fragmentation_3.csv')
    df = pd.read_csv(file_path)
    df.rename(columns={'iso': 'iso_a3'}, inplace=True)
    return df

def load_map():
    file_path = os.path.join(MAP_DIR, 'naturalearth_lowres.shp')
    world = gpd.read_file(file_path)
    world = world[(world.pop_est>0) & (world.name!="Antarctica")]
    return world

def fig1(df):
    fig, ax = plt.subplots()

    ax.set_aspect('equal')

    axins1 = inset_axes(
        ax,
        width='3%',
        height='50%',
        loc='upper right',
        borderpad=0.0
    )

    axins2 = inset_axes(
        ax,
        width='3%',
        height='50%',
        loc='lower right',
        borderpad=0.0
    )

    im1 = df.plot(
        column='W+',
        vmax=0.6,
        ax=ax,
        legend=True,
        cax=axins1,
        cmap='Blues',
        missing_kwds={'color': 'lightgrey'})

    im2 = df.plot(
        column='W-',
        vmin=-15,
        ax=ax,
        legend=True,
        cax=axins2,
        cmap='Reds_r')

    # Set tick marks and labels
    axins1.set_yticks([0.0, 0.2, 0.4, 0.6])
    axins1.set_yticklabels(['$0$', '$0.2$', '$0.4$', '$0.6$'], fontsize=8)
    axins2.set_yticks([-5, -10, -15])
    axins2.set_yticklabels(['$-5$', '$-10$', '$-15$'], fontsize=8)

    ax.set_axis_off()

    # Save
    fig.savefig(
        os.path.join(FIGDIR, 'fig2.png'),
        bbox_inches='tight',
        pad_inches=0
    )

    return fig



# %% Main

df = load_data()

# Separate into positive and negative values
df['W+'] = df['Welfare']
df.loc[df['W+']<=0, 'W+'] = np.nan
df['W-'] = df['Welfare']
df.loc[df['W-']>=0, 'W-'] = np.nan

# Merge data into world map
world = load_map()
df2 = world.merge(df, on='iso_a3')

# Plot the figure
fig1(df2)

# %%
