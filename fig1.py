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
SOURCE_DATA_DIR = os.path.join(ROOTDIR, 'source_data')
MAP_DIR = os.path.join(ROOTDIR, 'map')

# %% Functions

def load_data():
    file_path = os.path.join(SOURCE_DATA_DIR, 'def_blocs.csv')
    df = pd.read_csv(file_path, sep=";")
    df.rename(columns={'iso': 'iso_a3'}, inplace=True)
    return df

def load_map():
    file_path = os.path.join(MAP_DIR, 'naturalearth_lowres.shp')
    world = gpd.read_file(file_path)
    world = world[(world.pop_est>0) & (world.name!="Antarctica")]
    return world

def fig_blocs(df):
    fig, ax = plt.subplots()

    ax.set_aspect('equal')

    # axins1 = inset_axes(
    #     ax,
    #     width='3%',
    #     height='50%',
    #     loc='upper right',
    #     borderpad=0.0
    # )

    # axins2 = inset_axes(
    #     ax,
    #     width='3%',
    #     height='50%',
    #     loc='lower right',
    #     borderpad=0.0
    # )

    im1 = df.plot(
        column='Western',
        ax=ax,
        vmin=0,
        vmax=1.5,
        #legend=True,
        #cax=axins1,
        cmap='Blues',
        # missing_kwds={'color': 'lightgrey'}
        )

    im2 = df.plot(
        column='Eastern',
        ax=ax,
        vmin=0,
        vmax=2,
        #legend=True,
        #cax=axins2,
        cmap='Reds')

    im3 = df.plot(
        column='Neutral',
        ax=ax,
        vmin=0,
        vmax=1.2,
        #legend=True,
        #cax=axins2,
        cmap='spring')

    im3 = df.plot(
        column='Missing',
        ax=ax,
        vmin=0,
        vmax=5.0,
        #legend=True,
        #cax=axins2,
        cmap='Greys')

    # Set tick marks and labels
    # axins1.set_yticks([0.0, 0.2, 0.4, 0.6])
    # axins1.set_yticklabels(['$0$', '$0.2$', '$0.4$', '$0.6$'], fontsize=8)
    # axins2.set_yticks([-1.0, -2.0])
    # axins2.set_yticklabels(['$-1$', '$-2$'], fontsize=8)

    ax.set_axis_off()

    # Save
    fig.savefig(
        os.path.join(FIGDIR, 'fig1.png'),
        bbox_inches='tight',
        pad_inches=0
    )

    return fig



# %% Main

df = load_data()

# Merge data into world map
world = load_map()
df2 = world.merge(df, on='iso_a3', how='left')

df2['Western'] = 1 * (df2.bloc == 'Western')
df2.loc[df2['Western']==0, 'Western'] = np.nan
df2['Eastern'] = 1 * (df2.bloc == 'Eastern')
df2.loc[df2['Eastern']==0, 'Eastern'] = np.nan
df2['Neutral'] = 1 * (df2.bloc == 'Neutral')
df2.loc[df2['Neutral']==0, 'Neutral'] = np.nan
df2['Missing'] = 1 * (df2.bloc != 'Western') * (df2.bloc != 'Eastern') * (df2.bloc != 'Neutral')
df2.loc[df2['Missing']==0, 'Missing'] = np.nan
# Plot the figure
fig_blocs(df2)

# %%
