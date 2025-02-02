#!/usr/bin/env/python3

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def plot_distribution(df, columnname, output, label="Distribution", nbins=20):
    sns.set_theme()
    ax = sns.displot(df, x=columnname, bins=nbins)
    ax.savefig(output)
    plt.close()
