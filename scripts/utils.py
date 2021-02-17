import itertools as it
import math
import numpy as np
import os
import pandas as pd
from bokeh.palettes import Magma, d3
from bokeh.plotting import figure, output_file, show

def bokeh_xprs_distribution(dataset_name, outdir, xprs, groups):
    """
    Create Figure for the Distribution of Expression
    """
    fig = figure(
        plot_width=1280, plot_height=720,
        title=f'Distribution of Read Count ({dataset_name})',
        x_axis_location="below",
        tools="save",
    )

    palettes = d3['Category20'][20]
    palettes = list(it.islice(it.cycle(palettes), len(groups.unique())))
    for tissue, color in zip(groups.unique(), palettes):
        group = groups.loc[groups == tissue]
        xprs_group = xprs[group.index]
        density, edges = np.histogram(
            np.log2(xprs_group.mean(axis=1) + 1), 
            density=True, bins=50
        )
        edges_average = (edges[:-1] + edges[1:]) / 2

        _ = fig.line(
            x=edges_average,
            y=density,
            color=color,
            line_width=3,
            alpha=0.8,
            legend_label=tissue
        )

    fig.legend.location = "top_right"
    fig.legend.click_policy = "hide"
    fig.xaxis.axis_label = 'log₂(count + 1)'
    fig.yaxis.axis_label = 'Density'
    fig.axis.axis_label_text_font_size = "20px"
    fig.axis.major_label_text_font_size = "16px"
    fig.title.text_font_size = '20px'

    output_file(
        os.path.join(outdir, f'{dataset_name.lower()}_xprs_distribution.html'),
        title=f'{dataset_name} Expression Distribution'
    )
    show(fig)
    
    return

#%%
def bokeh_num_nonexpressed_genes(dataset_name, outdir, xprs, groups):
    """
    Create Figure for the Boxplot of the Number of Non-expressed Genes
    """
    def get_outliers(group):
        cat = group.name
        condition_a = (group['count'] > upper.loc[cat]['count'])
        condition_b = (group['count'] < lower.loc[cat]['count'])
        selector = condition_a | condition_b
        return group[selector]['count']

    palettes = d3['Category20'][20]
    palettes = list(it.islice(it.cycle(palettes), len(groups.unique())))
    tissue_zero_genes = []
    tissue_colors = []
    for tissue, color in zip(groups.unique(), palettes):
        group = groups.loc[groups == tissue]
        xprs_group = xprs[group.index]
        tissue_zero_genes.append(pd.DataFrame({
            'count': (xprs_group == 0).sum(),
            'tissue': tissue,
        }))
        tissue_colors.append([color])

    tissue_zero_genes_df = pd.concat(tissue_zero_genes)
    tissue_colors_df = pd.DataFrame(tissue_colors, index=groups.unique())
    

    fig = figure(
        title=f'{dataset_name} Number of Non-expressed Genes by Tissue', 
        background_fill_color="#fafafa",

        x_range=tissue_zero_genes_df['tissue'].unique(),
        plot_width=1280, plot_height=720
    )

    tissues = tissue_zero_genes_df.groupby('tissue')
    q1 = tissues.quantile(q=0.25)
    q2 = tissues.quantile(q=0.5)
    q3 = tissues.quantile(q=0.75)
    iqr = q3 - q1
    upper = q3 + 1.5 * iqr
    lower = q1 - 1.5 * iqr

    outliers = tissues.apply(get_outliers).dropna()

    if not outliers.empty:
        outx = []
        outy = []
        for keys in outliers.index:
            outx.append(keys[0])
            outy.append(outliers.loc[keys[0]].loc[keys[1]])
        fig.circle(x=outx, y=outy, size=6, color="gray", fill_alpha=0.6)

    qmin = tissues.quantile(q=0.00)
    qmax = tissues.quantile(q=1.00)
    upper['count'] = [min([x,y]) for (x, y) in zip(list(qmax.loc[:, 'count']), upper['count'])]
    lower['count'] = [max([x,y]) for (x, y) in zip(list(qmin.loc[:, 'count']), lower['count'])]

    fig.segment(x0=upper.index, y0=upper['count'], x1=upper.index, y1=q3['count'], line_color="black")
    fig.segment(x0=upper.index, y0=lower['count'], x1=upper.index, y1=q1['count'], line_color="black")

    fill_color = tissue_colors_df.loc[q2.index].values.reshape(-1)

    fig.rect(x=q2.index, y=(q3['count'] + q2['count']) / 2, width=0.7, height=(q3['count'] - q2['count']), fill_color=fill_color, line_color="black")
    fig.rect(x=q2.index, y=(q1['count'] + q2['count']) / 2, width=0.7, height=(q2['count'] - q1['count']), fill_color=fill_color, line_color="black")
    
    fig.rect(x=lower.index, y=lower['count'], width=0.2, height=0.01, line_color="black")
    fig.rect(x=upper.index, y=upper['count'], width=0.2, height=0.01, line_color="black")

    fig.xgrid.grid_line_color = None
    fig.ygrid.grid_line_color = "white"
    fig.grid.grid_line_width = 2
    fig.xaxis.major_label_orientation = math.pi / 4
    fig.xaxis.major_label_text_font_size="16px"

    output_file(
        os.path.join(outdir, f'{dataset_name.lower()}_number_of_non_expressed_genes.html'),
        title=f'{dataset_name} Boxplot of Number of Non-Expressed Genes'
    )
    show(fig)

    return

def extract_tissue_exclusive_gene(
    xprs,
    tissues,
    high_threshold=10,
    low_threshold=1
) -> pd.Index:
    index = {}
    for tissue_a in tissues.unique():
        xprs_tissue_a = xprs.loc[:, (tissues == tissue_a)]
        gene_filter_a = xprs_tissue_a.quantile(0.5, axis=1) >= high_threshold
        for tissue_b in tissues.unique():
            if tissue_a == tissue_b:
                continue

            xprs_tissue_b = xprs.loc[:, (tissues == tissue_b)]
            gene_filter_b = xprs_tissue_b.quantile(0.5, axis=1) <= low_threshold
            gene_filter_a = gene_filter_a & gene_filter_b

        if (gene_filter_a.sum() >= 1):
            filtered = xprs.loc[gene_filter_a]
            index[tissue_a] = list(filtered.index)

    return index