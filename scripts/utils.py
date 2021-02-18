import itertools as it
import math
import numpy as np
import os
import pandas as pd

from bokeh.models import ColorBar, ColumnDataSource, HoverTool, LabelSet, Legend
from bokeh.models import LinearColorMapper, FactorRange
from bokeh.palettes import Magma, d3
from bokeh.plotting import figure, output_file, show

def bokeh_correlation_heatmap(table, group, label, norm_method, outdir):
    for corr in ['spearman']:
        figure_table = table.T.corr(norm_method=corr).stack()
        figure_table.index.names = ['0', '1']
        figure_table = figure_table.reset_index()
        figure_table.columns = ['xname', 'yname', 'target']
        figure_table['target'] = np.round(figure_table['target'], 2)

        figure_table.loc[:, 'group_x'] = group[figure_table.xname].values
        figure_table.loc[:, 'group_y'] = group[figure_table.yname].values
        figure_table.loc[:, 'color'] = '#000000'
        figure_table.loc[figure_table['target'] >= 0.5, 'color'] = '#FFFFFF'

        x_data = [(data.group_x, data.xname) for i, data in figure_table.iterrows()]
        x_uniq = [(data.group_y, data.yname) for i, data in figure_table.loc[figure_table.xname == figure_table.xname[0]].iterrows()]
        
        y_data = [(data.group_y, data.yname) for i, data in figure_table.iterrows()]

        mapper = LinearColorMapper(palette=list(Magma[256])[256:128:-1], low=0, high=1)

        title = ''.join((
            f'{label.capitalize()} Lowly Expressed Genes {corr.capitalize()} ',
            f'Correlations ({norm_method.capitalize()})'
        ))

        data_source = ColumnDataSource(
            data=dict(
                xname=x_data,
                yname=y_data,
                target=figure_table.target
            )
        )

        fig = figure(
            #title=title,
            title='',
            x_axis_location="below",
            tools="hover, save",
            x_range=FactorRange(*x_uniq),
            y_range=FactorRange(*x_uniq),
            tooltips=[
                ('gene A', '@xname'),
                ('gene B', '@yname'),
                (corr, '@target')
            ])

        fig.plot_width = len(figure_table.xname.unique()) * 9
        fig.plot_height = len(figure_table.yname.unique()) * 9
        fig.grid.grid_line_color = None
        fig.axis.axis_line_color = None
        fig.axis.major_tick_line_color = None
        fig.axis.major_label_text_font_size = "6px"
        fig.axis.major_label_standoff = 0
        fig.axis.major_label_text_alpha = 0
        fig.xaxis.major_label_orientation = np.pi / 2
        #fig.xaxis.major_label_text_alpha = 0
        #fig.xaxis.major_label_text_font_size = '6px'
        #fig.axis.group_text_font_size = '32px'
        fig.axis.group_text_font_size = f'{int(len(figure_table.xname.unique())/4.5)}px'
        fig.axis.group_text_align = 'center'

        fig.axis.major_label_standoff = 0
        #fig.title.text_font_size = '36px'
        fig.title.text_font_size = f'{int(len(figure_table.xname.unique())/7)}px'

        fig.rect(
            x_data='xname',
            y_data='yname',
            width=0.9,
            height=0.9,
            source=data_source,
            line_color=None,
            fill_color={'field': 'target', 'transform': mapper},
            hover_line_color='black'
        )

        labels = LabelSet(
            x='xname',
            y='yname',
            text='target',
            text_align='center',
            text_font_size='8px',
            text_color='color',
            x_offset=0,
            y_offset=-5,
            source=data_source,
            render_mode='canvas'
        )

        color_bar = ColorBar(
            color_mapper=mapper,
            #major_label_text_font_size="24px",
            #label_standoff=16,
            width=int(len(figure_table.xname.unique())/5),
            major_label_text_font_size=f'{int(len(figure_table.xname.unique())/6)}px',
            label_standoff=int(len(figure_table.xname.unique())/15),
            border_line_color=None,
            location=(0, 0)
        )

        fig.add_layout(labels)
        fig.add_layout(color_bar, 'right')

        file_name = os.path.join(outdir, f"{label}_{corr}_heatmap_{norm_method.lower().replace(' ', '_')}.html")
        output_file(file_name, title=title)

        show(fig)

        return

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

def bokeh_spikein_lineplot(outdir, spikeins, title):
    num_genes = spikeins.shape[0]
    palettes = d3['Category20'][20]

    palettes = list(it.islice(it.cycle(palettes), num_genes))
    shapes = list(it.islice(it.cycle(['solid', 'dashed', 'dotted']), num_genes))

    data_source = ColumnDataSource(np.log2(spikeins.T + 1))
    data_source.data['file_accession'] = list(data_source.data['file_accession'])

    fig = figure(
        title=title,
        x_axis_location='below',
        tools='hover,save',
        x_range=list(spikeins.columns),
        tooltips=[
            ('sample', '$x'),
            ('expression', '$y'),
            ('spikein', '$name')
        ]
    )

    legends = []
    for gene, color, shape in zip(list(spikeins.index), palettes, shapes):
        line = fig.line(
            x='file_accession',
            y=gene,
            name=gene,
            line_width=2,
            line_dash=shape,
            color=color,
            source=data_source
        )
        legends.append((gene, [line]))

    legends = Legend(items=legends, label_text_font_size='8pt')
    legends.click_policy="hide"

    fig.add_layout(legends, 'left')

    fig.plot_width = 1920
    fig.plot_height = 1080
    fig.grid.grid_line_color = None
    fig.axis.axis_line_color = None
    fig.axis.major_tick_line_color = None
    fig.axis.major_label_text_font_size = "10px"
    fig.axis.major_label_standoff = 0
    fig.xaxis.major_label_orientation = np.pi / 3
    fig.yaxis.axis_label = 'log\u2082(Count + 1)'

    output_file(
        os.path.join(outdir, f'{title.lower().replace(" ", "_")}.html'), 
        title=title
    )

    show(fig)

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