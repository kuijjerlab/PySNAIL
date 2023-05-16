import itertools as it
import math
import numpy as np
import os
import pandas as pd
import scipy.stats as stats

from bokeh.models import ColorBar, ColumnDataSource, HoverTool, LabelSet, Legend
from bokeh.models import LinearColorMapper, FactorRange, Span
from bokeh.palettes import Magma, d3
from bokeh.plotting import figure, output_file, save
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from tqdm import tqdm

def pearson_pvalue(a, b, only_positive=True):
    corr, pvalue = stats.pearsonr(a, b)
    if only_positive:
        if corr <= 0:
            pvalue = 1
    return pvalue

def spearman_pvalue(a, b, only_positive=True):
    corr, pvalue = stats.spearmanr(a, b)
    if only_positive:
        if corr <= 0:
            pvalue = 1
    return pvalue

def bokeh_curve(
    outdir,
    file_prefix,
    corrs,
    ground_truth,
    significant_boundary=-1.0,
    legend_labels=['Qsmooth', 'SNAIL' 'RLE', 'TMM'],
    colors=['blueviolet', 'crimson', 'chocolate', 'navy'],
    method='ROC',
):

    if method.lower() == 'roc':
        auc_fn = roc_curve
        metric_x_label = 'False Positive Rate'
        metric_y_label = 'True Positive Rate'
    else:
        auc_fn = precision_recall_curve
        metric_x_label = 'Recall'
        metric_y_label = 'Precision'
    
    title = ''

    fig = figure(
        plot_width=800, plot_height=800,
        title=title,
        x_axis_location="below",
        tools=["pan", "wheel_zoom", "box_zoom", "reset", "hover", "save"],
        tooltips=[
            (metric_x_label, '$x'),
            (metric_y_label, '$y'),
            ('Normalization', '@label')
        ]
    )
    
    corr_truth_binarized = ground_truth['Ground Truth'].astype(int)

    auc_results = {}
    significant_boundary_points = []
    vertical_line = -1
    n_corrs = len(corrs)
    for corr, label, color in zip(corrs[0:n_corrs], legend_labels[0:n_corrs], colors[0:n_corrs]):
        auc_list = list()
        num_edges = list()
        cross_index = -1
        metric_a, metric_b, threshold = auc_fn(
            corr_truth_binarized.values, 
            np.abs(corr.values)
        )
        cross_index = np.min(np.where(threshold >= significant_boundary)[0])
        if method.lower() == 'roc':
            metric_a_label = metric_x_label
            metric_b_label = metric_y_label
        else:
            metric_a_label = metric_y_label
            metric_b_label = metric_x_label
        
        table = pd.DataFrame({
            metric_a_label: metric_a,
            metric_b_label: metric_b,
            'method': label
        })
        if method.lower() == 'roc':
            x = table[metric_a_label].values
            y = table[metric_b_label].values
        else:
            x = table[metric_b_label].values
            y = table[metric_a_label].values
        
        score = auc(x, y)
        #auc_list.append([
        #        threshold,
        #        score,
        #        label])
        #if threshold >= significant_boundary and not cross_significance:
        #    cross_score = score
        #    cross_significance = True
        #    vertical_line = threshold
        #    significant_boundary_points.append([label, color, threshold, score])
        #auc_df = pd.DataFrame(auc_list, columns=['threshold', 'auc', 'label'])
        #legend_label = label
        #if cross_score >= 0.9995:
        #    cross_score = 0.999
        legend_label = f'{label} (AUC: {score:>.3f})'
        
        data_source = ColumnDataSource(
            data=dict(
                x=x,
                y=y,
                label=table['method']
            )
        )

        fig.line(
            x='x',
            y='y',
            color=color, alpha=0.5, line_width=5,
            legend_label=legend_label,
            source=data_source,
        )
        if cross_index > 0:
            #vline = Span(location=vertical_line, dimension='height', line_color='grey', line_dash='dashed', line_width=2)
            #fig.renderers.extend([vline])
            fig.asterisk(x[cross_index], y[cross_index], size=20, color=color, alpha=1, line_width=3)
    #for point in significant_boundary_points:
    #vline = Span(location=vertical_line, dimension='height', line_color='grey', line_dash='dashed', line_width=2)
    #fig.renderers.extend([vline])
    #for point in significant_boundary_points:
    #    fig.asterisk([point[2]], [point[3]], size=20, color=point[1], alpha=1, line_width=3)
    
    fig.legend.location = 'bottom_left'
    fig.xaxis.axis_label = metric_x_label
    fig.yaxis.axis_label = metric_y_label
    fig.xaxis.axis_label_text_font_size = '20pt'
    fig.yaxis.axis_label_text_font_size = '20pt'
    fig.legend.label_text_font_size = '18pt'
    fig.axis.major_label_text_font_size = '12pt'
    fig.title.text_font_size = '20pt'
    fig.legend.click_policy="hide"

    outfile = os.path.join(
        outdir,
        f'{file_prefix.lower()}_{method.lower()}_curve.html'
    )
    output_file(outfile, title=f'{method.upper()} Curve')
    save(fig)

    return

def bokeh_area_under_curve(
    outdir,
    corrs,
    significant_boundary=0.5,
    legend_labels=['Validation', 'Qsmooth', 'SNAIL' 'RLE', 'TMM'],
    colors=['green', 'blueviolet', 'crimson', 'chocolate', 'navy'],
    method='ROC',
):

    if method.lower() == 'roc':
        auc_fn = roc_curve
        metric_x_label = 'FPR'
        metric_y_label = 'TPR'
    else:
        auc_fn = precision_recall_curve
        metric_x_label = 'Recall'
        metric_y_label = 'Precision'
    
    title = ''

    fig = figure(
        plot_width=800, plot_height=800,
        title=title,
        x_axis_location="below",
        tools=["pan", "wheel_zoom", "box_zoom", "reset", "hover", "save"],
        tooltips=[
            (metric_x_label, '$x'),
            (metric_y_label, '$y'),
            ('Normalization', '@label')
        ]
    )
    
    corr_truth = corrs[0]
    auc_results = {}
    significant_boundary_points = []
    vertical_line = -1
    n_corrs = len(corrs)
    for corr, label, color in zip(corrs[1:n_corrs], legend_labels[1:n_corrs], colors[1:n_corrs]):
        print(label)
        auc_list = list()
        num_edges = list()
        cross_significance = False
        cross_score = -1
        for threshold in np.linspace(0.2, 0.8, 101):
            corr_truth_binarized = corr_truth.copy()
            corr_truth_binarized[np.abs(corr_truth_binarized) >= threshold] = 1
            corr_truth_binarized[np.abs(corr_truth_binarized) < threshold] = 0

            metric_a, metric_b, _ = auc_fn(
                corr_truth_binarized.values, 
                np.abs(corr.values)
            )
            num_edges.append([
                threshold,
                (np.abs(corr_truth_binarized) >= threshold).sum(),
                (np.abs(corr_truth_binarized) < threshold).sum()
            ])
            if method.lower() == 'roc':
                metric_a_label = metric_x_label
                metric_b_label = metric_y_label
            else:
                metric_a_label = metric_y_label
                metric_b_label = metric_x_label
            
            table = pd.DataFrame({
                metric_a_label: metric_a,
                metric_b_label: metric_b,
                'method': label
            })
            if method.lower() == 'roc':
                score = auc(table[metric_a_label].values, table[metric_b_label].values)
            else:
                score = auc(table[metric_b_label].values, table[metric_a_label].values)
            auc_list.append([
                    threshold,
                    score,
                    label])
            if threshold >= significant_boundary and not cross_significance:
                cross_score = score
                cross_significance = True
                vertical_line = threshold
                significant_boundary_points.append([label, color, threshold, score])
        auc_df = pd.DataFrame(auc_list, columns=['threshold', 'auc', 'label'])
        legend_label = label
        if cross_score >= 0.9995:
            cross_score = 0.999
        legend_label = f'{label} (AUC*: {cross_score:>.3f})'
        
        data_source = ColumnDataSource(
            data=dict(
                threshold=auc_df.threshold.values,
                auc=auc_df.auc.values,
                label=auc_df.label
            )
        )

        fig.line(
            x='threshold',
            y='auc',
            color=color, alpha=0.5, line_width=5,
            legend_label=legend_label,
            source=data_source,
        )
        auc_results[label] = auc_df

    vline = Span(location=vertical_line, dimension='height', line_color='grey', line_dash='dashed', line_width=2)
    fig.renderers.extend([vline])
    for point in significant_boundary_points:
        fig.asterisk([point[2]], [point[3]], size=20, color=point[1], alpha=1, line_width=3)
    
    fig.legend.location = 'bottom_left'
    fig.xaxis.axis_label = 'Threshold applied to the validation dataset'
    fig.yaxis.axis_label = f'AU{method.upper()}'
    fig.xaxis.axis_label_text_font_size = '20pt'
    fig.yaxis.axis_label_text_font_size = '20pt'
    fig.legend.label_text_font_size = '18pt'
    fig.axis.major_label_text_font_size = '12pt'
    fig.title.text_font_size = '20pt'
    fig.legend.click_policy="hide"

    outfile = os.path.join(
        outdir,
        f'au{method.lower()}.html'
    )
    output_file(outfile, title=f'AU{method.upper()}')
    save(fig)

    columns = [
        'Threshold',
        'Number of Postive Associations',
        'Number of Negative Associations'
    ]
    return pd.DataFrame(num_edges, columns=columns)

def bokeh_correlation_heatmap(xprs, gene_to_group, target_label, norm_label, outdir):
    n_genes = xprs.shape[0]
    n_pairs = n_genes * (n_genes - 1) / 2 + n_genes
    for i, method in tqdm(enumerate(['pearson', 'spearman', pearson_pvalue, spearman_pvalue])):
        if i == 0 or i == 2:
            continue
        
        figure_table = xprs.T.corr(method=method).stack()
        figure_table.index.names = ['0', '1']
        figure_table = figure_table.reset_index()
        figure_table.columns = ['xname', 'yname', 'target']
        if i <= 1:
            figure_table['target'] = np.round(figure_table['target'], 2)
            method_label = method.capitalize()
        else:
            figure_table['target'] = figure_table['target'] * n_pairs
            figure_table.loc[:, 'target'].loc[figure_table['target'] > 1] = 1
            if i == 2:
                method_label = 'Pearson Pvalue'
            else:
                method_label = 'Spearman Pvalue'

        figure_table.loc[:, 'group_x'] = gene_to_group[figure_table.xname].values
        figure_table.loc[:, 'group_y'] = gene_to_group[figure_table.yname].values
        if i <= 1:
            mapper = LinearColorMapper(palette=list(Magma[256])[256:128:-1], low=0, high=1)
        else:
            mapper = LinearColorMapper(palette=list(Magma[256])[128:256], low=0, high=0.05)

        x_data = [(data.group_x, data.xname) for i, data in figure_table.iterrows()]
        x_uniq = [(data.group_y, data.yname) for i, data in figure_table.loc[figure_table.xname == figure_table.xname[0]].iterrows()]
        y_data = [(data.group_y, data.yname) for i, data in figure_table.iterrows()]
        num_genes = len(x_uniq)

        title = ''.join((
            f'{target_label.capitalize()} Lowly Expressed Genes {method_label} ',
            f'Correlations ({norm_label.capitalize()})'
        ))

        data_source = ColumnDataSource(
            data=dict(
                xname=x_data,
                yname=y_data,
                target=figure_table.target
            )
        )

        fig = figure(
            title='',
            x_axis_location="below",
            tools=['hover', 'save'],
            x_range=FactorRange(factors=x_uniq),
            y_range=FactorRange(factors=x_uniq),
            tooltips=[
                ('Gene X', '@xname'),
                ('Gene Y', '@yname'),
                (method_label, '@target')
            ]
        )

        figsize = max(num_genes * 3 + int(num_genes / 24), 1024)
        font_size_large = max(int(num_genes / 12), 28)
        font_size_medium = max(int(num_genes / 24), 14)

        fig.plot_width = figsize
        fig.plot_height = figsize
        fig.grid.grid_line_color = None
        fig.axis.axis_line_color = None
        fig.axis.major_tick_line_color = None
        fig.axis.major_label_text_font_size = f'{font_size_medium}px'  #'0.05vh'
        fig.axis.major_label_standoff = font_size_medium
        fig.axis.major_label_text_alpha = 0
        fig.xaxis.major_label_orientation = np.pi / 2
        fig.axis.group_text_font_size = f'{font_size_large}px' #'2.2vh'
        fig.axis.group_text_align = 'center'
        fig.title.text_font_size = f'{font_size_large}px' #'2.5vh'

        fig.rect(
            x='xname',
            y='yname',
            width=1.0,
            height=1.0,
            source=data_source,
            line_color=None,
            fill_color={'field': 'target', 'transform': mapper},
            hover_line_color='black'
        )

        color_bar = ColorBar(
            color_mapper=mapper,
            width=int(font_size_large),
            major_label_text_font_size= f'{font_size_large}px',
            label_standoff=int(font_size_large),
            border_line_color=None,
            location=(0, 0)
        )

        fig.add_layout(color_bar, 'right')

        file_name = os.path.join(
            outdir,
            f"{method_label.lower().replace(' ', '_')}_heatmap_{norm_label.lower().replace(' ', '_')}.html"
        )
        output_file(file_name, title=title)

        save(fig)

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
    fill_color = list(map(lambda x: tuple(int(x.lstrip('#')[i:i+2], 16) for i in (0, 2, 4)), list(fill_color)))

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
    save(fig)

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

    save(fig)

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
    save(fig)
    
    return

def compute_correlation_coefficients(data, method='spearman', stack=True):
    corr = data.T.corr(method=method).fillna(0)
    if stack:
        index = np.triu(np.ones(corr.shape)).astype(bool)
        index[np.diag_indices_from(index)] = False
        corr = corr.where(index).stack()
        corr.columns = ['gene_a', 'gene_b', method]

    return corr

def extract_tissue_exclusive_genes(
    xprs,
    tissues,
    xprs_high_threshold=10,
    xprs_low_threshold=1,
    #num_high_threshold=1000,
    #num_low_threshold=5
):
    index = {}
    gene_count = []
    for tissue_a in tissues.unique():
        xprs_tissue_a = xprs.loc[:, (tissues == tissue_a)]
        gene_filter_a = xprs_tissue_a.quantile(0.5, axis=1) >= xprs_high_threshold
        intersection = gene_filter_a.copy()
        for tissue_b in tissues.unique():
            if tissue_a == tissue_b:
                continue
            xprs_tissue_b = xprs.loc[:, (tissues == tissue_b)]
            gene_filter_b = xprs_tissue_b.quantile(0.5, axis=1) <= xprs_low_threshold
            intersection = intersection & gene_filter_b
        
        gene_count.append(intersection.sum())

        #if num_high_threshold >= intersection.sum() >= num_low_threshold:
        filtered = xprs.loc[intersection]
        index[tissue_a] = list(filtered.index)
        
    gene_count = pd.Series(gene_count, index=tissues.unique())

    return index, gene_count

def compute_rmse(true, pred):
    return np.sqrt(((true - pred) ** 2).mean())

def rename_tissue(tissue_to_gene, tissue_exclusive_count):
    tissue_to_gene_rename = {}
    for tissue, genes_of_tissue in tissue_to_gene.items():
        tissue_name = tissue.capitalize()
        if tissue_name == 'Cells - ebv-transformed lymphocytes':
            tissue_name = 'LCL'
        if tissue_name == 'Brain - cerebellum':
            tissue_name = 'Cerebellum'
        if tissue_name == 'Embryonic facial prominence':
            tissue_name = 'EF'
        if tissue_name == 'Stomach':
            tissue_name = 'Stom.'
        if tissue_name == 'Stomach':
            tissue_name = 'Stom.'
        if tissue_name == 'Intestine':
            tissue_name = 'Inst.'

        tissue_to_gene_rename[tissue_name] = genes_of_tissue
    index = []
    for tissue in tissue_exclusive_count.index:
        tissue_name = tissue.capitalize()
        if tissue_name == 'Cells - ebv-transformed lymphocytes':
            tissue_name = 'LCL'
        if tissue_name == 'Brain - cerebellum':
            tissue_name = 'Cerebellum'
        if tissue_name == 'Embryonic facial prominence':
            tissue_name = 'EF'
        if tissue_name == 'Stomach':
            tissue_name = 'Stom.'
        if tissue_name == 'Stomach':
            tissue_name = 'Stom.'
        if tissue_name == 'Intestine':
            tissue_name = 'Inst.'
        index.append(tissue_name)
    tissue_exclusive_count.index = index

    return tissue_to_gene_rename, tissue_exclusive_count

def sort_xprs_samples(xprs, tissues, tissue_exclusive_genes):
    # TODO: archive this function
    tissue_of_gene = []
    sample_order = []
    for tissue, genes_of_tissue in tissue_exclusive_genes.items():
        tissue_name = tissue.capitalize()
        if tissue_name == 'Cells - EBV-transformed lymphocytes':
            tissue_name = 'LCL'
        if tissue_name == 'Brain - Cerebellum':
            tissue_name = 'Cerebellum'
        if tissue_name == 'Embryonic facial prominence':
            tissue_name = 'EFP'
        
        tissue_of_gene.append(pd.Series(
            [tissue_name] * len(genes_of_tissue), 
            index=genes_of_tissue
        ))
        samples_of_tissue = (tissues == tissue)
        xprs_by_tissue = xprs.loc[genes_of_tissue, samples_of_tissue]
        xprs_sorted_samples = xprs_by_tissue.iloc[:, xprs_by_tissue.max().argsort()]
        sample_order.extend(xprs_sorted_samples.columns.values)

    tissue_of_gene = pd.concat(tissue_of_gene)
    
    return tissue_of_gene, sample_order
# %%
