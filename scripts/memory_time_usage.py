# %%
import argparse
import os
import tracemalloc
import yaml

from bokeh.models import ColorBar, ColumnDataSource, HoverTool, LabelSet, Legend
from bokeh.models import LinearColorMapper, FactorRange
from bokeh.palettes import Magma, d3
from bokeh.plotting import figure, output_file, show
import pandas as pd

from caiman import Analysis

#%%
def run_analysis(filename):
    xprs = pd.read_csv(filename, sep='\t', index_col=0)
    data_memory = [];
    for sample_size in [1, 10, 50, 100, 150, 200, 250]:
        tracemalloc.start()
        xprs_target = xprs.iloc[0:30000, 0:10]
        tmp = []
        for i in range(0, sample_size):
            tmp.append(xprs_target.add_prefix(f'{i}_'))
        xprs_target = pd.concat(tmp, axis=1)
        data_memory.append([sample_size * 10, tracemalloc.get_traced_memory()[1] * 9.53 * 1e-7])
        tracemalloc.stop()
        analysis = Analysis(xprs_target, **{'index_col': 0, 'sep': '\t'})
        corrected = analysis.correct(
            method='filter',
            adaptive_num_components=False,
            verbose=False,
            monitor=True
        )
    data_memory = pd.DataFrame(data_memory)
    data_memory.columns = ['sample_size', 'peak_memory']
    return data_memory

#%%
def output_figures(outdir, data_memory):
    result = pd.read_csv('./tmp/correction.log', sep='\t', header=None)
    result = result.dropna()
    result.columns = ['event', 'task', 'elapsed_time', 'cpu_time', 'peak_memory', 'sample_size']

    fitting = result.loc[result.task == 'fitting']
    correct = result.loc[result.task == 'correct']

    fig = figure(
        plot_width=800, plot_height=800,
        title='',
        x_axis_location="below",
        y_range=(50, 6500),
        tools="hover, save",
        tooltips=[
            ('Number of flanking components', '$x'),
            ('Peak memory (MiB)', '$y'),
        ]
    )
    fig.circle(fitting.sample_size, fitting.peak_memory, color="navy", alpha=0.5, size=20)
    fig.line(fitting.sample_size, fitting.peak_memory, color="navy", alpha=0.5, line_width=5, legend_label=f'Fitting and correction')

    fig.circle(data_memory.sample_size, data_memory.peak_memory, color="crimson", alpha=0.5, size=20)
    fig.line(data_memory.sample_size, data_memory.peak_memory, color="crimson", alpha=0.5, line_width=5, legend_label=f'Loading of the dataset')

    fig.xaxis.axis_label = 'Number of flanking components'
    fig.yaxis.axis_label = 'Peak memory (MiB)'
    fig.xaxis.axis_label_text_font_size = '20pt'
    fig.yaxis.axis_label_text_font_size = '20pt'
    fig.axis.major_label_text_font_size = '12pt'
    fig.legend.location = 'top_left'
    fig.legend.label_text_font_size = '18pt'

    output_file(os.path.join(outdir, 'memory_usage.html'), title='Peak Memory (MiB)')
    show(fig)

    fig = figure(
        plot_width=800, plot_height=800,
        title='',
        x_axis_location="below",
        tools="hover, save",
        tooltips=[
            ('Number of samples', '$x'),
            ('Computational time (sec)', '$y'),
        ]
    )
    fig.circle(fitting.sample_size, fitting.elapsed_time, color="crimson", alpha=0.5, size=20)
    fig.circle(correct.sample_size, correct.elapsed_time, color="navy", alpha=0.5, size=20)
    fig.line(fitting.sample_size, fitting.elapsed_time, color="crimson", alpha=0.5, line_width=5, legend_label=f'Fitting')
    fig.line(correct.sample_size, correct.elapsed_time, color="navy", alpha=0.5, line_width=5, legend_label=f'Correction')
    fig.legend.location = 'top_left'
    fig.xaxis.axis_label = 'Number of samples'
    fig.yaxis.axis_label = 'Elapsed Time (sec)'
    fig.xaxis.axis_label_text_font_size = '20pt'
    fig.yaxis.axis_label_text_font_size = '20pt'
    fig.legend.label_text_font_size = '18pt'
    fig.axis.major_label_text_font_size = '12pt'

    output_file(os.path.join(outdir, 'elapsed_time.html'), title='Elapsed Time (sec)')
    show(fig)

    fig = figure(
        plot_width=800, plot_height=800,
        title='',
        x_axis_location="below",
        tools="hover, save",
        tooltips=[
            ('Number of samples', '$x'),
            ('Computational time (sec)', '$y'),
        ]
    )
    fig.circle(fitting.sample_size, fitting.cpu_time, color="crimson", alpha=0.5, size=20)
    fig.circle(correct.sample_size, correct.cpu_time, color="navy", alpha=0.5, size=20)
    fig.line(fitting.sample_size, fitting.cpu_time, color="crimson", alpha=0.5, line_width=5, legend_label=f'Fitting')
    fig.line(correct.sample_size, correct.cpu_time, color="navy", alpha=0.5, line_width=5, legend_label=f'Correction')
    fig.legend.location = 'top_left'
    fig.xaxis.axis_label = 'Number of samples'
    fig.yaxis.axis_label = 'CPU time usage (sec)'
    fig.xaxis.axis_label_text_font_size = '20pt'
    fig.yaxis.axis_label_text_font_size = '20pt'
    fig.legend.label_text_font_size = '18pt'
    fig.axis.major_label_text_font_size = '12pt'

    output_file(os.path.join(outdir, 'cpu_time_usage.html'), title='CPU Time Usage (sec)')
    show(fig)

def main():
    description = """Compute Memory and CPU Time Usage for CAIMAN"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-x', '--xprs',
        metavar='[path]',
        type=str,
        default=None,
        help='Path to the count data.',
        required=True,
    )
    parser.add_argument(
        '-c', '--config',
        metavar='[config.yaml]',
        type=str,
        default=None,
        help='Path to the config data.',
        required=True,
    )

    args = parser.parse_args()
    config = yaml.load(open(args.config, 'r'), Loader=yaml.FullLoader)

    os.makedirs(config['out_dir'], exist_ok=True)

    data_memory = run_analysis(args.xprs)
    output_figures(config['out_dir'], data_memory)

if __name__ == '__main__':
    main()