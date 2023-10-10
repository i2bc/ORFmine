#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 11:25:53 2020

@author: christospapadopoulos
"""

import argparse
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path
import pkg_resources
import random
from scipy.stats import ks_2samp
from typing import List, Union

from orfmine import DOCKER_IMAGE
from orfmine.orfold.lib.utils import parse_orfold_tab, generate_the_plot
from orfmine.utilities.container import ContainerCLI, add_container_args


def get_args():
    """
    Returns:
        Parameters
    """
    parser = argparse.ArgumentParser(description='ORF Foldability Distribution Plot')

    parser.add_argument(
        "--tab", "-T",
        type=str,
        required=True, 
        nargs="*",
        help="Table file(s) of foldability scores calculated by ORFold"
    )

    parser.add_argument(
        "--labels", "-L",
        type=str,
        required=False,
        default=[],
        nargs="*",
        help="Label(s) to be used for the legend. Label is set according to the given tab file(s) basename by default."
    )

    parser.add_argument(
        "--out", "-O",
        required=False,
        nargs="?",
        default='./',
        type=str,
        help="Output directory ('./' by default)."
    )

    parser = add_container_args(parser=parser)

    args = parser.parse_args()

    return args



data_file_path = pkg_resources.resource_filename('orfmine.orfold', 'data/globular.tab')
REF_HCA, REF_IUPRED, REF_TANGO = parse_orfold_tab(tab=data_file_path)


cmap = plt.get_cmap('tab20')
COLORS_BANK = cmap([0., 0.05263158, 0.10526316, 0.15789474, 0.21052632,
                    0.26315789, 0.31578947, 0.36842105, 0.42105263, 0.47368421,
                    0.52631579, 0.57894737, 0.63157895, 0.68421053, 0.73684211,
                    0.78947368, 0.84210526, 0.89473684, 0.94736842, 1.])


def get_kolmogorv_smirnov_test(scores: List):
    KS = []

    try:
        for j in range(len(scores)):
            KS_tmp = []
            if j > 0:
                for t in range(1, 1000):
                    my_sample_values = []
                    my_sample_idx = random.sample(k=559, population=range(len(scores[j])))
                    for l in my_sample_idx :
                        my_sample_values.append(scores[j][l])

                    ks_test = ks_2samp(scores[0], my_sample_values)
                    KS_tmp.append(ks_test[1])

                mean_KS = sum(KS_tmp)/len(KS_tmp)
                if mean_KS > 0.05:
                    KS.append("")
                elif mean_KS <= 0.05 and mean_KS > 0.01:
                    KS.append("*")
                elif mean_KS <= 0.01 and mean_KS > 0.001:
                    KS.append("**")
                elif mean_KS <= 0.001:
                    KS.append("***")
                else:
                    KS.append("PROBLEM")
    except:
        print("Warning: something went wrong with the Kolmogorov-Smirnov test... more scored sequences might be needed.")

    return KS


def run_orfplot(tabfiles: List[Union[str, Path]], outpath: Union[str, Path], labels: List[str]=[]):

    scores_to_plot = [REF_HCA]  # to_plot = [dis_ref_hca , REF_HCA ,tra_ref_hca]
    colors = ['grey']  # colors = ['forestgreen','grey','mediumvioletred']
    bins = [17]  # bins = [10,17,20]
    ref_labels = ["Globular proteins"]  # ref_labels = ["Disorder","Foldable","Aggregation"]

    # set outpath
    out_path = Path(outpath)
    out_path.mkdir(exist_ok=True)

    # set basename based on the first file basename
    basename = str(Path(tabfiles[0]).stem)
    out_basename = out_path / basename

    for idx, _file in enumerate(tabfiles):
        scores_hca , scores_iupred , scores_tango = parse_orfold_tab(tab=_file)
        scores_to_plot.append(scores_hca)
        colors.append(list(COLORS_BANK[idx]))  # colors.append([random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)])

        if not labels:
            ref_labels.append(basename)
        else:
            ref_labels.append(labels[idx])

        bins.append(15)

    my_plot = generate_the_plot(to_plot=scores_to_plot, colors=colors, bins=bins, labels=ref_labels)

    # Add the statistical KS-test
    KS = get_kolmogorv_smirnov_test(scores=scores_to_plot)

    my_plot = plt.text(x=-9.8, y=0.37, s=str('\n'.join(KS)), verticalalignment='top', horizontalalignment='center', fontsize="large", linespacing=1.35)

    fig = my_plot.get_figure()
    fig.savefig(str(out_basename) + '.png', transparent=False)
    fig.savefig(str(out_basename) + '_alpha.png', transparent=True)
    fig.savefig(str(out_basename) + '.pdf', transparent=True, dpi=300)

    fig2, ax = plt.subplots(figsize=(6, 1))
    fig2.subplots_adjust(bottom=0.5)
    norm = mpl.colors.Normalize(vmin=-10, vmax=10)
    cb1 = mpl.colorbar.ColorbarBase(
        ax,
        cmap=mpl.cm.coolwarm,
        orientation='horizontal',
        norm=norm
    )

    fig2.savefig(str(out_basename) + '_scale.pdf', transparent=True, dpi=300)



def run_orfplot_containerized(parameters: argparse.Namespace):
    # instantiate containerCLI handler
    cli = ContainerCLI(
            input_args=["--tab"],
            output_arg="--out",
            args=parameters,
            image_base=DOCKER_IMAGE,
            prog="orfplot",
            container_type="docker" if parameters.docker else "singularity",
            dev_mode=parameters.dev,
            package_binding={"orfmine": "/home/orfuser/orfmine/orfmine"}
        )

    cli.show()
    if not parameters.dry_run:
        cli.run()


def main():

    parameters = get_args()

    if parameters.docker or parameters.singularity:
        run_orfplot_containerized(parameters=parameters)
    else:
        start_time = datetime.now()

        run_orfplot(tabfiles=parameters.tab, outpath=parameters.out, labels=parameters.labels)

        end_time = datetime.now()
        print('\nDuration: {}'.format(end_time - start_time))


if __name__ == "__main__":
    main()
