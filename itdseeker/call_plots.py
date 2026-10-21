"""Final sequence-call markers shared by discovery diagnostic plots."""
import numpy as np
import seaborn as sns


def reported_calls(validated, reported_aliases=None):
    """Only reporting-filter survivors; stable order also separates same-size calls."""
    return [(str(alias), info) for alias, info in sorted(
        (validated or {}).items(), key=lambda item: (item[1]['itd_len'], str(item[0])))
        if str(alias).upper() != 'WT' and info['itd_len'] > 0
        and (reported_aliases is None or str(alias) in reported_aliases)]


def draw_call_panel(ax, calls, *, wt_length=0, full_allele=False):
    """One row per reported sequence, with no implied Gaussian or read counts."""
    colors = sns.color_palette('husl', max(len(calls), 1))
    ax.axvline(wt_length, color='black', lw=1, ls=':', label=f'WT reference: {wt_length:g} bp')
    labels, positions = [], []
    for index, ((alias, info), color) in enumerate(zip(calls, colors)):
        position = wt_length + info['itd_len']
        positions.append(position)
        ax.plot(position, index, marker='D', color=color, ms=7, linestyle='none')
        source = 'observed rescue' if alias.startswith('ITD_OBS_') else 'peak-derived'
        labels.append(f"{alias}\n{info['itd_len']:g} bp | AF {info['af']*100:.2f}% | {info['reads']:,} reads | {source}")
    ax.set_yticks(range(len(labels)), labels, fontsize=7.5)
    ax.set_ylim(max(len(labels)-.3, .7), -.7)
    if not calls:
        ax.text(.5,.5,'No ITD passed the reporting threshold',ha='center',va='center',transform=ax.transAxes)
    bounds = [wt_length] + positions
    padding = max((max(bounds)-min(bounds))*.07,5)
    ax.set_xlim(min(bounds)-padding,max(bounds)+padding)
    if full_allele:
        ax.set_title('Final reported alleles: WT reference length + insertion size (sequence calls, not GMM fits)',fontsize=9)
        ax.set_xlabel(f'Expected complete allele length (bp); WT reference = {wt_length:g} bp')
    else:
        ax.set_title('Final reported insertion sizes (sequence calls; markers are not histogram counts)',fontsize=9)
        ax.set_xlabel('Inserted bases relative to WT (bp); WT = 0 bp')
    ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=7)
    ax.grid(axis='x',alpha=.25)
    return np.asarray(positions)
