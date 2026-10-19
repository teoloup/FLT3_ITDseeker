"""Main report plots use final passed read assignments, never discovery clusters."""
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .call_plots import reported_calls


def final_support_lengths(validation, read_lengths, calls, wt_reference_length):
    """Join final reported/WT support to measured trimmed read lengths.

    Net read length is not a WT-alignment insertion event: trimming and all
    sequencing indels contribute. No missing value is replaced by a call size.
    """
    valid=validation['valid_support'].astype(str).str.lower().eq('true')
    keep=valid & validation['filter_reason'].eq('pass') & (
        (validation['support_call'].eq('ITD-supporting') & validation['ref_alias'].isin(calls)) |
        (validation['support_call'].eq('WT-supporting') & validation['ref_alias'].eq('WT')))
    frame=validation.loc[keep,['read_id','ref_alias']].copy()
    frame['read_id']=frame['read_id'].astype(str)
    if frame['read_id'].duplicated().any():
        raise ValueError('Final plot input has multiple passed assignments for a read')
    frame['read_length']=frame['read_id'].map(read_lengths)
    if frame['read_length'].isna().any():
        raise ValueError('Actual trimmed read lengths are missing for final supporting reads')
    frame['read_length']=frame['read_length'].astype(int)
    frame['net_read_length']=frame['read_length']-wt_reference_length
    frame['called_insertion_size']=frame['ref_alias'].map({a:c['itd_len'] for a,c in calls.items()}).fillna(0).astype(int)
    counts=frame.groupby('ref_alias')['read_id'].nunique()
    for alias,info in calls.items():
        if int(counts.get(alias,0)) != int(info['reads']):
            raise ValueError(f'Final plot count disagrees with report for {alias}')
    return frame.sort_values(['ref_alias','read_id']).reset_index(drop=True)


def _labels(calls):
    totals=Counter(info['itd_len'] for _,info in calls)
    seen=Counter()
    labels={}
    for alias,info in calls:
        size=info['itd_len']
        suffix=' '+chr(ord('A')+seen[size]) if totals[size]>1 else ''
        seen[size]+=1
        labels[alias]=f"{size:g} bp{suffix} | AF {info['af']*100:.3f}% | {info['reads']:,} reads"
    return labels


def _bins(values):
    if not len(values):
        return np.array([-.5,.5])
    return np.arange(np.floor(min(values))-.5,np.ceil(max(values))+1.5,1)


def _histogram(ax, values, bins, color, label, linestyle='-'):
    # Every histogram input is the measured value for an individual final read.
    ax.hist(values,bins=bins,color=color,histtype='stepfilled',alpha=.12)
    ax.hist(values,bins=bins,color=color,histtype='step',linewidth=1.8,
            linestyle=linestyle,label=label)


def _call_markers(ax,calls,offset):
    for size in sorted({info['itd_len'] for _,info in calls}):
        pos=offset+size
        ax.axvline(pos,color='#555555',ls=':',lw=1)
        ax.text(pos, .98,f'{size:g} bp',transform=ax.get_xaxis_transform(),
                rotation=90,ha='right',va='top',fontsize=7,
                bbox=dict(facecolor='white',alpha=.75,edgecolor='none',pad=1))


def plot_final_call_distributions(frame, validated, out_dir, sample_name,
                                  wt_reference_length=336):
    """Write actual final read-length and net-length profiles at report paths.

    File names retain compatibility with the report loader. The first is no
    longer a GMM plot; the genuine GMM diagnostic has its own earlier-stage file.
    """
    calls=reported_calls(validated)
    labels=_labels(calls)
    colors=dict(zip((a for a,_ in calls),sns.color_palette('colorblind',max(1,len(calls)))))
    profiles={alias:frame.loc[frame.ref_alias.eq(alias)] for alias,_ in calls}
    wt=frame.loc[frame.ref_alias.eq('WT'),'read_length'].to_numpy()
    out_dir=Path(out_dir)
    out_dir.mkdir(parents=True,exist_ok=True)
    paths={}

    fig,axes=plt.subplots(2,1,figsize=(12,8),dpi=150)
    full_bins=_bins(frame.read_length.to_numpy())
    if len(wt):
        _histogram(axes[0],wt,full_bins,'#888888',f'WT | n={len(wt):,}')
    itd_values=frame.loc[frame.ref_alias.ne('WT'),'read_length'].to_numpy()
    zoom_bins=_bins(itd_values)
    for index,(alias,_) in enumerate(calls):
        values=profiles[alias].read_length.to_numpy()
        style=['-','--','-.',':'][index%4]
        _histogram(axes[0],values,full_bins,colors[alias],labels[alias],style)
        _histogram(axes[1],values,zoom_bins,colors[alias],labels[alias],style)
    axes[0].axvline(wt_reference_length,color='black',ls=':',lw=1,
                    label=f'WT reference length: {wt_reference_length:g} bp')
    axes[0].set_title(f'{sample_name}: WT and reported ITDs')
    axes[1].set_title('ITD reads — zoom')
    for ax in axes:
        _call_markers(ax,calls,wt_reference_length)
        ax.set_xlabel('Trimmed read length (bp)')
        ax.set_ylabel('Reads per 1 bp')
        ax.grid(alpha=.2)
        if ax is axes[0] or calls:
            ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=7)
    if not calls:
        axes[1].text(.5,.5,'No ITD passed the reporting threshold',transform=axes[1].transAxes,
                     ha='center',va='center')
    fig.tight_layout()
    paths['read_length_plot']=str(out_dir/f'{sample_name}_itd_gmm_fit_plot.png')
    fig.savefig(paths['read_length_plot'],dpi=150,bbox_inches='tight')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(12,5),dpi=150)
    net_values=frame.loc[frame.ref_alias.ne('WT'),'net_read_length'].to_numpy()
    bins=_bins(net_values)
    for index,(alias,_) in enumerate(calls):
        _histogram(ax,profiles[alias].net_read_length.to_numpy(),bins,colors[alias],
                   labels[alias],['-','--','-.',':'][index%4])
    _call_markers(ax,calls,0)
    ax.set_title(f'{sample_name}: observed net read-length changes of final reported allele reads')
    ax.set_xlabel(f'Actual trimmed read length − WT reference length ({wt_reference_length:g} bp)')
    ax.set_ylabel('Final supporting reads per 1 bp bin')
    ax.grid(alpha=.2)
    if calls:
        ax.legend(loc='upper left',bbox_to_anchor=(1.01,1),fontsize=7)
    else:
        ax.text(.5,.5,'No ITD passed the reporting threshold',transform=ax.transAxes,ha='center',va='center')
    fig.text(.5,.01,'Dotted markers: exact called insertion sizes. Histograms: measured net read lengths, affected by trimming and sequencing indels.',
             ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1))
    paths['net_read_length_plot']=str(out_dir/f'{sample_name}_itd_size_distribution.png')
    fig.savefig(paths['net_read_length_plot'],dpi=150,bbox_inches='tight')
    plt.close(fig)
    paths['profiles']={alias:dict(reads=len(profile),called_size=int(validated[alias]['itd_len']),
        read_length_counts={int(k):int(v) for k,v in profile.read_length.value_counts().sort_index().items()},
        net_read_length_counts={int(k):int(v) for k,v in profile.net_read_length.value_counts().sort_index().items()})
        for alias,profile in profiles.items()}
    return paths
