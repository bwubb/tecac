#Brad Wubbenhorst

#import scipy
import argparse
import json
import csv
import numpy
from scipy import stats
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
from matplotlib import pyplot as plt
import seaborn as sns
import os
from datetime import datetime



# SN, Summary numbers:
# SN    [2]id   [3]key  [4]value
#SN	0	number of SNPs:	36054
#SN	0	number of indels:	2990
def sn(data):
    try:
        if 'number of SNPs' in data[2]:
            return {'#SNVs':int(data[3])}
        elif 'number of indels' in data[2]:
            return {'#INDELs':int(data[3])}
        else:
            return{}
    except (IndexError,ValueError):
        return{}

# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
#TSTV	0	23581	12483	1.89	23579	12475	1.89
def tstv(data):
    try:
        return {'Ti:Tv':float(data[4])}
    except (IndexError,ValueError):
        return{}

# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
def depth(data):
    pass

# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.
# PSC   [2]id   [3]sample       [4]nRefHom      [5]nNonRefHom   [6]nHets        [7]nTransitions [8]nTransversions       [9]nIndels      [10]average depth        [11]nSingletons [12]nHapRef     [13]nHapAlt     [14]nMissing
#PSC	0	sample_id	10632	10117	15662	19044	6735	879	65.0	26658	0	0	1755
def psc(data):
    try:
        return {'SM':data[2],
                'nRefHom_SNVs':int(data[3]),
                'nNonRefHom_SNVs':int(data[4]),
                'nHET_SNVs':int(data[5]),
                'nTransitions':int(data[6]),
                'nTransversions':int(data[7])}
    except (IndexError,ValueError):
        return{}


# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets.
# PSI   [2]id   [3]sample       [4]in-frame     [5]out-frame    [6]not applicable       [7]out/(in+out) ratio   [8]nInsHets     [9]nDelHets  [10]nInsAltHoms [11]nDelAltHoms
#PSI	0	sample_id	0	0	0	0.00	245	259	203	176
def psi(data):
    try:
        return {'nInsHet':int(data[7]),'nDelHet':int(data[8]),'nHET_INDELs':int(data[7])+int(data[8]),'nInsHom':int(data[9]),'nDelHom':int(data[10]),'nHOM_INDELs':int(data[9])+int(data[10]),'nIns':int(data[7])+int(data[9]),'nDel':int(data[8])+int(data[10])}
    except (IndexError,ValueError):
        return{}

def next_line(data):
    return {}

def read_stats(file,pull_from):
    all_samples={}
    current_sample=None
    sample_data=defaultdict(dict)
    
    with open(file,'r') as stats_file:
        for line in stats_file:
            if line.startswith('#'):
                continue
            else:
                data=line.rstrip().split('\t')
                section=data[0]
                
                # Handle per-sample sections (PSC, PSI)
                if section in ['PSC','PSI']:
                    sample_id=data[2]
                    if sample_id not in sample_data:
                        sample_data[sample_id]={}
                    sample_data[sample_id].update(pull_from.get(section,next_line)(data))
                
                # Handle global sections (SN, TSTV)
                # elif section in ['SN','TSTV']:
                #     global_data=pull_from.get(section,next_line)(data)
                #     for sample_id in sample_data:
                #         sample_data[sample_id].update(global_data)
    
    # Calculate derived metrics for each sample
    for sample_id,metrics in sample_data.items():
        try:
            # Het:Hom ratio - SNVs only (INDELs commented out for now)
            metrics['Het:Hom']=(metrics.get('nHET_SNVs',0))/(metrics.get('nNonRefHom_SNVs',0))
            # Original calculation with INDELs:
            # metrics['Het:Hom']=(metrics.get('nHET_SNVs',0)+metrics.get('nHET_INDELs',0))/(metrics.get('nNonRefHom_SNVs',0)+metrics.get('nHOM_INDELs',0))
        except ZeroDivisionError:
            metrics['Het:Hom']=0
        try:
            metrics['Ins:Del']=metrics.get('nIns',0)/metrics.get('nDel',0)
        except ZeroDivisionError:
            metrics['Ins:Del']=0
        
        # Calculate per-sample SNV and INDEL counts from PSC/PSI sections
        metrics['#SNVs']=metrics.get('nNonRefHom_SNVs',0)+metrics.get('nHET_SNVs',0)
        metrics['#INDELs']=metrics.get('nHOM_INDELs',0)+metrics.get('nHET_INDELs',0)
        
        # Calculate per-sample Ti:Tv ratio from PSC section
        try:
            metrics['Ti:Tv']=metrics.get('nTransitions',0)/metrics.get('nTransversions',0)
        except ZeroDivisionError:
            metrics['Ti:Tv']=0
        
        all_samples[sample_id]=metrics
    
    return all_samples

def plot_values(label,values,mdn,MAD,outliers,mad_threshold,upper,lower,output_dir='vcf_stats'):
    plt.figure(figsize=(10,6))
    sns.displot(values,kind='kde',color='b',legend=False,warn_singular=False)
    plt.axvline(x=mdn,c='k',label=f'Median: {mdn:.2f}')
    plt.axvline(x=upper,c='tab:green',linestyle='dashed',label=f'+{mad_threshold} MAD: {upper:.2f}')
    plt.axvline(x=lower,c='tab:olive',linestyle='dashed',label=f'-{mad_threshold} MAD: {lower:.2f}')
    
    # Highlight outliers
    if outliers:
        plt.scatter(outliers,[0]*len(outliers),c='red',s=50,alpha=0.7,label=f'Outliers ({len(outliers)})')
    
    plt.title(f'{label} Distribution (n={len(values)})')
    plt.xlabel(label)
    plt.ylabel('Density')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Density_{label}.png',dpi=300,bbox_inches='tight')
    plt.close()

def plot_values_by_site(label,values_by_site,output_dir='vcf_stats'):
    plt.figure(figsize=(10,6))
    colors=['b','r','g','orange','purple','brown','pink','gray','olive','cyan']
    for i,(site,values) in enumerate(values_by_site.items()):
        if values:
            color=colors[i%len(colors)]
            sns.kdeplot(values,color=color,label=f'{site} (n={len(values)})',warn_singular=False)
    
    plt.title(f'{label} Distribution by Site')
    plt.xlabel(label)
    plt.ylabel('Density')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Density_{label}_by_site.png',dpi=300,bbox_inches='tight')
    plt.close()

def create_side_by_side_plots(label,values,mdn,MAD,outliers,mad_threshold,upper,lower,values_by_site,output_dir='vcf_stats'):
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(20,6))
    
    # Left plot - original with MAD
    sns.kdeplot(values,color='b',ax=ax1,warn_singular=False)
    ax1.axvline(x=mdn,c='k',label=f'Median: {mdn:.2f}')
    ax1.axvline(x=upper,c='tab:green',linestyle='dashed',label=f'+{mad_threshold} MAD: {upper:.2f}')
    ax1.axvline(x=lower,c='tab:olive',linestyle='dashed',label=f'-{mad_threshold} MAD: {lower:.2f}')
    
    if outliers:
        ax1.scatter(outliers,[0]*len(outliers),c='red',s=50,alpha=0.7,label=f'Outliers ({len(outliers)})')
    
    ax1.set_title(f'{label} Distribution (n={len(values)})')
    ax1.set_xlabel(label)
    ax1.set_ylabel('Density')
    ax1.legend()
    
    # Right plot - by site
    colors=['b','r','g','orange','purple','brown','pink','gray','olive','cyan']
    for i,(site,site_values) in enumerate(values_by_site.items()):
        if site_values:
            color=colors[i%len(colors)]
            sns.kdeplot(site_values,color=color,label=f'{site} (n={len(site_values)})',ax=ax2,warn_singular=False)
    
    ax2.set_title(f'{label} Distribution by Site')
    ax2.set_xlabel(label)
    ax2.set_ylabel('Density')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/Density_{label}_combined.png',dpi=300,bbox_inches='tight')
    plt.close()

def these_stats():
    return {'PSC':psc,'PSI':psi}
    # return {'SN':sn,'TSTV':tstv,'PSC':psc,'PSI':psi}

def process_input_files(argv):
    """Process all input files and collect sample data"""
    input_values={}
    pull_from=these_stats()
    
    for n,file in enumerate(argv['input_fp']):
        print(f"Processing file {n+1}: {file}")
        if not os.path.exists(file):
            print(f"Warning: File {file} does not exist, skipping")
            continue
        sample_data=read_stats(file,pull_from)
        input_values.update(sample_data)
    
    print(f"Total samples processed: {len(input_values)}")
    return input_values

def check_missing_metrics(input_values):
    """Check for missing metrics and warn"""
    metrics=['#SNVs','#INDELs','Ti:Tv','Het:Hom','Ins:Del']
    for sample_id,metrics_dict in input_values.items():
        for label in metrics:
            if label not in metrics_dict:
                print(f"Warning: {sample_id} is missing {label}")

def calculate_statistics(input_values,argv,metadata=None,output_dir=None):
    """Calculate statistics and identify outliers"""
    output_values=defaultdict(dict)
    outliers_by_metric=defaultdict(list)
    csv_header=['sample_id']
    metrics=['#SNVs','#INDELs','Ti:Tv','Het:Hom','Ins:Del']
    
    # Add raw count columns to CSV header
    raw_columns=['nRefHom_SNVs','nNonRefHom_SNVs','nHET_SNVs','nTransitions','nTransversions','nHOM_INDELs','nHET_INDELs','nIns','nDel']
    csv_header.extend(raw_columns)
    
    with open(f'{output_dir}/stats_filter.txt','w') as filter_file:
        for label in metrics:
            values=[input_values[k][label] for k in input_values.keys() if label in input_values[k]]
            if not values:
                print(f"Warning: No data for metric {label}")
                continue
                
            mdn=numpy.median(values)
            MAD=stats.median_abs_deviation(values)
            upper=mdn+(argv['mad_threshold']*MAD)
            lower=mdn-(argv['mad_threshold']*MAD)
            
            csv_header.extend([f'{label}',f'mdn_{label}',f'MAD_{label}',f'upper_{label}',f'lower_{label}'])
            
            # Find outliers (two-sided)
            outlier_values=[]
            for sample_id,value in [(k,input_values[k][label]) for k in input_values.keys() if label in input_values[k]]:
                output_values[sample_id].update({
                    'sample_id':sample_id,
                    f'{label}':value,
                    f'mdn_{label}':mdn,
                    f'MAD_{label}':MAD,
                    f'upper_{label}':upper,
                    f'lower_{label}':lower
                })
                
                if value>upper or value<lower:
                    outliers_by_metric[label].append(sample_id)
                    outlier_values.append(value)
                    filter_file.write(f'{sample_id} {label} {value:.2f} (upper: {upper:.2f}, lower: {lower:.2f})\n')
            
            # Organize values by site if metadata available
            values_by_site=defaultdict(list)
            if metadata:
                for sample_id,value in [(k,input_values[k][label]) for k in input_values.keys() if label in input_values[k]]:
                    if sample_id in metadata:
                        site=metadata[sample_id]['site']
                        values_by_site[site].append(value)
            
            # Generate plots - always create the original MAD plot
            plot_values(label,values,mdn,MAD,outlier_values,argv['mad_threshold'],upper,lower,output_dir)
            
            # Generate site-specific plot if metadata available
            if metadata and values_by_site:
                plot_values_by_site(label,values_by_site,output_dir)
    
    # Add raw values to output
    for sample_id,metrics_dict in input_values.items():
        for col in raw_columns:
            output_values[sample_id][col]=metrics_dict.get(col,0)
    
    return output_values,outliers_by_metric,csv_header

def write_output_files(output_values,csv_header,output_dir):
    """Write CSV and JSON output files"""
    with open(f'{output_dir}/qc_stats.json','w') as outfile:
        json.dump(dict(output_values),outfile,indent=4)

    with open(f'{output_dir}/qc_stats.csv','w') as outfile:
        writer=csv.DictWriter(outfile,delimiter=',',fieldnames=csv_header)
        writer.writeheader()
        for k,v in output_values.items():
            writer.writerow(v)

def print_summary(outliers_by_metric,all_samples):
    """Print summary statistics"""
    total_samples=len(all_samples)
    failed_samples=set()
    for outliers in outliers_by_metric.values():
        failed_samples.update(outliers)
    total_failed=len(failed_samples)
    total_passed=total_samples-total_failed
    
    print(f"\n=== QC Summary ===")
    print(f"Total samples processed: {total_samples}")
    print(f"Samples passed QC: {total_passed}")
    print(f"Samples failed QC: {total_failed}")
    print(f"Pass rate: {total_passed/total_samples*100:.1f}%")
    
    print(f"\n=== Failures by Metric ===")
    for metric,outliers in outliers_by_metric.items():
        print(f"{metric}: {len(outliers)} failures ({len(outliers)/total_samples*100:.1f}% of samples)")
    
    if failed_samples:
        print(f"\n=== Failed Samples ===")
        for sample_id in sorted(failed_samples):
            failure_reasons=[]
            for metric,outliers in outliers_by_metric.items():
                if sample_id in outliers:
                    failure_reasons.append(metric)
            print(f"{sample_id}: failed {', '.join(failure_reasons)}")
    else:
        print(f"\n=== All samples passed QC ===")



def generate_html_report(output_values,outliers_by_metric,argv,all_samples,metadata=None,output_dir=None):
    """Generate comprehensive HTML QC report with detailed failure tracking"""
    
    # Calculate summary statistics
    total_samples=len(all_samples)
    failed_samples=set()
    for outliers in outliers_by_metric.values():
        failed_samples.update(outliers)
    total_failed=len(failed_samples)
    total_passed=total_samples-total_failed
    
    # Calculate case/control and site breakdowns
    case_control_stats={'Case':{'total':0,'passed':0,'failed':0},'Control':{'total':0,'passed':0,'failed':0}}
    site_stats={}
    
    if metadata:
        for sample_id in all_samples:
            if sample_id in metadata:
                case_control=metadata[sample_id]['case_control']
                site=metadata[sample_id]['site']
                
                # Update case/control stats
                if case_control in case_control_stats:
                    case_control_stats[case_control]['total']+=1
                    if sample_id in failed_samples:
                        case_control_stats[case_control]['failed']+=1
                    else:
                        case_control_stats[case_control]['passed']+=1
                
                # Update site stats
                if site not in site_stats:
                    site_stats[site]={'total':0,'passed':0,'failed':0}
                site_stats[site]['total']+=1
                if sample_id in failed_samples:
                    site_stats[site]['failed']+=1
                else:
                    site_stats[site]['passed']+=1
    
    # Create detailed failure summary
    failure_summary={}
    for sample_id in failed_samples:
        failure_reasons=[]
        for metric,outliers in outliers_by_metric.items():
            if sample_id in outliers:
                value=output_values[sample_id].get(metric,0)
                upper=output_values[sample_id].get(f'upper_{metric}',0)
                lower=output_values[sample_id].get(f'lower_{metric}',0)
                # Format value properly based on type
                if isinstance(value,float):
                    value_str=f"{value:.3f}"
                else:
                    value_str=str(value)
                if isinstance(upper,float):
                    upper_str=f"{upper:.3f}"
                else:
                    upper_str=str(upper)
                if isinstance(lower,float):
                    lower_str=f"{lower:.3f}"
                else:
                    lower_str=str(lower)
                failure_reasons.append(f"{metric}: {value_str} (upper: {upper_str}, lower: {lower_str})")
        failure_summary[sample_id]=failure_reasons
    
    # Add metadata to failure summary if available
    if metadata:
        for sample_id in failure_summary:
            if sample_id in metadata:
                site=metadata[sample_id]['site']
                case_control=metadata[sample_id]['case_control']
                failure_summary[sample_id].insert(0,f"Site: {site}, Type: {case_control}")
    
    # Count failures by metric
    metric_failure_counts={}
    for metric,outliers in outliers_by_metric.items():
        metric_failure_counts[metric]=len(outliers)
    
    html_content=f"""
<!DOCTYPE html>
<html>
<head>
    <title>VCF QC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
        .metric {{ margin: 10px 0; padding: 10px; background-color: #f9f9f9; }}
        .outlier {{ color: red; font-weight: bold; }}
        .summary {{ background-color: #e8f4f8; padding: 15px; border-radius: 5px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .stat-box {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; text-align: center; }}
        .stat-number {{ font-size: 2em; font-weight: bold; color: #007bff; }}
        .stat-label {{ color: #6c757d; margin-top: 5px; }}
        .failure-box {{ background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .success-box {{ background-color: #d4edda; border: 1px solid #c3e6cb; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .plot {{ text-align: center; margin: 20px 0; }}
        .plot img {{ max-width: 600px; width: 100%; height: auto; border: 1px solid #ddd; }}
        .failure-reasons {{ margin-left: 20px; }}
        .failure-reasons li {{ margin: 5px 0; }}
        .metadata-info {{ font-style: italic; color: #666; }}
        .breakdown-table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        .breakdown-table th, .breakdown-table td {{ border: 1px solid #ddd; padding: 8px; text-align: center; }}
        .breakdown-table th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>VCF Quality Control Report</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Input files:</strong> {', '.join(argv['input_fp'])}</p>
        <p><strong>MAD threshold:</strong> {argv['mad_threshold']}</p>
        {f"<p><strong>Bed file:</strong> {argv.get('bed')}</p>" if argv.get('bed') else ""}
        {f"<p><strong>Metadata file:</strong> {argv.get('metadata')}</p>" if argv.get('metadata') else ""}
    </div>

    <div class="section">
        <h2>Summary Statistics</h2>
        <div class="stats-grid">
            <div class="stat-box">
                <div class="stat-number">{total_samples}</div>
                <div class="stat-label">Total Samples</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{total_passed}</div>
                <div class="stat-label">Passed QC</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{total_failed}</div>
                <div class="stat-label">Failed QC</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{total_passed/total_samples*100:.1f}%</div>
                <div class="stat-label">Pass Rate</div>
            </div>
        </div>
        
        <div class="summary">
            <h3>QC Results Overview</h3>
            <div class="success-box">
                <strong>✅ {total_passed} samples passed all QC metrics</strong>
            </div>
            <div class="failure-box">
                <strong>❌ {total_failed} samples failed one or more QC metrics</strong>
            </div>
"""
    
    # Add case/control breakdown if metadata available
    if metadata and any(case_control_stats[cat]['total']>0 for cat in case_control_stats):
        html_content+="""
            <h4>Case/Control Breakdown</h4>
            <table class="breakdown-table">
                <tr>
                    <th>Type</th>
                    <th>Total</th>
                    <th>Passed</th>
                    <th>Failed</th>
                    <th>Pass Rate</th>
                </tr>
"""
        for case_control,stats in case_control_stats.items():
            if stats['total']>0:
                pass_rate=stats['passed']/stats['total']*100
                html_content+=f"""
                <tr>
                    <td>{case_control}</td>
                    <td>{stats['total']}</td>
                    <td>{stats['passed']}</td>
                    <td>{stats['failed']}</td>
                    <td>{pass_rate:.1f}%</td>
                </tr>
"""
        html_content+="</table>\n"
    
    # Add site breakdown if metadata available
    if metadata and site_stats:
        html_content+="""
            <h4>Site Breakdown</h4>
            <table class="breakdown-table">
                <tr>
                    <th>Site</th>
                    <th>Total</th>
                    <th>Passed</th>
                    <th>Failed</th>
                    <th>Pass Rate</th>
                </tr>
"""
        for site,stats in sorted(site_stats.items()):
            pass_rate=stats['passed']/stats['total']*100
            html_content+=f"""
                <tr>
                    <td>{site}</td>
                    <td>{stats['total']}</td>
                    <td>{stats['passed']}</td>
                    <td>{stats['failed']}</td>
                    <td>{pass_rate:.1f}%</td>
                </tr>
"""
        html_content+="</table>\n"
    
    html_content+="""
        </div>
    </div>

    <div class="section">
        <h2>Distribution Plots</h2>
        <p>Distribution plots for each QC metric with outlier detection thresholds.</p>
"""
    
    # Add plots at the top
    for metric in ['#SNVs','#INDELs','Ti:Tv','Het:Hom','Ins:Del']:
        mad_image_path=f'{output_dir}/Density_{metric}.png'
        site_image_path=f'{output_dir}/Density_{metric}_by_site.png'
        
        mad_embedded=embed_image_as_base64(mad_image_path)
        site_embedded=embed_image_as_base64(site_image_path)
        
        html_content+=f"""
        <div class="metric">
            <h3>{metric} Distribution</h3>
            <div style="display: flex; justify-content: space-between; gap: 20px; align-items: flex-start;">
"""
        
        if metadata and site_embedded:
            html_content+=f"""
                <div style="flex: 1; max-width: 50%;">
                    <h4>Distribution by Site</h4>
                    <img src="{site_embedded}" alt="{metric} site distribution plot" style="width: 100%; height: auto;">
                </div>
"""
        
        html_content+=f"""
                <div style="flex: 1; max-width: 50%;">
                    <h4>MAD Outlier Detection</h4>
"""
        
        if mad_embedded:
            html_content+=f'<img src="{mad_embedded}" alt="{metric} MAD distribution plot" style="width: 100%; height: auto;">'
        else:
            html_content+=f'<p style="color: red;">Image not found: {mad_image_path}</p>'
        
        html_content+="""
                </div>
            </div>
        </div>
"""
    
    html_content+="""
    </div>

    <div class="section">
        <h2>Failure Analysis by Metric</h2>
        <div class="summary">
            <h3>Number of Failures per Metric</h3>
            <ul>
"""
    
    for metric,count in metric_failure_counts.items():
        html_content+=f"<li><strong>{metric}:</strong> {count} failures ({count/total_samples*100:.1f}% of samples)</li>\n"
    
    html_content+="""
            </ul>
        </div>
    </div>

    <div class="section">
        <h2>Data Sources and Metrics</h2>
        <div class="summary">
            <h3>Metrics Description</h3>
            <ul>
                <li><strong>#SNVs:</strong> Total number of single nucleotide variants per sample = nHOM_SNVs + nHET_SNVs (from PSC section)</li>
                <li><strong>#INDELs:</strong> Total number of insertions/deletions per sample = nHOM_INDELs + nHET_INDELs (from PSI section)</li>
                <li><strong>Ti:Tv:</strong> Transition/Transversion ratio per sample = nTransitions / nTransversions (from PSC section)</li>
                <li><strong>Het:Hom:</strong> Heterozygous to homozygous ratio = (nHET_SNVs + nHET_INDELs) / (nNonRefHom_SNVs + nHOM_INDELs)</li>
                <li><strong>Ins:Del:</strong> Insertion to deletion ratio = nIns / nDel (from PSI section)</li>
            </ul>
            <p><strong>Data parsing:</strong> Extracted from bcftools stats output using PSC (per-sample counts) and PSI (per-sample indels) sections</p>
            <p><strong>Outlier detection:</strong> Samples with values > median + {argv['mad_threshold']} × MAD are flagged as outliers</p>
        </div>
    </div>

    <div class="section">
        <h2>Failed Samples Details</h2>
"""
    
    if failure_summary:
        html_content+="<h3>Samples with QC Failures</h3>\n"
        for sample_id,reasons in failure_summary.items():
            html_content+=f"<div class='failure-box'>\n"
            html_content+=f"<strong>Sample: {sample_id}</strong>\n"
            html_content+=f"<div class='failure-reasons'>\n<ul>\n"
            for reason in reasons:
                if reason.startswith('Site:'):
                    html_content+=f"<li class='metadata-info'>{reason}</li>\n"
                else:
                    html_content+=f"<li>{reason}</li>\n"
            html_content+=f"</ul>\n</div>\n</div>\n"
    else:
        html_content+="<div class='success-box'><strong>No samples failed QC!</strong></div>\n"
    
    html_content+="""
    </div>
</body>
</html>
"""
    
    with open(f'{output_dir}/qc_report.html','w') as f:
        f.write(html_content)
    
    print("HTML report generated: vcf_stats/qc_report.html")

def generate_pdf_report(html_file):
    """Generate PDF from HTML report using wkhtmltopdf"""
    pdf_file=html_file.replace('.html','.pdf')
    
    try:
        import subprocess
        result=subprocess.run(['wkhtmltopdf',html_file,pdf_file],capture_output=True,text=True)
        if result.returncode==0:
            print(f"PDF report: {pdf_file}")
            return True
        else:
            print(f"wkhtmltopdf failed: {result.stderr}")
            return False
    except FileNotFoundError:
        print("wkhtmltopdf not found. Install with: conda install -c conda-forge wkhtmltopdf")
        return False
    except Exception as e:
        print(f"PDF generation failed: {e}")
        return False

def embed_image_as_base64(image_path):
    """Embed image as base64 data in HTML"""
    try:
        import base64
        with open(image_path,'rb') as img_file:
            img_data=base64.b64encode(img_file.read()).decode('utf-8')
        # Get file extension to determine MIME type
        ext=image_path.split('.')[-1].lower()
        mime_type=f"image/{ext}" if ext in ['png','jpg','jpeg','gif'] else "image/png"
        return f"data:{mime_type};base64,{img_data}"
    except Exception as e:
        print(f"Warning: Could not embed image {image_path}: {e}")
        return None

def load_metadata(metadata_file):
    """Load sample metadata from CSV file"""
    metadata={}
    if metadata_file and os.path.exists(metadata_file):
        with open(metadata_file,'r') as f:
            reader=csv.DictReader(f)
            for row in reader:
                sample_id=row.get('SampleID','')
                if sample_id:
                    status=row.get('Status','Unknown')
                    # Convert status to case/control
                    case_control='Case' if status=='1' else 'Control' if status=='0' else status
                    metadata[sample_id]={
                        'site':row.get('Site','Unknown'),
                        'case_control':case_control
                    }
        print(f"Loaded metadata for {len(metadata)} samples")
    return metadata

def get_args():
    p=argparse.ArgumentParser(description='QC analysis for multi-sample VCF stats')
    p.add_argument('-b','--bed',help='Bed file used in variant calling.')
    p.add_argument('-o','--output_prefix',default='vcf_stats',help='Output prefix for files and directories (default: vcf_stats)')
    p.add_argument('--mad-threshold',type=float,default=4.0,help='MAD multiplier for outlier detection (default: 4.0)')
    p.add_argument('--metadata',help='Sample metadata file (CSV with SampleID, Site, Status columns)')
    p.add_argument('--pdf',action='store_true',help='Generate PDF report (requires wkhtmltopdf to be installed)')
    p.add_argument('input_fp',nargs='+',help='One or more bcftools stats input files')
    argv=p.parse_args()
    return vars(argv)

def main(argv=None):
    if argv is None:
        argv=get_args()
    
    output_prefix=argv.get('output_prefix','vcf_stats')
    output_dir=output_prefix
    
    # Create output directory
    os.makedirs(output_dir,exist_ok=True)
    
    # Process input files
    all_samples=process_input_files(argv)
    
    # Check for missing metrics
    check_missing_metrics(all_samples)
    
    # Load metadata first
    metadata=load_metadata(argv.get('metadata'))
    
    # Calculate statistics and identify outliers
    output_values,outliers_by_metric,csv_header=calculate_statistics(all_samples,argv,metadata,output_dir)
    
    # Print summary
    print_summary(outliers_by_metric,all_samples)
    
    # Write output files
    write_output_files(output_values,csv_header,output_dir)
    
    # Generate HTML report
    generate_html_report(output_values,outliers_by_metric,argv,all_samples,metadata,output_dir)
    
    # Generate PDF report only if requested
    if argv.get('pdf'):
        html_file=f'{output_dir}/qc_report.html'
        pdf_success=generate_pdf_report(html_file)
        if not pdf_success:
            print("PDF generation failed - install wkhtmltopdf to generate PDFs")
    
    print(f"\nQC analysis complete. Results saved to {output_dir}/")
    print(f"HTML report: {output_dir}/qc_report.html")

if __name__=='__main__':
    main(get_args())
