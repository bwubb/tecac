#Brad Wubbenhorst

import argparse
import csv
import json
from collections import defaultdict
from datetime import datetime
import os

def read_variant_csv(csv_file):
    """Read variant CSV and extract key metrics"""
    variants=[]
    filter_counts=defaultdict(int)
    lof_by_filter=defaultdict(lambda: defaultdict(int))
    consequence_counts=defaultdict(int)
    impact_counts=defaultdict(int)
    gene_counts=defaultdict(int)
    variant_type_counts=defaultdict(int)
    
    with open(csv_file,'r') as f:
        reader=csv.DictReader(f)
        for row in reader:
            variants.append(row)
            
            # Count FILTER values
            filter_val=row.get('FILTER','.')
            filter_counts[filter_val]+=1
            
            # Count LoF levels by filter
            lof_level=row.get('Variant.LoF_level','.')
            lof_by_filter[filter_val][lof_level]+=1
            
            # Count consequences
            consequence=row.get('Variant.Consequence','.')
            consequence_counts[consequence]+=1
            
            # Count impacts
            impact=row.get('IMPACT','.')
            impact_counts[impact]+=1
            
            # Count genes
            gene=row.get('Gene','.')
            if gene!='.':
                gene_counts[gene]+=1
            
            # Count variant types (SNV vs indel)
            ref_len=len(row.get('REF',''))
            alt_len=len(row.get('ALT',''))
            if ref_len==1 and alt_len==1:
                variant_type_counts['SNV']+=1
            else:
                variant_type_counts['INDEL']+=1
    
    return {
        'variants':variants,
        'filter_counts':dict(filter_counts),
        'lof_by_filter':dict(lof_by_filter),
        'consequence_counts':dict(consequence_counts),
        'impact_counts':dict(impact_counts),
        'gene_counts':dict(gene_counts),
        'variant_type_counts':dict(variant_type_counts)
    }

def calculate_lof_loss_analysis(lof_by_filter):
    """Calculate how many LoF variants would be lost if only keeping PASS"""
    total_lof1=0
    total_lof2=0
    pass_lof1=0
    pass_lof2=0
    
    for filter_val,lof_counts in lof_by_filter.items():
        lof1_count=lof_counts.get('1',0)
        lof2_count=lof_counts.get('2',0)
        
        total_lof1+=lof1_count
        total_lof2+=lof2_count
        
        if filter_val=='PASS':
            pass_lof1=lof1_count
            pass_lof2=lof2_count
    
    lost_lof1=total_lof1-pass_lof1
    lost_lof2=total_lof2-pass_lof2
    
    return {
        'total_lof1':total_lof1,
        'total_lof2':total_lof2,
        'pass_lof1':pass_lof1,
        'pass_lof2':pass_lof2,
        'lost_lof1':lost_lof1,
        'lost_lof2':lost_lof2
    }

def generate_html_report(data,argv):
    """Generate comprehensive HTML variant QC report"""
    
    lof_analysis=calculate_lof_loss_analysis(data['lof_by_filter'])
    total_variants=len(data['variants'])
    
    # Get top genes
    top_genes=sorted(data['gene_counts'].items(),key=lambda x: x[1],reverse=True)[:20]
    
    # Get top consequences
    top_consequences=sorted(data['consequence_counts'].items(),key=lambda x: x[1],reverse=True)[:10]
    
    # Calculate percentages for LoF analysis
    lost_lof1_pct = lof_analysis['lost_lof1'] / lof_analysis['total_lof1'] * 100 if lof_analysis['total_lof1'] > 0 else 0
    lost_lof2_pct = lof_analysis['lost_lof2'] / lof_analysis['total_lof2'] * 100 if lof_analysis['total_lof2'] > 0 else 0
    
    html_content=f"""
<!DOCTYPE html>
<html>
<head>
    <title>Variant QC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
        .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
        .summary {{ background-color: #e8f4f8; padding: 15px; border-radius: 5px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .stat-box {{ background-color: #f8f9fa; padding: 15px; border-radius: 5px; text-align: center; }}
        .stat-number {{ font-size: 2em; font-weight: bold; color: #007bff; }}
        .stat-label {{ color: #6c757d; margin-top: 5px; }}
        .warning-box {{ background-color: #fff3cd; border: 1px solid #ffeaa7; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .success-box {{ background-color: #d4edda; border: 1px solid #c3e6cb; padding: 15px; border-radius: 5px; margin: 10px 0; }}
        .table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        .table th, .table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        .table th {{ background-color: #f2f2f2; }}
        .lof-highlight {{ background-color: #ffebee; }}
        .pass-highlight {{ background-color: #e8f5e8; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Variant Quality Control Report</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Input file:</strong> {argv['input_csv']}</p>
        <p><strong>Total variants:</strong> {total_variants:,}</p>
    </div>

    <div class="section">
        <h2>Summary Statistics</h2>
        <div class="stats-grid">
            <div class="stat-box">
                <div class="stat-number">{total_variants:,}</div>
                <div class="stat-label">Total Variants</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{data['filter_counts'].get('PASS',0):,}</div>
                <div class="stat-label">PASS Variants</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{lof_analysis['total_lof1']:,}</div>
                <div class="stat-label">LoF Level 1</div>
            </div>
            <div class="stat-box">
                <div class="stat-number">{lof_analysis['total_lof2']:,}</div>
                <div class="stat-label">LoF Level 2</div>
            </div>
        </div>
    </div>

    <div class="section">
        <h2>FILTER Field Analysis</h2>
        <div class="summary">
            <h3>Filter Distribution</h3>
            <table class="table">
                <tr>
                    <th>Filter</th>
                    <th>Count</th>
                    <th>Percentage</th>
                    <th>LoF Level 1</th>
                    <th>LoF Level 2</th>
                </tr>
"""
    
    for filter_val,count in sorted(data['filter_counts'].items(),key=lambda x: x[1],reverse=True):
        percentage=count/total_variants*100
        lof1_count=data['lof_by_filter'][filter_val].get('1',0)
        lof2_count=data['lof_by_filter'][filter_val].get('2',0)
        row_class='pass-highlight' if filter_val=='PASS' else ''
        
        html_content+=f"""
                <tr class="{row_class}">
                    <td>{filter_val}</td>
                    <td>{count:,}</td>
                    <td>{percentage:.1f}%</td>
                    <td>{lof1_count:,}</td>
                    <td>{lof2_count:,}</td>
                </tr>
"""
    
    html_content+=f"""
            </table>
        </div>
    </div>

    <div class="section">
        <h2>Loss of Function Analysis</h2>
        <div class="summary">
            <h3>Impact of PASS-only Filtering</h3>
            <div class="warning-box">
                <strong>⚠️ LoF variants that would be lost if only keeping PASS:</strong>
                <ul>
                    <li><strong>LoF Level 1:</strong> {lof_analysis['lost_lof1']:,} variants ({lost_lof1_pct:.1f}% of total LoF1)</li>
                    <li><strong>LoF Level 2:</strong> {lof_analysis['lost_lof2']:,} variants ({lost_lof2_pct:.1f}% of total LoF2)</li>
                </ul>
            </div>
            <div class="success-box">
                <strong>✅ LoF variants that would be retained with PASS-only:</strong>
                <ul>
                    <li><strong>LoF Level 1:</strong> {lof_analysis['pass_lof1']:,} variants</li>
                    <li><strong>LoF Level 2:</strong> {lof_analysis['pass_lof2']:,} variants</li>
                </ul>
            </div>
        </div>
    </div>

    <div class="section">
        <h2>Variant Type Distribution</h2>
        <div class="summary">
            <table class="table">
                <tr>
                    <th>Type</th>
                    <th>Count</th>
                    <th>Percentage</th>
                </tr>
"""
    
    for var_type,count in data['variant_type_counts'].items():
        percentage=count/total_variants*100
        html_content+=f"""
                <tr>
                    <td>{var_type}</td>
                    <td>{count:,}</td>
                    <td>{percentage:.1f}%</td>
                </tr>
"""
    
    html_content+="""
            </table>
        </div>
    </div>

    <div class="section">
        <h2>Impact Distribution</h2>
        <div class="summary">
            <table class="table">
                <tr>
                    <th>Impact</th>
                    <th>Count</th>
                    <th>Percentage</th>
                </tr>
"""
    
    for impact,count in sorted(data['impact_counts'].items(),key=lambda x: x[1],reverse=True):
        percentage=count/total_variants*100
        html_content+=f"""
                <tr>
                    <td>{impact}</td>
                    <td>{count:,}</td>
                    <td>{percentage:.1f}%</td>
                </tr>
"""
    
    html_content+="""
            </table>
        </div>
    </div>

    <div class="section">
        <h2>Top Genes</h2>
        <div class="summary">
            <table class="table">
                <tr>
                    <th>Gene</th>
                    <th>Variant Count</th>
                    <th>Percentage</th>
                </tr>
"""
    
    for gene,count in top_genes:
        percentage=count/total_variants*100
        html_content+=f"""
                <tr>
                    <td>{gene}</td>
                    <td>{count:,}</td>
                    <td>{percentage:.1f}%</td>
                </tr>
"""
    
    html_content+="""
            </table>
        </div>
    </div>

    <div class="section">
        <h2>Top Consequences</h2>
        <div class="summary">
            <table class="table">
                <tr>
                    <th>Consequence</th>
                    <th>Count</th>
                    <th>Percentage</th>
                </tr>
"""
    
    for consequence,count in top_consequences:
        percentage=count/total_variants*100
        html_content+=f"""
                <tr>
                    <td>{consequence}</td>
                    <td>{count:,}</td>
                    <td>{percentage:.1f}%</td>
                </tr>
"""
    
    html_content+=f"""
            </table>
        </div>
    </div>

    <div class="section">
        <h2>Data Summary</h2>
        <div class="summary">
            <h3>Analysis Details</h3>
            <ul>
                <li><strong>Total variants analyzed:</strong> {total_variants:,}</li>
                <li><strong>PASS variants:</strong> {data['filter_counts'].get('PASS',0):,} ({data['filter_counts'].get('PASS',0)/total_variants*100:.1f}%)</li>
                <li><strong>Non-PASS variants:</strong> {total_variants-data['filter_counts'].get('PASS',0):,} ({(total_variants-data['filter_counts'].get('PASS',0))/total_variants*100:.1f}%)</li>
                <li><strong>LoF Level 1 variants:</strong> {lof_analysis['total_lof1']:,}</li>
                <li><strong>LoF Level 2 variants:</strong> {lof_analysis['total_lof2']:,}</li>
                <li><strong>Unique genes:</strong> {len(data['gene_counts']):,}</li>
            </ul>
        </div>
    </div>
</body>
</html>
"""
    
    with open('variant_qc_report.html','w') as f:
        f.write(html_content)
    
    print("Variant QC report generated: variant_qc_report.html")

def write_json_output(data):
    """Write JSON output for programmatic access"""
    with open('variant_qc_stats.json','w') as f:
        json.dump(data,f,indent=4)
    print("JSON stats written: variant_qc_stats.json")

def get_args():
    p=argparse.ArgumentParser(description='Variant QC analysis for VEP parser output')
    p.add_argument('input_csv',help='Input CSV file from VEP parser')
    p.add_argument('--json',action='store_true',help='Also output JSON stats file')
    argv=p.parse_args()
    return vars(argv)

def main(argv=None):
    if not os.path.exists(argv['input_csv']):
        print(f"Error: Input file {argv['input_csv']} not found")
        return
    
    print(f"Processing variant QC for: {argv['input_csv']}")
    
    # Read and analyze data
    data=read_variant_csv(argv['input_csv'])
    
    # Generate HTML report
    generate_html_report(data,argv)
    
    # Write JSON if requested
    if argv.get('json'):
        write_json_output(data)
    
    print("Variant QC analysis complete.")

if __name__=='__main__':
    main(get_args()) 