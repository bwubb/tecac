import csv
import argparse
import re
import sys
from collections import defaultdict

def extract_aa_residue(hgvsp):
    """Extract amino acid residue from HGVSp notation (e.g., p.Glu415Lys -> Glu415, p.Ala259= -> Ala259)"""
    if not hgvsp or hgvsp=='.':
        return None
    # Match patterns like p.Glu415Lys, p.Ala259=, p.Pro2Thr
    match=re.search(r'p\.([A-Z][a-z]{2})(\d+)',hgvsp)
    if match:
        aa_code=match.group(1)
        residue_num=match.group(2)
        return f"{aa_code}{residue_num}"
    return None

def extract_c_position(hgvsc):
    """Extract cDNA position number from HGVSc notation (e.g., c.4C>A -> 4, c.1243G>A -> 1243)"""
    if not hgvsc or hgvsc=='.':
        return None
    # Match patterns like c.4C>A, c.1243G>A, c.-5G>A, c.*5G>A
    match=re.search(r'c\.(-?\d+|\*?\d+)',hgvsc)
    if match:
        try:
            pos_str=match.group(1)
            # Handle * notation (UTR positions)
            if pos_str.startswith('*'):
                return int(pos_str[1:])
            return int(pos_str)
        except ValueError:
            return None
    return None

def get_alt_frac(value):
    """Convert alt fraction string to float, return None if invalid"""
    if not value or value=='.':
        return None
    try:
        return float(value)
    except (ValueError,TypeError):
        return None

def get_depth(value):
    """Convert depth string to int, return None if invalid"""
    if not value or value=='.':
        return None
    try:
        return int(float(value))
    except (ValueError,TypeError):
        return None

def are_alt_fracs_close(frac1,frac2,depth1=None,depth2=None,threshold=0.10,min_depth=20):
    """
    Assess if two alt fractions are close enough to suggest cis (same chromosome).
    Returns confidence metrics for manual filtering - does NOT filter variants.
    
    Args:
        frac1: First alt fraction (float or None)
        frac2: Second alt fraction (float or None)
        depth1: First variant depth (int or None)
        depth2: Second variant depth (int or None)
        threshold: Maximum absolute difference to consider "close" (default 0.10)
        min_depth: Minimum depth for high confidence (default 20)
    
    Returns:
        tuple: (is_close: bool, difference: float or None, confidence: str, reason: str)
    """
    if frac1 is None or frac2 is None:
        return (False,None,"Low","Missing alt fraction data")
    
    diff=abs(frac1-frac2)
    
    # Assess confidence based on depth
    has_sufficient_depth=(depth1 is not None and depth2 is not None and 
                         depth1>=min_depth and depth2>=min_depth)
    
    if not has_sufficient_depth:
        # Low depth - lower confidence
        if diff<=threshold:
            return (True,diff,"Low",f"Alt fractions close but low depth (depth1={depth1},depth2={depth2})")
        else:
            return (False,diff,"Low",f"Alt fractions differ, low depth (depth1={depth1},depth2={depth2})")
    
    # Sufficient depth - higher confidence
    if diff<=threshold:
        return (True,diff,"High","Alt fractions are close (likely cis)")
    else:
        return (False,diff,"Medium","Alt fractions differ (possibly trans)")

def find_mnp_variants(input_csv,output_file,mode,regions_path=None,id_file_path=None):
    """Find variants affecting the same amino acid residue within 1-3 bases"""
    variant_groups=defaultdict(list)
    
    # Read input CSV and detect column structure
    with open(input_csv,'r') as infile:
        reader=csv.DictReader(infile)
        fieldnames=reader.fieldnames
        
        has_sample_id='Sample.ID' in fieldnames
        has_tumor_id='Tumor.ID' in fieldnames
        has_normal_id='Normal.ID' in fieldnames
        
        # Determine mode if not specified
        if mode=='auto':
            if has_sample_id:
                mode='cohort'
            elif has_tumor_id and has_normal_id:
                mode='tumor_normal'
            else:
                mode='no_sample'
        
        for row in reader:
            gene=row.get('Gene','')
            hgvsp=row.get('HGVSp','')
            hgvsc=row.get('HGVSc','')
            variant_class=row.get('Variant.Class','')
            
            if not gene or gene=='.':
                continue
            
            # Restrict to SNV only
            if variant_class!='SNV':
                continue
            
            aa_residue=extract_aa_residue(hgvsp)
            c_pos=extract_c_position(hgvsc)
            
            if not aa_residue:
                continue
            
            # Build grouping key based on mode
            if mode=='cohort' and has_sample_id:
                sample_id=row.get('Sample.ID','')
                if not sample_id:
                    continue
                key=(sample_id,gene,aa_residue)
            elif mode=='tumor_normal' and has_tumor_id:
                tumor_id=row.get('Tumor.ID','')
                if not tumor_id:
                    continue
                # Group by tumor only - always looking for MNPs in tumor
                key=(tumor_id,gene,aa_residue)
            else:
                # no_sample mode - group by gene and residue only
                key=(gene,aa_residue)
            
            variant_groups[key].append({
                'row':row,
                'c_pos':c_pos,
                'aa_residue':aa_residue
            })
    
    # Find groups with multiple variants (potential MNPs)
    mnp_pairs=[]
    for key,variants in variant_groups.items():
        if len(variants)>1:
            for i in range(len(variants)):
                for j in range(i+1,len(variants)):
                    var1=variants[i]
                    var2=variants[j]
                    if var1['row'].get('Chr')==var2['row'].get('Chr') and var1['row'].get('Start')==var2['row'].get('Start'):
                        continue
                    # Check if cDNA positions are within 1-3 bases
                    distance=None
                    if var1['c_pos'] is not None and var2['c_pos'] is not None:
                        distance=abs(var1['c_pos']-var2['c_pos'])
                        if distance>3:
                            continue
                    
                    # Check alt fractions if available (for cis/trans determination)
                    # Only check if columns exist
                    var1_alt_frac=None
                    var2_alt_frac=None
                    var1_depth=None
                    var2_depth=None
                    var1_tumor_frac=None
                    var2_tumor_frac=None
                    var1_tumor_depth=None
                    var2_tumor_depth=None
                    sample_cis_check=(False,None,"Low","No sample data")
                    tumor_cis_check=(False,None,"Low","No sample data")
                    
                    if has_sample_id:
                        var1_alt_frac=get_alt_frac(var1['row'].get('Sample.AltFrac','.'))
                        var2_alt_frac=get_alt_frac(var2['row'].get('Sample.AltFrac','.'))
                        var1_depth=get_depth(var1['row'].get('Sample.Depth','.'))
                        var2_depth=get_depth(var2['row'].get('Sample.Depth','.'))
                        sample_cis_check=are_alt_fracs_close(var1_alt_frac,var2_alt_frac,var1_depth,var2_depth)
                    
                    if has_tumor_id:
                        var1_tumor_frac=get_alt_frac(var1['row'].get('Tumor.AltFrac','.'))
                        var2_tumor_frac=get_alt_frac(var2['row'].get('Tumor.AltFrac','.'))
                        var1_tumor_depth=get_depth(var1['row'].get('Tumor.Depth','.'))
                        var2_tumor_depth=get_depth(var2['row'].get('Tumor.Depth','.'))
                        tumor_cis_check=are_alt_fracs_close(var1_tumor_frac,var2_tumor_frac,var1_tumor_depth,var2_tumor_depth)
                    
                    mnp_pairs.append({
                        'var1':var1['row'],
                        'var2':var2['row'],
                        'distance':distance,
                        'var1_alt_frac':var1_alt_frac,
                        'var2_alt_frac':var2_alt_frac,
                        'var1_depth':var1_depth,
                        'var2_depth':var2_depth,
                        'var1_tumor_frac':var1_tumor_frac,
                        'var2_tumor_frac':var2_tumor_frac,
                        'var1_tumor_depth':var1_tumor_depth,
                        'var2_tumor_depth':var2_tumor_depth,
                        'sample_cis_likely':sample_cis_check[0],
                        'sample_alt_frac_diff':sample_cis_check[1],
                        'sample_confidence':sample_cis_check[2],
                        'tumor_cis_likely':tumor_cis_check[0],
                        'tumor_alt_frac_diff':tumor_cis_check[1],
                        'tumor_confidence':tumor_cis_check[2]
                    })
    
    # Build output columns based on mode
    base_columns=['Chr','Start','REF','ALT','FILTER','ID',
                 'Gene','Variant.LoF_level','Variant.Category','Variant.Class',
                 'Variant.Consequence','HGVSc','HGVSp']
    
    # Add cis/trans assessment columns for filtering
    cis_columns=['AltFrac_Diff','Cis_Likely','Confidence']
    
    if mode=='cohort' and has_sample_id:
        output_columns=['Sample.ID']+base_columns+['Sample.Zyg','Sample.Depth','Sample.AltDepth','Sample.AltFrac']+cis_columns
    elif mode=='tumor_normal' and has_tumor_id and has_normal_id:
        output_columns=['Tumor.ID','Normal.ID']+base_columns+['Tumor.Zyg','Tumor.Depth','Tumor.AltDepth','Tumor.AltFrac',
                                                              'Normal.Zyg','Normal.Depth','Normal.AltDepth','Normal.AltFrac',
                                                              'Tumor_AltFrac_Diff','Tumor_Cis_Likely','Tumor_Confidence']
    else:
        output_columns=base_columns
    
    # Write output
    if output_file:
        outfile=open(output_file,'w')
    else:
        outfile=sys.stdout

    if output_file and mnp_pairs:
        pairs_path=output_file.replace('.csv','')+'.pairs.txt'
        with open(pairs_path,'w') as pf:
            for pair in mnp_pairs:
                id1=(pair['var1'].get('ID') or '').strip()
                id2=(pair['var2'].get('ID') or '').strip()
                if id1 and id2:
                    pf.write(f"{id1}\t{id2}\n")
        print(f"Pairs written to {pairs_path}",file=sys.stderr)

    if regions_path or id_file_path:
        positions=set()
        ids=set()
        for pair in mnp_pairs:
            for var in (pair['var1'],pair['var2']):
                chrom=(var.get('Chr') or '').strip()
                start_s=(var.get('Start') or '').strip()
                id_=(var.get('ID') or '').strip()
                if chrom and start_s and id_:
                    try:
                        positions.add((chrom,int(start_s)))
                        ids.add(id_)
                    except ValueError:
                        pass
        if regions_path:
            with open(regions_path,'w') as f:
                for (chrom,pos) in sorted(positions):
                    f.write(f"{chrom}\t{pos}\n")
            print(f"Regions written to {regions_path}",file=sys.stderr)
        if id_file_path:
            with open(id_file_path,'w') as f:
                for id_ in sorted(ids):
                    f.write(f"{id_}\n")
            print(f"IDs written to {id_file_path}",file=sys.stderr)
    
    writer=csv.DictWriter(outfile,fieldnames=output_columns,delimiter=',',
                         restval='.',extrasaction='ignore',quoting=csv.QUOTE_NONNUMERIC,
                         dialect='excel')
    writer.writeheader()
    
    for pair in mnp_pairs:
        var1=pair['var1']
        var2=pair['var2']
        
        row1={col:var1.get(col,'.') for col in output_columns if col not in cis_columns+['Tumor_AltFrac_Diff','Tumor_Cis_Likely','Tumor_Confidence']}
        row2={col:var2.get(col,'.') for col in output_columns if col not in cis_columns+['Tumor_AltFrac_Diff','Tumor_Cis_Likely','Tumor_Confidence']}
        
        # Add cis/trans assessment columns for manual filtering
        if mode=='cohort' and has_sample_id:
            row1['AltFrac_Diff']=f"{pair['sample_alt_frac_diff']:.4f}" if pair['sample_alt_frac_diff'] is not None else '.'
            row1['Cis_Likely']='Yes' if pair['sample_cis_likely'] else 'No'
            row1['Confidence']=pair.get('sample_confidence','.')
            row2['AltFrac_Diff']=f"{pair['sample_alt_frac_diff']:.4f}" if pair['sample_alt_frac_diff'] is not None else '.'
            row2['Cis_Likely']='Yes' if pair['sample_cis_likely'] else 'No'
            row2['Confidence']=pair.get('sample_confidence','.')
        elif mode=='tumor_normal' and has_tumor_id and has_normal_id:
            row1['Tumor_AltFrac_Diff']=f"{pair['tumor_alt_frac_diff']:.4f}" if pair['tumor_alt_frac_diff'] is not None else '.'
            row1['Tumor_Cis_Likely']='Yes' if pair['tumor_cis_likely'] else 'No'
            row1['Tumor_Confidence']=pair.get('tumor_confidence','.')
            row2['Tumor_AltFrac_Diff']=f"{pair['tumor_alt_frac_diff']:.4f}" if pair['tumor_alt_frac_diff'] is not None else '.'
            row2['Tumor_Cis_Likely']='Yes' if pair['tumor_cis_likely'] else 'No'
            row2['Tumor_Confidence']=pair.get('tumor_confidence','.')
        
        writer.writerow(row1)
        writer.writerow(row2)
    
    if output_file:
        outfile.close()
        print(f"Found {len(mnp_pairs)} pairs of variants affecting the same amino acid residue",file=sys.stderr)
        print(f"Output written to {output_file}",file=sys.stderr)
    else:
        print(f"Found {len(mnp_pairs)} pairs of variants affecting the same amino acid residue",file=sys.stderr)

def main():
    parser=argparse.ArgumentParser(description='Find MNP variants affecting the same amino acid residue')
    parser.add_argument('-i','--input_csv',required=True,help='Input CSV file from vep_vcf_parser2.py')
    parser.add_argument('-o','--output_csv',required=False,default=None,help='Output CSV file with MNP variant pairs (default: stdout)')
    parser.add_argument('-m','--mode',default='auto',choices=['auto','cohort','tumor_normal','no_sample'],
                       help='Input file mode: auto (detect), cohort, tumor_normal, or no_sample')
    parser.add_argument('-r','--regions',default=None,help='Output regions file (chrom\\tpos per line)')
    parser.add_argument('-I','--id-file',dest='id_file',default=None,help='Output ID file (one ID per line)')
    args=parser.parse_args()

    find_mnp_variants(args.input_csv,args.output_csv,args.mode,regions_path=args.regions,id_file_path=args.id_file)

if __name__=='__main__':
    main()

