#import scipy
import argparse
import json
import csv
import numpy
import seaborn as sns
from scipy import stats
from collections import defaultdict
from matplotlib import pyplot as plt



# SN, Summary numbers:
# SN    [2]id   [3]key  [4]value
#SN	0	number of SNPs:	36054
#SN	0	number of indels:	2990
def sn(data):
    if 'number of SNPs' in data[2]:
        return {'#SNVs':int(data[3])}
    elif 'number of indels' in data[2]:
        return {'#INDELs':int(data[3])}
    else:
        return{}

# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
#TSTV	0	23581	12483	1.89	23579	12475	1.89
def tstv(data):
    return {'Ti:Tv':float(data[4])}

# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
def depth(data):
    pass

# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.
# PSC   [2]id   [3]sample       [4]nRefHom      [5]nNonRefHom   [6]nHets        [7]nTransitions [8]nTransversions       [9]nIndels  [10]average depth        [11]nSingletons [12]nHapRef     [13]nHapAlt     [14]nMissing
#PSC	0	sample_id	10632	10117	15662	19044	6735	879	65.0	26658	0	0	1755
def psc(data):
    return {'SM':data[2],'nHOM_SNVs':int(data[4]),'nHET_SNVs':int(data[5])}


# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets.
# PSI   [2]id   [3]sample       [4]in-frame     [5]out-frame    [6]not applicable       [7]out/(in+out) ratio   [8]nInsHets     [9]nDelHets  [10]nInsAltHoms [11]nDelAltHoms
#PSI	0	sample_id	0	0	0	0.00	245	259	203	176
def psi(data):
    return {'nInsHet':int(data[7]),'nDelHet':int(data[8]),'nHET_INDELs':int(data[7])+int(data[8]),'nInsHom':int(data[9]),'nDelHom':int(data[10]),'nHOM_INDELs':int(data[9])+int(data[10]),'nIns':int(data[7])+int(data[9]),'nDel':int(data[8])+int(data[10])}

def next_line(data):
    return {}

def read_stats(file,pull_from):
    with open(file,'r') as stats_file:
        sample_values={}
        for line in stats_file:
            if line.startswith('#'):
                continue
            else:
                data=line.rstrip().split('\t')
                sample_values.update(pull_from.get(data[0],next_line)(data))
    try:
        sample_values['Het:Hom']=(sample_values.get('nHET_SNVs',0)+sample_values.get('nHET_INDELs',0))/(sample_values.get('nHOM_SNVs',0)+sample_values.get('nHOM_INDELs',0))
    except ZeroDivisionError:
        sample_values['Het:Hom']=0
    try:
        sample_values['Ins:Del']=sample_values.get('nIns',0)/sample_values.get('nDel',0)
    except ZeroDivisionError:
        sample_values['Ins:Del']=0

    return {sample_values.get('SM','NULL'):sample_values}

def plot_values(label,values,mdn,MAD):
    #MAD=stats.median_abs_deviation(values)
    #mdn=numpy.median(values)
    sns.displot(values,kind='kde',color='b',legend=False)
    plt.axvline(x=mdn,c='k')
    plt.axvline(x=mdn+MAD,c='tab:green',linestyle='dashed')
    plt.axvline(x=mdn+(2*MAD),c='tab:olive',linestyle='dashed')
    plt.axvline(x=mdn+(3*MAD),c='tab:orange',linestyle='dashed')
    plt.axvline(x=mdn+(4*MAD),c='tab:red',linestyle='dashed')
    plt.title('Density Plot with MAD')
    plt.xlabel(label)
    plt.ylabel('Density')
    plt.savefig(f'vcf_stats/Density_{label}.png')

def these_stats():
    return {'SN':sn,'TSTV':tstv,'PSC':psc,'PSI':psi}

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('-o','--output_fp',help='Output file base')
    p.add_argument('input_fp',nargs=argparse.REMAINDER,help='One or more input files, but since you are calculating medians and sdv you should probably have a lot of files.')
    argv=p.parse_args()
    return vars(argv)

def main(argv=None):
    input_values={}
    output_values=defaultdict(dict)
    pull_from=these_stats()
    for n,file in enumerate(argv['input_fp']):
        #print(n+1,file)
        input_values.update(read_stats(file,pull_from))
    #print(input_values)
    csv_header=['sample_id']
    with open('vcf_stats/stats_filter.txt','w') as filter_file:
        for label in ['#SNVs','#INDELs','Ti:Tv','Het:Hom','Ins:Del']:
            values=[input_values[k][label] for k in input_values.keys()]
            mdn=numpy.median(values)
            MAD=stats.median_abs_deviation(values)
            csv_header.extend([f'{label}',f'mdn_{label}',f'MAD_{label}'])
            for k in input_values.keys():
                output_values[k].update({'sample_id':k,f'{label}':input_values[k][label],f'mdn_{label}':mdn,f'MAD_{label}':MAD})
                #print(f'{input_values[k][label]} {mdn+(4*MAD)}')
                if input_values[k][label]>mdn+(4*MAD):
                    filter_file.write(f'{k} {label}\n')
            plot_values(label,values,mdn,MAD)

    with open('vcf_stats/qc_stats.json','w') as outfile:
        json.dump(dict(output_values),outfile,indent=4)

    with open('vcf_stats/qc_stats.csv','w') as outfile:
        writer=csv.DictWriter(outfile,delimiter=',',fieldnames=csv_header)
        writer.writeheader()
        for k,v in output_values.items():
            writer.writerow(v)


if __name__=='__main__':
    main(get_args())
