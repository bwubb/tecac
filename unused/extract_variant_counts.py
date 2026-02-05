
import os
import argparse
import pandas as pd

def read_data(all_files):
    frames=[]
    for file in all_files:
        df=pd.read_csv(file,header=None,delim_whitespace=True)
        #not sure if csv is better output, and if bcftools query will do it.
        frames.append(df)
    df=pd.concat(frames,axis=0,ignore_index=True)
    return df

def count_filter(df,x=2):
    u=df['id'].value_counts()<x
    return df.loc[df['id'].isin(u[u.values].index)].reset_index()

def get_args(argv):
    p=argparse.ArgumentParser()
    p.add_argument('-o','--out_fp',help='Path to output file.')
    p.add_argument('-c','--count',type=int,help='Passing variants count less than this number.')
    p.add_argument('all_files',nargs=argparse.REMAINDER,help='Input files')
    return p.parse_args(argv)

def main(argv=None):
    args=get_args(argv)
    df=count_filter(read_data(args.all_files),args.count)
    df.to_csv(args.out_fp,sep=' ',header=False,index=False)

if __name__=='__main__':
    main()
