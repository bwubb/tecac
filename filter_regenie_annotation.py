import argparse
import os

p=argparse.ArgumentParser(description="Drop variant IDs from regenie.annotation.txt style lines (ID GENE mask).")
p.add_argument("--annotation",required=True,help="regenie.annotation.txt")
p.add_argument("--exclude-ids",required=True,help="One variant ID per line")
p.add_argument("--output",required=True,help="Filtered annotation output")
args=p.parse_args()

with open(args.exclude_ids,'r') as f:
    ex={line.strip() for line in f if line.strip()}
os.makedirs(os.path.dirname(args.output) or '.',exist_ok=True)
kept=0
dropped=0
with open(args.annotation,'r') as inf,open(args.output,'w') as outf:
    for line in inf:
        line=line.rstrip('\n')
        if not line.strip():
            continue
        vid=line.split()[0]
        if vid in ex:
            dropped+=1
            continue
        outf.write(line+'\n')
        kept+=1
print(f"kept={kept} dropped={dropped} wrote {args.output}")
