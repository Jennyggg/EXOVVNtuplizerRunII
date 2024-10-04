
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--inputfile',required=True,help="The input file")
parser.add_argument('--outputfile',required=True,help="The output file")
args = parser.parse_args()


inputlist=[]
with open(args.inputfile,'r') as fin:
  line = fin.readline()
  while line:
    if "'file:///pnfs/psi.ch/cms/trivcat/store/user/kcormier/instantons/instantons/miniAODs/" in line:
      inputlist.append(line.replace("'","").replace(",","").strip())
    line=fin.readline()
with open(args.outputfile,'w') as fout:
  fout.write("\n".join(inputlist))


