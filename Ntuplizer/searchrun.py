def getFileListDAS(dataset,instance="prod/global",run=-1):

   cmd = 'dasgoclient --query="file dataset=%s instance=%s" --limit=10000' %(dataset,instance)
   if run != -1: cmd = 'dasgoclient --query="file run=%i dataset=%s instance=%s" --limit=10000' %(run,dataset,instance)
   print cmd
   cmd_out = commands.getoutput( cmd )
   tmpList = cmd_out.split(os.linesep)
   files = []
   for l in tmpList:
      if l.find(".root") != -1:
         files.append(l)

   return files
