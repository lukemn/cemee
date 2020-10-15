args = commandArgs(trailingOnly=F)
print(args)
f=args[grep("file=", args)]
print(strsplit(f, '=')[[1]][2])


