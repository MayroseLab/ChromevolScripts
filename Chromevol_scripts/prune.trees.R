library(ape)

# args: working_dir, trees_file, remove_species_list

### PARSE INPUT ARGUMENTS
args <- commandArgs( TRUE )
for( i in 1:length(args) ){
  eval( parse( text = args[[i]] ) )
}

setwd(working_dir)
trees=try(read.tree(trees_file,keep.multi=TRUE),silent = TRUE)
if (inherits(trees,"try-error")){
	trees = read.nexus(trees_file,force.multi=TRUE)
}
drop.sp = read.table(remove_species_list)

if (length(drop.sp$V1)>0){
  for (i in 1:length(trees)){
    trees[[i]] = drop.tip(trees[[i]],as.character(drop.sp$V1))
  }
}

write.tree(trees,outfile)
