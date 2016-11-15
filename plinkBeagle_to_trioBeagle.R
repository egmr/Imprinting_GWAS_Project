## Authour: Edwardo Reynolds
## 21st October 2016
## Reads is family files (trios, pairs, nopars)
## Reads in unphased nopar Beagle 3.3.2 files (made by PLINK)
## Rearranges Beagle 3 ready files into trio, pair, and unphased (nopar) formats
## NOTE:  A pair indicates a father-offspring relationship, 
## and code will need adjustment if one wants mother-offspring relationships to also be formatted

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2){
  stop("2 arguments required:  chrStart#, chrrEnd#")
}
chrStart = as.numeric(args[1])
chrEnd = as.numeric(args[2])

directory = ""
setwd(directory)
#father-mother-offspring trios
trios_filename = ""
trios = read.table(trios_filename, header = TRUE, skipNul = TRUE)
#father-offspring pairs
pairs_filename = ""
pairs =read.table(pairs_filename, header = TRUE, skipNul = TRUE)
#no parents
nopar_filename = ""
nopar =read.table(nopar_filename, header = TRUE, skipNul = TRUE)


#Create trio data_frames depending on what pedigrees we have.
#get's file of father, mother, child trios (unphased)
#trios.final = dataframe of indivs with both parents genotypes available.
#these trios are ordered as father(p1), mother(p2), indiv. 

get_trio_file = function(fam, ped.data){
  ## create zipped order
  idx = order(c(seq_along(fam$P1), seq_along(fam$P2), seq_along(fam$IID)))
  x = unlist(c(fam$P1,fam$P2,fam$IID))[idx] # merges  p1,p2, iid like a zip
  # need to double them
  x = rep(x, each = 2)
  ## add .x, .y to them
  x = c("I.x", "IID.y", paste0(x, c(".x",".y")))
  ## Create column names of ped.data
  colnames(ped.data) = paste0(ped.data[2,], c(".x", ".y"))
  colnames(ped.data)[1:2] = c("I.x","IID.y")
  # get the right order
  trios = ped.data[,x]
  return(trios)
  
}

get_pair_file = function(fam, ped.data){
  idx = order(c(seq_along(fam$P1),  seq_along(fam$IID)))
  x = unlist(c(fam$P1,fam$IID))[idx]
  x = rep(x, each = 2)
  x = c("I.x", "IID.y", paste0(x, c(".x",".y")))
  ## Create column names of ped.data
  colnames(ped.data) = paste0(ped.data[2,], c(".x", ".y"))
  colnames(ped.data)[1:2] = c("I.x","IID.y")
  # get the right order
  pairs = ped.data[,x]
  return(pairs)
}
get_nopar_file = function(fam, ped.data){
  idx = order(seq_along(fam$IID))
  x = unlist(fam$IID)[idx]
  x = rep(x, each = 2)
  x = c("I.x", "IID.y", paste0(x, c(".x",".y")))
  ## Create column names of ped.data
  colnames(ped.data) = paste0(ped.data[2,], c(".x", ".y"))
  colnames(ped.data)[1:2] = c("I.x","IID.y")
  # get the right order
  nopar = ped.data[,x]
  return(nopar)
}

####  For each chromosome, read in the chromosomes .dat file
####  Generate 3 data.frames (trios, par_1, no_parent) for each chromosome
####  Write those data.frames to file.
directory = ""
base_filename = ""
output_filename = ""
setwd(directory)
for(i in chrStart:chrEnd){
  chr.filename = paste0(base_filename, i, ".dat") 
  chr.data = read.table(chr.filename, skipNul = TRUE, header = FALSE, stringsAsFactors = FALSE)
  #Create separated ped files.
  chr.trios = get_trio_file(trios, chr.data)
  chr.pairs = get_pair_file(pairs, chr.data)
  chr.nopar = get_nopar_file(nopar, chr.data)
  #Write these datasets to file for beagle
  write.table(chr.trios, paste0(output_filename, "_chr", i, ".trio"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(chr.pairs, paste0(output_filename, "_chr", i, ".pair"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(chr.nopar, paste0(output_filename, "_chr", i, ".nopar"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  ##Not sure if necessary, since they will be overwritten quickly anyway, but removing them anyway.
  rm(chr.data, chr.trios, chr.pairs, chr.nopar, chr.filename)
}
