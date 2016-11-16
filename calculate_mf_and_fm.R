# Authour: Edwardo Reynolds
# 8th November 2016

## Read in GT.FORMAT file (from vcf)
## Calculate the number of each heterozygote phasing orientiation (FATHER|MOTHER, or MOTHER|FATHER)

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2){
  stop("2 argument required:  chr#, filename")
}
chr = as.numeric(args[1])
gt_filename = as.character(args[2])

directory = ""
setwd(directory)
gt = read.table(gt_filename, header = TRUE, stringsAsFactors = FALSE)
### Read in family file with each individuals parent IDs
fam_filename = ""
fam = read.table(fam_filename)

### ALLOCATES EACH TRIVIAL HETEROZYGOTE OF EACH OFFSPRING TRIO TO MOTHER|FATHER OR FATHER|MOTHER
# Make sure to adjust these
indiv_start = 4000
indiv_end = 5000
# adding breaks after allocation = TRUE lines allows indivs to be put in a single vector so they can be separated
mother_father_ind = vector()
father_mother_ind = vector()
unknown_ind = vector()

for(i in indiv_start:indiv_end){
  id = paste0("X",fam[i,2])
  p1 = paste0("X",fam[i,3])
  p2 = paste0("X",fam[i,4])
  family = cbind(gt[,c(1,2)],gt[,c(id,p1,p2)])
  ## find all the rows the offspring is 0|1 
  pos01 = which(family[,3] == "0|1")
  pos10 = which(family[,3] == "1|0")
  if(length(pos01) == 0  & length(pos10) == 0){
    unknown_ind = c(unknown_ind, i)
    next
  }
  for(j in pos01){
    row = family[j,]
    if(row[4] == "0|0" & grepl(row[5], '1')){
      father_mother_ind = c(father_mother_ind, i)
    }else if (row[4] == "1|1" & grepl(row[5],'0')){
      mother_father_ind = c(mother_father_ind, i)
    } else if (row[5] == "0|0" & grepl(row[4], '1')){
      mother_father_ind = c(mother_father_ind, i)
    } else if (row[5] == "1|1" & grepl(row[4], '0')){
      father_mother_ind = c(father_mother_ind, i)
    }
  }
  for(j in pos10){
    row = family[j,]
    if(row[4] == "0|0" & grepl(row[5], '1')){
      mother_father_ind = c(mother_father_ind, i)
    }else if (row[4] == "1|1" & grepl(row[5],'0')){
      father_mother_ind = c(father_mother_ind, i)
    } else if (row[5] == "0|0" & grepl(row[4], '1')){
      father_mother_ind = c(father_mother_ind, i)
    } else if (row[5] == "1|1" & grepl(row[4], '0')){
      mother_father_ind = c(mother_father_ind, i)
    }
  }
}

### check all indivs have been allocated
## Adds unallocated individuals to FATHER|MOTHER since they are completely homozygote it doesn't matter which (for my analysis)
if (length(mother_father_ind) + length(father_mother_ind) != dim(fam)[1]){
  print("Warning:  Not all individuals have been allocated")
  father_mother_ind = c(father_mother_ind, unknown_ind)
}
## by tabling up the indices we can determine those indivs solely father|mother or mother|father,
## and those indivs  which are a mix (and to what degree)
fm = table(father_mother_ind)
mf = table(mother_father_ind)

a = which(names(fm) %in% names(mf))
b = which(names(mf) %in% names(fm))

## PLOTS indivs who are a mix
par(mai= c(1, 1, 0.5, 0.5), cex.main = 1.8, cex.lab = 1.6)
plot(as.vector(fm[a]), as.vector(mf[b]), xlab = "Sire | Dam  frequency", ylab = "Dam | Sire  frequency", main = "Frequencies of phasing orientation - Post-phased dataset", pch = 19)
text(as.vector(fm[a]), as.vector(mf[b]), labels = names(fm[a]))
legend(300, 700, "225/1000 individuals have at least one heterozygote \n of each phasing orientation across 2306 markers", cex = 1.5, bty = "n")



     