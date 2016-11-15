## Authour: Edwardo Reynolds
## 8th November 2016

## Overcoming incorrect phasing orientations 
## Read in phased GT.FORMAT file (genotypes from a vcf file)
## Read in a family file
## Assign a phasing orientation to each offspring

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1){
  stop("1 argument required:  chr#")
}
chr = as.numeric(args[1])

directory = ""
setwd(directory)
gt_filename = "...GT.FORMAT"
gt = read.table(gt_filename, header = TRUE, stringsAsFactors = FALSE)
fam_filename = "...fam"
fam = read.table(fam_filename)

## DETERMINE WHICH FAMLIES ARE PHASED MOTHER | FATHER AND THOSE PHASED FATHER | MOTHER
mother_father_ind = vector()
father_mother_ind = vector()
unknown_ind = vector()

for(i in 1:dim(fam)[1]){
  id = paste0("X",fam[i,2])
  p1 = paste0("X",fam[i,3])
  p2 = paste0("X",fam[i,4])
  family = cbind(gt[,c(1,2)],gt[,c(id,p1,p2)])
  ## find all the rows the offspring is 0|1 (there must be at least 1 in the chromosome)
  pos01 = which(family[,3] == "0|1")
  pos10 = which(family[,3] == "1|0")
  if(length(pos01) == 0  & length(pos10) == 0){
    unknown_ind = c(unknown_ind, i)
    next
  }
  allocation  = FALSE
  for(j in pos01){
    row = family[j,]
    if(row[4] == "0|0" & grepl(row[5], '1')){
      father_mother_ind = c(father_mother_ind, i)
      allocation = TRUE
      break
    }else if (row[4] == "1|1" & grepl(row[5],'0')){
      mother_father_ind = c(mother_father_ind, i)
      allocation = TRUE
      break
    } else if (row[5] == "0|0" & grepl(row[4], '1')){
      mother_father_ind = c(mother_father_ind, i)
      allocation = TRUE
      break
    } else if (row[5] == "1|1" & grepl(row[4], '0')){
      father_mother_ind = c(father_mother_ind, i)
      allocation = TRUE
      break
    }
  }
  if(!allocation){
    print("There were no 0|1's that could be decided")
    for(j in pos10){
      row = family[j,]
      if(row[4] == "0|0" & grepl(row[5], '1')){
        mother_father_ind = c(mother_father_ind, i)
        allocation = TRUE
        break
      }else if (row[4] == "1|1" & grepl(row[5],'0')){
        father_mother_ind = c(father_mother_ind, i)
        allocation = TRUE
        break
      } else if (row[5] == "0|0" & grepl(row[4], '1')){
        father_mother_ind = c(father_mother_ind, i)
        allocation = TRUE
        break
      } else if (row[5] == "1|1" & grepl(row[4], '0')){
        mother_father_ind = c(mother_father_ind, i)
        allocation = TRUE
        break
      }
    }
  }
}

### check all indivs have been allocated
if (length(mother_father_ind) + length(father_mother_ind) != dim(fam)[1]){
  print("Warning:  Not all individuals have been allocated")
  print(unknown_ind)
  print("Individuals are Homozygous")
  father_mother_ind = c(father_mother_ind, unknown_ind)
}

## split fam into two files, and then split gt into two files of individuals
mf_fam = fam[mother_father_ind,]
fm_fam = fam[father_mother_ind,]

mf_id = paste0("X", mf_fam$V2)
fm_id = paste0("X", fm_fam$V2)

mf_gt = cbind(gt[,c(1,2)], gt[,c(mf_id)])
fm_gt = cbind(gt[,c(1,2)], gt[,c(fm_id)])

## CONVERT TO A D I FORMATS
## create regular gt_offspring for A and D
gt_id = paste0("X", fam$V2)
gt_offspring = cbind(gt[,c(1,2)], gt[,c(gt_id)])

## CREATE A  ( DOES NOT NEED DIFFERENT MF/FM DATAFRAMES )
A = gt_offspring
A[A == "0|0"] = -1
A[A == "0|1"] = 0
A[A == "1|0"] = 0
A[A == "1|1"] = 1
## CREATE D
D = gt_offspring
D[D == "0|0"] = 0
D[D == "0|1"] = 1
D[D == "1|0"] = 1
D[D == "1|1"] = 0

### CREATE I
# MOTHER | FATHER
I_mf = mf_gt
I_mf[I_mf == "0|0"] = 0
I_mf[I_mf == "0|1"] = 1
I_mf[I_mf == "1|0"] = -1
I_mf[I_mf == "1|1"] = 0

# FATHER | MOTHER
I_fm = fm_gt
I_fm[I_fm == "0|0"] = 0
I_fm[I_fm == "0|1"] = -1
I_fm[I_fm == "1|0"] = 1
I_fm[I_fm == "1|1"] = 0

## Concatenate I matrices
I = cbind(I_mf, I_fm[,-c(1,2)])

## need to transpose A, D, I into indiv x marker instead of marker x indiv
## create column of IIDs (colnames)
# I
I[,-1] = apply(I[,-1], 2, as.numeric)
I = I[,-1]
I_t_matrix = t(I)
I_frame = as.data.frame(I_t_matrix)
colnames(I_frame) = as.vector(I_frame[1,])
I_frame = I_frame[-1,]
I_frame = cbind(substr(rownames(I_frame), 2, 9), I_frame)
colnames(I_frame)[1] = "IID"

# D
D[,-1] = apply(D[,-1], 2, as.numeric)
D = D[,-1]
D_t_matrix = t(D)
D_frame = as.data.frame(D_t_matrix)
colnames(D_frame) = as.vector(D_frame[1,])
D_frame = D_frame[-1,]
D_frame = cbind(substr(rownames(D_frame),2,9), D_frame)
colnames(D_frame)[1] = "IID"

# A
A[,-1] = apply(A[,-1], 2, as.numeric)
A = A[,-1]
A_t_matrix = t(A)
A_frame = as.data.frame(A_t_matrix)
colnames(A_frame) = as.vector(A_frame[1,])
A_frame = A_frame[-1,]
A_frame = cbind(substr(rownames(A_frame), 2, 9), A_frame)
colnames(A_frame)[1] = "IID"

### WRITE A, D, and I to file
info_filename = paste0("data_Beagle_chr", chr, ".txt")
write.table(A_frame, paste0("A", info_filename), quote = FALSE, row.names = FALSE)
write.table(D_frame, paste0("D", info_filename), quote = FALSE, row.names = FALSE)
write.table(I_frame, paste0("I", info_filename), quote = FALSE, row.names = FALSE)


## Add a check to make sure each chromosome has all offspring assigned.
if(dim(I_frame)[1] != dim(fam)[1]){
  print("Incorrect number of individuals assigned")
}





