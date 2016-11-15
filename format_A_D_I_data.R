
# Authour: Edwardo Reynolds
# 13 September 2016

### Reads in trio and pair phased genotype output from Beagle 3.3.2 
### Creates A_data, D_data, I_data for a chromosome and writes it to file,
### Writes data to file it can be used in linear model and mixed model analyses

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1){
  stop("Incorrect number of arguments.  1 required: Chromosome number")
}
chr = args[1]

### Read in phased trios and pairs files, remove excess columns for each, concatenate the trios and pairs.

###TRIOS
trio_directory = ""
trio_filename = paste0("chr",chr,".chr", chr, ".trio.phased")
setwd(trio_directory)
data.trio = read.table(trio_filename)  
##Remove rows 1,3 and col 1 as uninformative.
data.trio = data.trio[-c(1,3),-1]
##Removes parents, leaving only phased individuals
indiv_index_trios = c(T, rep(c(F, F, F, F, T, T), length = dim(data.trio)[2] - 1))
data.trio = data.trio[,indiv_index_trios]  

###PAIRS
pair_directory = ""
pair_filename = paste0("chr",chr,".chr",chr,".pair.phased")
setwd(pair_direcotry)
data.pair = read.table(pair_filename)
data.pair = data.pair[-c(1,3),-1]
##Removes parents, leaving only phased individuals
indiv_index_pair = c(T, rep(c(F, F, T, T), length = dim(data.pair)[2] - 1)) 
data.pair = data.pair[,indiv_index_pair] 

data = cbind(data.trio, data.pair[,-1])
rm(data.pair, data.trio, indiv_index_pair, indiv_index_trios)



### Set imputed values to missings (0s)
### Import missing mendel data
mendel_directory = ""
mendel_filename = "mendel_errors.mendel"
setwd(mendel_directory)
mendel = read.table(mendel_filename, header  =TRUE)
mendel_chr = mendel[mendel$CHR == as.numeric(chr),]
mendel_chr = mendel_chr[,c(2,4)]
mendel_chr$SNP = as.character(mendel_chr$SNP)
##Remove all inconsistencies in mendel file
for(i in 1:dim(mendel_chr)[1]){
  rows = which(data[,1] == mendel_chr$SNP[i])
  cols = which(data[1,] == mendel_chr$KID[i])
  data[rows, cols] = 0
}
rm(mendel, mendel_chr, cols, rows, i) 


### Translate genotypes to ADI files, create A, D, I, makes missings the averages 
data_geno = as.matrix(data[-1,-1])
#Change 1,2's to 0.5, -0.5's respectively'
data_geno[data_geno == 1] = 0.5
data_geno[data_geno == 2] = -0.5
data_geno[data_geno == 0] = NA

##Now this matrix can be used to calculate the different A, D, I.

##So create A:
## A + B
A = matrix(0,nrow(data_geno), ncol(data_geno)/2)  
i = seq(from=1, to=(ncol(data_geno)/2), by=1) ##number to new columns
j = seq(from=1, to=ncol(data_geno), by=2) ##every second column index
A[,i] = data_geno[,j]+data_geno[,j+1] 

##Create D:
## abs(A - B) , use same i, j as above
D = matrix(0,nrow(data_geno), ncol(data_geno)/2)  
D[,i] = abs(data_geno[,j]-data_geno[,j+1]) 

##Create I:
## A - B
I = matrix(0,nrow(data_geno), ncol(data_geno)/2) 
I[,i] = data_geno[,j]-data_geno[,j+1] 

rm(data_geno)

## NA correction can likely become a function
## Give NA's mean values
for (i in 1:ncol(A)){
  n = 2*nrow(A)
  A[,i][is.na(A[,i])] = ((((2*as.numeric(table((A[,i]==1))[2])) + as.numeric(table((A[,i]==0))[2]))/n) - 0.5)
}
for (i in 1:ncol(D)){
  n = 2*nrow(D)
  D[,i][is.na(D[,i])] = ((((2*as.numeric(table((D[,i]==1))[2])) + as.numeric(table((D[,i]==0))[2]))/n) - 0.5)
}

for (i in 1:ncol(I)){
  n = 2*nrow(I)
  I[,i][is.na(I[,i])] = ((((2*as.numeric(table((I[,i]==1))[2])) + as.numeric(table((I[,i]==0))[2]))/n) - 0.5)
}

rm(i, j, n)


### Transpose matrices
A = data.frame(t(A))
D = data.frame(t(D))
I = data.frame(t(I))

## Reassign IDs and marker IDs to matrices
ids_double = as.vector(data[1,-1])
every_second_id = rep(c(T,F), length = length(ids_double))
ids = ids_double[every_second_id]  ##Individual IDs
ids = t(ids)

marker_ids = as.vector(data[,1])

##bind IIDs
A_data = cbind(ids, A)
D_data = cbind(ids, D)
I_data = cbind(ids, I)
##rename cols as Markers
colnames(A_data) = marker_ids 
colnames(D_data) = marker_ids
colnames(I_data) = marker_ids

rm(data, A, D, I, marker_ids, every_second_id, ids, ids_double)

output_directory = ""
setwd(output_directory)

#### Write A_data, D_data, I_data to file
write.table(A_data, paste0("A_data.chr",chr), row.names=FALSE, quote = FALSE)
write.table(D_data, paste0("D_data.chr",chr), row.names=FALSE, quote = FALSE)
write.table(I_data, paste0("I_data.chr",chr), row.names=FALSE, quote = FALSE)
