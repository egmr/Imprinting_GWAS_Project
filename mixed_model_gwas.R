## Authour: Edwardo Reynolds
##  6 September 2016

## Mixed models using Additive, Dominance, and Imprinting effects
## Imports SNP effects, Fixed effects, Random Effects 
## Runs mixed model analysis on markers as part of GWAS, using lme4

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 3){
  stop("3 arguments required:  chr#, markerStart#, markerEnd#")
}
chr = args[1]
markerStart = as.numeric(args[2])
markerEnd = as.numeric(args[3])

#### READ IN A_DATA, D_DATA_I_DATA FOR THE CHROMOSOME
A_data = read.table(paste0("A_data.chr", chr), header = TRUE)
D_data = read.table(paste0("D_data.chr", chr), header = TRUE)
I_data = read.table(paste0("I_data.chr", chr), header = TRUE)

### Phenotypes
pheno_filename = ""
phenos = read.table(pheno_filename, header=TRUE)
## get the phenotype of interest, and the IID, removing rows of missing phenoypes

## Get other fixed effects
fixed_filename = ""
cohorts = read.table(fixed_filename, header = TRUE)


#### MIXED MODELS
## Adapted from: Timothée Flutre, https://gist.github.com/timflutre/43daacf2c8868f609489

suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(lme4))

#Data that is constant per marker
marker = 1000
data_pheno_cohort = merge(phenos, cohorts, by="IID")
data_p_c_a = merge(data_pheno_cohort, A_data[,c(1, marker)], by = "IID")  
data_p_c_a_d = merge(data_p_c_a, D_data[,c(1,marker)], by = "IID") 
data_p_c_a_d_i = merge(data_p_c_a_d, I_data[,c(1,marker)], by="IID") 
colnames(data_p_c_a_d_i)[4:6] = c("A","D","I")

#Formulas
formulaA = Pheno ~ 1 + cohort + A + (1|IID)
formulaAD = Pheno ~ 1 + cohort + A + D + (1|IID)
formulaADI = Pheno ~ 1 + cohort + A + D + I + (1|IID)

## Read in relationship matrix (variance component of random effect term)
rel_mat  = read.table("....rel")
relmat_ids = read.table("...rel.id")

## Be careful to ensure the IIDs of the relationship matrix are in the same order as those fixed effects

rel_mat = as.matrix(rel_mat)
rel_mat = forceSymmetric(rel_mat)
relmat <- list(IID = rel_mat)
relfac <- relmat

## Setting up 
parsedFormula <- lFormula(formula=formulaADI, data=data_p_c_a_d_i, REML = FALSE, control=lmerControl(check.nobs.vs.nlev="ignore", check.nobs.vs.nRE="ignore"))
flist <- parsedFormula$reTrms[["flist"]] # list of grouping factors
Ztlist <- parsedFormula$reTrms[["Ztlist"]] # list of transpose of the sparse model matrices
relmat[[1]] <- Matrix(relmat[[1]], sparse=TRUE)
relfac[[1]] <- chol(relmat[[1]])
Ztlist[[1]] <- relfac[[1]] %*% Ztlist[[1]]


#The Following are required for each marker in loop
# Run the models on each marker
markerNumber = markerEnd - markerStart + 1
A_lik = rep(0, markerNumber)
D_lik = rep(0, markerNumber)
I_lik = rep(0, markerNumber)
AA_ttest = rep(0, markerNumber)
D_ttest = data.frame(rep(0,markerNumber), rep(0, markerNumber))
I_ttest = data.frame(rep(0,markerNumber), rep(0, markerNumber), rep(0, markerNumber))
colnames(D_ttest) = c("DA_ttest", "DD_ttest")
colnames(I_ttest) = c("IA_ttest", "ID_ttest", "II_ttest")

for(i in markerStart:markerEnd){
  data_p_c_a = merge(data_pheno_cohort, A_data[,c(1, i)], by = "IID")  
  data_p_c_a_d = merge(data_p_c_a, D_data[,c(1, i)], by = "IID") 
  data_p_c_a_d_i = merge(data_p_c_a_d, I_data[,c(1, i)], by="IID") 
  colnames(data_p_c_a_d_i)[4:6] = c("A","D","I")
  rm(data_p_c_a, data_p_c_a_d)

  ##parse formulas for each formula
  parsedFormulaA <- lFormula(formula=formulaA, data=data_p_c_a_d_i, REML = FALSE,
                             control=lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.nRE="ignore"))  
  parsedFormulaAD <- lFormula(formula=formulaAD, data=data_p_c_a_d_i, REML = FALSE,
                              control=lmerControl(check.nobs.vs.nlev="ignore",
                                                  check.nobs.vs.nRE="ignore"))
  parsedFormulaADI <- lFormula(formula=formulaADI, data=data_p_c_a_d_i, REML = FALSE,
                               control=lmerControl(check.nobs.vs.nlev="ignore",
                                                   check.nobs.vs.nRE="ignore"))
  
  parsedFormulaA$reTrms[["Ztlist"]] <- Ztlist
  parsedFormulaA$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
  parsedFormulaAD$reTrms[["Ztlist"]] <- Ztlist
  parsedFormulaAD$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
  parsedFormulaADI$reTrms[["Ztlist"]] <- Ztlist
  parsedFormulaADI$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
  
  devianceFunctionA <- do.call(mkLmerDevfun, parsedFormulaA)
  devianceFunctionAD <- do.call(mkLmerDevfun, parsedFormulaAD)
  devianceFunctionADI <- do.call(mkLmerDevfun, parsedFormulaADI)
  
  
  optimizerOutputA <- optimizeLmer(devianceFunctionA)
  optimizerOutputAD <- optimizeLmer(devianceFunctionAD)
  optimizerOutputADI <- optimizeLmer(devianceFunctionADI)
  
  fit.A <- mkMerMod(rho=environment(devianceFunctionA), opt=optimizerOutputA, reTrms=parsedFormulaA$reTrms, fr=parsedFormulaA$fr)
  fit.AD <- mkMerMod(rho=environment(devianceFunctionAD), opt=optimizerOutputAD, reTrms=parsedFormulaAD$reTrms, fr=parsedFormulaAD$fr)
  fit.ADI <- mkMerMod(rho=environment(devianceFunctionADI), opt=optimizerOutputADI, reTrms=parsedFormulaADI$reTrms, fr=parsedFormulaADI$fr)
  
  frame_position = i - markerStart + 1
  A_lik[frame_position] = summary(fit.A)$logLik
  D_lik[frame_position] = summary(fit.AD)$logLik
  I_lik[frame_position] = summary(fit.ADI)$logLik
  AA_ttest[frame_position] = data.frame(summary(fit.A)$coefficients)[3,3]
  D_ttest[frame_position,] = data.frame(summary(fit.AD)$coefficients)[3:4,3]
  I_ttest[frame_position,] = data.frame(summary(fit.ADI)$coefficients)[3:5,3]
  
}


output_table = cbind(A_lik, D_lik, I_lik, AA_ttest, D_ttest, I_ttest)
output_name = paste0("Chr", chr, "_Marker", markerStart, "_", markerEnd, "_mixedModeloutput.txt")
write.table(output_table, output_name, quote = FALSE, row.names = FALSE)


