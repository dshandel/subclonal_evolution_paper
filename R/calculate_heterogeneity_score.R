
# ---------------------------------------------------------------------------------------------------------------------
#                                           calculate_heterogeneity_score
# ---------------------------------------------------------------------------------------------------------------------
# Organization:         ERIBA (CRIPSR/iPSC facility)
# Programmer:           René Wardenaar
# Starting date:        06-12-18
# Last modified:        07-12-18
# Version:              1.0
# ---------------------------------------------------------------------------------------------------------------------
# 
# ----------------------
# INPUT
# ----------------------
# <cn_matrix> :         A dataframe containing all the copy number calls of different cells / libraries.
#                       Columns: 01) seqnames
#                                02) start
#                                03) end
#                                04) width
#                                05) strand
#                                06) [copy.number cell 1]
#                                07) [copy.number cell 2], etc.
# <min_num> :           Minimum number of values for calculating heterogeneity. For some cells NAs are allowed but the
#                       number of values for a given bin should be at least this high.
# ----------------------
# PROCEDURE
# ----------------------
# [01] :                Select bins for which enough cells have copy numbers (not NA).
# [02] :                Determine for each number of cells with values the total number of pairwise combinations.
# [03] :                Calculate heterogeneity scores for each bin in three different ways:
#                         01) Implementation in AneuFinder.
#                         02) New method 1: Proportion of pairwise combinations for which both cells have a different
#                                           state. Not looking at the the size of difference itself.
#                         03) New method 2: Same as 2) but this time also looking of the size of the difference.
#                                           (1-somy and 3-somy: difference of 2) 
# [04] :                Calculate weighted average using the binsize as weight.
# [05] :                Calculate weighted average for each chromosome.
# ----------------------
# OUTPUT
# ----------------------
# <results> :           A list containing all the results.
#   <genomewide> :      A dataframe containing the genome-wide heterogeneity scores.
#                       Columns: 01) Heterogeneity_0     [Implementation AneuFinder]
#                                02) Heterogeneity_1     [New method 1]
#                                03) Heterogeneity_2     [New method 2]
#  <per.chromosome> :   A dataframe containing the heterogeneity scores for each chromosome separate.
#                       Columns: The same as <genomewide>
# ---------------------------------------------------------------------------------------------------------------------



calculate_heterogeneity_score <- function(cn_matrix,min_num){
  cat("Start | calculate_heterogeneity_score","\n")
  chr_width_het         <- cn_matrix[,c(1,4)]
  chr_width_het[,1]     <- as.character(chr_width_het[,1])
  chr_width_het$num.val <- apply(cn_matrix[,6:ncol(cn_matrix)],MARGIN=1,function(x){length(which(!is.na(x)))})
  rows_keep             <- which(chr_width_het$num.val >= min_num)
  chr_width_het         <- chr_width_het[rows_keep,]
  cn_matrix             <- cn_matrix[rows_keep,]
  pair_num              <- 0
  for(num in 1:(ncol(cn_matrix)-6)){
    pair_num            <- c(pair_num,sum(1:num))
  }
  chr_width_het$num.pair <- pair_num[chr_width_het$num.val]
  cn_tables             <- apply(cn_matrix[,6:ncol(cn_matrix)],MARGIN=1,function(x){sort(table(x), decreasing = TRUE)})
  het_score_0           <- unlist(lapply(cn_tables, function(x) {sum(x * 0:(length(x) - 1))}))/chr_width_het$num.val
  het_score_1           <- NULL
  for(i in 1:length(cn_tables)){
    if(length(cn_tables[[i]]) > 1){
      new_val           <- sum(cn_tables[[i]][1:(length(cn_tables[[i]])-1)]*
                                 rev(cumsum(rev(cn_tables[[i]][2:length(cn_tables[[i]])]))))/chr_width_het$num.pair[i]
      het_score_1       <- c(het_score_1,new_val)
    }else{
      het_score_1       <- c(het_score_1,0)
    }
  }
  het_score_2           <- NULL
  for(i in 1:length(cn_tables)){
    if(length(cn_tables[[i]]) > 1){
      diff_sum          <- 0
      cn_num            <- as.numeric(names(cn_tables[[i]]))
      for(j in 1:(length(cn_tables[[i]])-1)){
        diff_abs        <- abs(cn_num[(j+1):length(cn_num)] - cn_num[j])
        new_val         <- (sum(cn_tables[[i]][(j+1):length(cn_tables[[i]])]*diff_abs)*cn_tables[[i]][j])
        diff_sum        <- diff_sum + new_val
      }
      diff_ave          <- diff_sum/chr_width_het$num.pair[i]
      het_score_2       <- c(het_score_2,diff_ave)
    }else{
      het_score_2       <- c(het_score_2,0)
    }
  }
  chr_width_het$het_0   <- het_score_0
  chr_width_het$het_1   <- het_score_1
  chr_width_het$het_2   <- het_score_2
  gw_het_0              <- stats::weighted.mean(chr_width_het$het_0,chr_width_het$width)
  gw_het_1              <- stats::weighted.mean(chr_width_het$het_1,chr_width_het$width)
  gw_het_2              <- stats::weighted.mean(chr_width_het$het_2,chr_width_het$width)
  result                <- list()
  result[["genomewide"]] <- data.frame(Heterogeneity_0=gw_het_0,Heterogeneity_1=gw_het_1,Heterogeneity_2=gw_het_2)
  chr_het_0             <- NULL                                                                                        # N: Per chromosome
  chr_het_1             <- NULL
  chr_het_2             <- NULL
  for(i in unique(chr_width_het$seqnames)){
    chr_rows            <- which(chr_width_het$seqnames == i)
    chr_het_0           <- c(chr_het_0,stats::weighted.mean(chr_width_het$het_0[chr_rows],
                                                            chr_width_het$width[chr_rows]))
    chr_het_1           <- c(chr_het_1,stats::weighted.mean(chr_width_het$het_1[chr_rows],
                                                            chr_width_het$width[chr_rows]))
    chr_het_2           <- c(chr_het_2,stats::weighted.mean(chr_width_het$het_2[chr_rows],
                                                            chr_width_het$width[chr_rows]))
  }
  result[["per.chromosome"]] <- data.frame(Heterogeneity_0=chr_het_0,Heterogeneity_1=chr_het_1,
                                           Heterogeneity_2=chr_het_2,row.names=unique(chr_width_het$seqnames))
  cat("End   | calculate_heterogeneity_score","\n")
  return(result)
}