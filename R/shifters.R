#############
# Important #
#############
# These scripts were by René Wardenaar and kindly provided by the Kops group, Huybrecht, Utrecht.

# ---------------------------------------------------------------------------------------------------------------------
#                                                  create_cn_matrix
# ---------------------------------------------------------------------------------------------------------------------
# Organization:         ERIBA (CRIPSR/iPSC facility)
# Programmer:           René Wardenaar
# Starting date:        06-12-18
# Last modified:        17-12-18
# Version:              1.0
# ---------------------------------------------------------------------------------------------------------------------
# 
# ----------------------
# INPUT
# ----------------------
# <file_names> :        A vector with file names. The files (.RData) should contain model results for which the copy
#                       number calls are goining to be combined. Note: all files should contain copy number calls.
# ----------------------
# PROCEDURE
# ----------------------
# [01] :                Read first file and select chromosome name, the start and end position of the bins, the width
#                       of the bins, the strand of the bins and the copy number of the bins.
# [02] :                Load other files, select the same columns and merge with first.
# ----------------------
# OUTPUT
# ----------------------
# <model_all_calls> :   A dataframe containing all the copy number calls of all the files.
#                       Columns: 01) seqnames
#                                02) start
#                                03) end
#                                04) width
#                                05) strand
#                                06) [copy.number cell 1]
#                                07) [copy.number cell 2], etc.
# ---------------------------------------------------------------------------------------------------------------------

create_cn_matrix <- function(file_names){
  # relocate cells if at home

  cat("Start | create_cn_matrix","\n")
  cat("      | number of files: ",length(file_names),"\n",sep="")
  load(file_names[1])                                                                                                  # N: Object is called 'model'
  model_all_calls       <- as.data.frame(model$bins)
  dfnames               <- names(model_all_calls)
  if(any(dfnames == "copy.number")){
    type                <- "wcn"
    model_all_calls     <- model_all_calls[,c("seqnames","start","end","width","strand","copy.number")]
  }else if(any(dfnames == "state")){
    type                <- "ncn"
    model_all_calls     <- model_all_calls[,c("seqnames","start","end","width","strand","state")]
    model_all_calls$state <- as.character(model_all_calls$state)
    model_all_calls$state <- substr(model_all_calls$state,1,nchar(model_all_calls$state)-5)
    model_all_calls$state <- as.numeric(model_all_calls$state)
  }else{
    stop("There is no 'copy number' or 'state' column!")
  }
  model_all_calls$seqnames <- as.character(model_all_calls$seqnames)
  model_all_calls$strand <- as.character(model_all_calls$strand)
  names(model_all_calls)[ncol(model_all_calls)] <- basename(file_names[1])
  for(i in 2:length(file_names)){
    load(file_names[i])                                                                                                # N: Object is called 'model'
    model_cur_calls     <- as.data.frame(model$bins)
    dfnames             <- names(model_cur_calls)
    if(any(dfnames %in% c("copy.number","state"))){
      if(type == "wcn"){
        model_cur_calls <- model_cur_calls[,c("seqnames","start","end","width","strand","copy.number")]
      }else if(type == "ncn"){
        model_cur_calls <- model_cur_calls[,c("seqnames","start","end","width","strand","state")]
        model_cur_calls$state <- as.character(model_cur_calls$state)
        model_cur_calls$state <- substr(model_cur_calls$state,1,nchar(model_cur_calls$state)-5)
        model_cur_calls$state <- as.numeric(model_cur_calls$state)
      }
      model_cur_calls$seqnames <- as.character(model_cur_calls$seqnames)
      model_cur_calls$strand <- as.character(model_cur_calls$strand)
      names(model_cur_calls)[ncol(model_cur_calls)] <- basename(file_names[i])
      model_all_calls   <- merge(model_all_calls,model_cur_calls)
    }else{
      stop("There is no 'copy number' or 'state' column!")
    }
  }
  model_all_calls       <- model_all_calls[order(model_all_calls$seqnames,model_all_calls$start,model_all_calls$end),]
  cat("End   | create_cn_matrix","\n")
  return(model_all_calls)
}

# ---------------------------------------------------------------------------------------------------------------------
#                                            cluster_and_shift_transitions
# ---------------------------------------------------------------------------------------------------------------------
# Organization:         ERIBA (CRIPSR/iPSC facility)
# Programmer:           René Wardenaar
# Starting date:        06-12-18
# Last modified:        06-12-18
# Version:              1.0
# ---------------------------------------------------------------------------------------------------------------------
# 
# ----------------------
# INPUT
# ----------------------
# <model_all_calls> :   A dataframe containing all the copy number calls of different cells / libraries.
#                       Columns: 01) seqnames
#                                02) start
#                                03) end
#                                04) width
#                                05) strand
#                                06) [copy.number cell 1]
#                                07) [copy.number cell 2], etc.
#
# <max_dist> :          The maximum distance for clustering copy number state transitions between and within cells.
# ----------------------
# PROCEDURE
# ----------------------
# [01] :                Clustering of any transition within a maximum distance (max_dist)
# [02] :                When a cell has more than one transition within the same cluster; take the middle of the first
#                       and last transition as the transition point of this cell.
# [03] :                Determine the middle transition position of all transitions within the cluster. All transitions
#                       are shifted to this transition position.
# [04] :                The flanking states of the cluster are the new states.
# ----------------------
# OUTPUT
# ----------------------
# <results> :           A list containing all the results.
#   <matrix.shift> :    A dataframe similar to the input (<model_all_calls>) but now the copy number state transitions
#                       have been shifted to one shared position.
#   <segments> :        A dataframe similar to <matrix.shift> but now the bins have been collapsed into segments with
#                       the same copy number.
#   <cluster.info> :    A dataframe with information about the clusters.
#                       Note: It will return an empty list when there are no clusters.
#                       Columns: 01) chrom
#                                02) start       [start first bin after first transition inside cluster]
#                                03) end         [start first bin after last transition inside cluster]
#                                04) width       [difference between start and end]
#                                05) num.single  [number of cells with one transition]
#                                06) num.mult    [number of cells with more than one transition]
#                                07) cl.id       [cluster id]
#   <multiple.breakpoints> : A dataframe with information about cases when cells have multiple transitions.
#                            Note: It will return an empty list when there are no cases with multiple transitions.
#                            Columns: 01) cell   [name of the cell; file name]
#                                     02) chrom
#                                     03) start
#                                     04) end
#                                     05) width
#                                     06) cl.id  [cluster id]
# ---------------------------------------------------------------------------------------------------------------------



cluster_and_shift_transitions <- function(model_all_calls, max_dist){
  cat("Start | cluster_and_shift_transitions","\n")
  model_all_calls[,1]   <- as.character(model_all_calls[,1])
  model_all_calls       <- model_all_calls[order(model_all_calls[,1],model_all_calls[,2],model_all_calls[,3]),]
  model_all_calls_shift <- list()
  segment_info          <- list()
  mult_info             <- list()
  clust_info            <- list()
  for(chr in unique(model_all_calls$seqnames)){
    mat_chr             <- model_all_calls[which(model_all_calls$seqnames == chr),]
    seg_chr             <- mat_chr[1,]
    transpos            <- mat_chr$start[2:nrow(mat_chr)]
    state_diff          <- as.matrix(abs(mat_chr[1:(nrow(mat_chr)-1),6:ncol(mat_chr)] -
                                           mat_chr[2:nrow(mat_chr),6:ncol(mat_chr)]))
    pos_diff            <- which(state_diff > 0)
    if(length(pos_diff) > 0){
      state_diff[pos_diff] <- 1
      row_sum           <- apply(state_diff,MARGIN=1,sum)
      trans_rows        <- which(row_sum > 0)
      cluster           <- list()
      cluster[[1]]      <- c(transpos[trans_rows[1]],transpos[trans_rows[1]])
      if(length(trans_rows) > 1){
        for(i in 2:length(trans_rows)){
          if(transpos[trans_rows[i]] <= (cluster[[length(cluster)]][2]+max_dist)){                                     # N: Within cluster range
            cluster[[length(cluster)]][2] <- transpos[trans_rows[i]]
          }else{                                                                                                       # N: Start new cluster
            cluster[[(length(cluster)+1)]] <- c(transpos[trans_rows[i]],transpos[trans_rows[i]])
          }
        }
      }
      for(i in 1:length(cluster)){
        ri_pre          <- nrow(seg_chr) 
        ri_new          <- nrow(seg_chr) + 1
        start_ind       <- which(mat_chr$start == cluster[[i]][1]) - 1                                                 # N: Start index
        end_ind         <- which(mat_chr$start == cluster[[i]][2]) - 1                                                 # N: End index -> Minus one to get index of difference matrix (state_diff)
        num_trans       <- colSums(state_diff[start_ind:end_ind,,drop=FALSE])
        ind_one         <- which(num_trans == 1)
        ind_mult        <- which(num_trans > 1)                                                                        # N: Multiple transitions for one cell in the same cluster
        cp_one          <- NULL
        cp_two          <- NULL
        if(length(ind_one) > 0){ 
          cp_one        <- apply(state_diff[start_ind:end_ind,ind_one,drop=FALSE],MARGIN=2,which.max)                  # N: Cluster position one (location transition within cluster)
        }
        if(length(ind_mult) > 0){ 
          cp_two        <- apply(state_diff[start_ind:end_ind,ind_mult,drop=FALSE],MARGIN=2,                           # N: Cluster position mult (location middle transition within cluster)
                                 function(x){floor((min(which(x == 1)) + max(which(x == 1)))/2)})
        }
        median_ind      <- floor(median(c(cp_one,cp_two)))
        seg_chr         <- rbind(seg_chr,mat_chr[(end_ind+1),])
        seg_chr$end[ri_pre] <- mat_chr$end[(start_ind + median_ind - 1)]
        seg_chr$width[ri_pre] <- (seg_chr$end[ri_pre] - seg_chr$start[ri_pre] + 1)
        seg_chr$start[ri_new] <- mat_chr$start[(start_ind + median_ind)]
        num_before      <- median_ind
        num_after       <- (end_ind - start_ind + 2 - median_ind)
        for(j in c(ind_one,ind_mult)){
          shifted_states <- c(rep(mat_chr[start_ind,(j+5)],num_before),rep(mat_chr[(end_ind+1),(j+5)],num_after))
          mat_chr[start_ind:(end_ind+1),(j+5)] <- shifted_states
        }
        if(length(ind_mult) > 0){
          for(j in ind_mult){                                                                                          # N: Information of cases where one cell has multiple transitions in the same cluster
            mult_info[[length(mult_info)+1]] <- data.frame(cell=names(mat_chr)[(j+5)],chrom=chr,
                                                           start=mat_chr$start[(start_ind+1)],end=mat_chr$start[(end_ind+1)],   # N: Start: start bin of bin after first transition    End: start bin of bin after last transition
                                                           width=(mat_chr$start[(end_ind+1)]-mat_chr$start[(start_ind+1)]),
                                                           cl.id=paste0(chr,"_",i))
          }
        }
        clust_info[[length(clust_info)+1]] <- data.frame(chrom=chr,start=mat_chr$start[(start_ind+1)],                 # N: Start: start bin of bin after first transition    End: start bin of bin after last transition
                                                         end=mat_chr$start[(end_ind+1)],
                                                         width=(mat_chr$start[(end_ind+1)]-mat_chr$start[(start_ind+1)]),
                                                         num.single=length(ind_one),num.mult=length(ind_mult),
                                                         cl.id=paste0(chr,"_",i))
      }
    }
    model_all_calls_shift[[length(model_all_calls_shift)+1]] <- mat_chr
    seg_chr$end[nrow(seg_chr)] <- mat_chr$end[nrow(mat_chr)]
    seg_chr$width[nrow(seg_chr)] <- (seg_chr$end[nrow(seg_chr)] - seg_chr$start[nrow(seg_chr)] + 1)
    segment_info[[length(segment_info)+1]] <- seg_chr
  }
  model_all_calls_shift <- do.call("rbind",model_all_calls_shift)
  segment_info          <- do.call("rbind",segment_info)
  if(length(mult_info) > 0){
    mult_info           <- do.call("rbind",mult_info)
    mult_info           <- mult_info[order(mult_info$cell,mult_info$chrom,mult_info$start,mult_info$end),]
  }
  if(length(clust_info) > 0){
    clust_info          <- do.call("rbind",clust_info)
    clust_info          <- clust_info[order(clust_info$chrom,clust_info$start,clust_info$end),]
  }
  results               <- list()
  results$matrix.shift  <- model_all_calls_shift
  results$segments      <- segment_info
  results$cluster.info  <- clust_info
  results$multiple.breakpoints <- mult_info
  tot_cells_cl          <- sum(clust_info$num.single) + sum(clust_info$num.mult)
  per_single            <- round((sum(clust_info$num.single) / tot_cells_cl)*100,2)
  per_mult              <- round((sum(clust_info$num.mult) / tot_cells_cl)*100,2)
  cat("      | number of clusters: ",nrow(clust_info),"\n",sep="")
  cat("      | number of cells single transitions: ",sum(clust_info$num.single)," (",per_single,")","\n",sep="")
  cat("      | number of cells multiple transitions: ",sum(clust_info$num.mult)," (",per_mult,")","\n",sep="")
  cat("End   | cluster_and_shift_transitions","\n")
  return(results)
}

