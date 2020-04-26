# 11 April 2020 / Thinh
# clustering based on non parametric quantile regression for heterogeneous samples

library(cluster)
library(quantreg)

###############################################################################
# functions for detecting outlier based on E(X2) ~ degree of freedom
###############################################################################
get_error = function(Y, EST, n_sample){
  error = sapply(1:length(Y), function(i){
    y = Y[[i]]
    y = y[ , (ncol(y) - n_sample + 1):ncol(y)]
    est = EST[[i]]
    yhat = est$X %*% t(est$BETA)
    e = (y - yhat)/sqrt(yhat+0.1)
    return(e)
  })
  return(error)
}

get_degree = function(Y, EST, n_sample){
  degree = sapply(1:length(Y), function(i){
    est = EST[[i]]
    degree = n_sample * nrow(est$X) - n_sample*ncol(est$X) - nrow(est$X)*ncol(est$X)
    return(degree)
  })

  degree = ifelse(degree < 0, 1, degree)
  return(degree)
}

get_chi_square = function(error){
  chi_square_statistic = sapply(1:length(error), function(i){
    sum(error[[i]]^2)
  })
  return(chi_square_statistic)
}

exec_quant_reg = function(chi_square_statistic, degree, tau){
  nonpar_quanreg = rqss(chi_square_statistic~qss(degree, constraint="V"), tau=tau)$coef
  temp = sort(unique(degree))
  quan_res = data.frame(degree = temp[-1], chi_square = nonpar_quanreg[1]+nonpar_quanreg[-1])
  outlier_ids = c()
  for(i in 1:length(degree)){
    if(degree[i] > 1 & chi_square_statistic[i] > 0){
      maximum_chi_square = quan_res[which(quan_res$degree == degree[i]), ]$chi_square
      if(chi_square_statistic[i] > maximum_chi_square){
        outlier_ids = c(outlier_ids, i) 
      }
    }
  }
  return(outlier_ids)
}

###############################################################################
# functions for clustering on CRPs of outliers
###############################################################################

kmeans_func = function(y, k = 2){
  kmeans_results = kmeans(y, centers = k, nstart = 2)
  return(kmeans_results$cluster)
}

# silhouette method
# first training kmeans for different number of clusters k
# then calculating silhouette value (mean of all data points)
average_silhouette = function(y, k){
  kmeans_results = kmeans_func(y, k)
  cor_mat = suppressWarnings(cor(t(y), method="spearman"))
  diss = 1 - cor_mat
  diss[is.na(diss)] = 1
  silhouette_value = silhouette(kmeans_results, diss)
  return(mean(silhouette_value[, "sil_width"]))
}

silhouette_func = function(y, n_sample){
  max_n_cluster = round(n_sample/5)
  silhouette_values = sapply(2:max_n_cluster, function(k) average_silhouette(y, k))
  optimal_k = which(max(silhouette_values) == silhouette_values) + 1
  if(length(optimal_k) > 1){
    optimal_k = optimal_k[1]
  }
  # print(optimal_k)
  return(kmeans_func(y, optimal_k))
}

get_clusters = function(y, func, error, error_dis = FALSE, n_sample){
  # input: Ycount for each CRP
  # output: splitted Ycount for each cluster
  
  # backup y
  y_backup = list(y)
  
  # split isoform and sample from Ycount
  # using sample_equiv to cluster
  isoform_equiv = y[, 1:(ncol(y)-n_sample), drop=FALSE]
  y = y[, (ncol(y)-n_sample+1):ncol(y)]
  
  error_matrix = t(error)
  
  isoform_equiv = t(isoform_equiv)
  y = t(y)

  y_new = tryCatch(
    {
      if(error_dis){
        clusters = func(error_matrix, n_sample)
      } else {
        clusters = func(log2(y+1), n_sample)
      }
      cluster_per_samples = table(clusters)
      n_cluster = length(cluster_per_samples)
      
      #if exists a cluster has less than 3 samples, run kmeans(k=2)
      if(sum(cluster_per_samples < 3) >= 1){
        # print("Silhouette: < 3")
        clusters = kmeans_func(error_matrix)
        cluster_per_samples = table(clusters)
        n_cluster = length(cluster_per_samples)
        if(sum(cluster_per_samples < 3) >= 1){
          return(y_backup)
        }
      } 
      # if not, create new y by clusters and concat with isoform_equiv
      return(lapply(1:n_cluster, 
              function(x) t(rbind(isoform_equiv, y[which(clusters == x), , drop=FALSE]))))
      
    },
    error = function(cond){
      print(cond)
      return(y_backup)
    }
  )
  
}

###############################################################################
# main function
###############################################################################

exec_cluster = function(Y, EST, tau){
  n_sample = nrow(EST[[1]]$BETA)
  errors = get_error(Y, EST, n_sample)
  chi_square_statistic = get_chi_square(errors)
  degrees = get_degree(Y, EST, n_sample)
  outliers = exec_quant_reg(chi_square_statistic, degrees, tau)

  cluster_Y = lapply(outliers, function(i){
    y = Y[[i]]
    error = errors[[i]]
    res = get_clusters(y, silhouette_func, error, error_dis = TRUE, n_sample)
    return(res)
  })
  
  cluster_Y = unlist(cluster_Y, recursive = FALSE)

  original_Y = lapply(1:length(Y), function(i){
    if((i %in% outliers) == FALSE){
      return(Y[[i]])
    }
  })
  original_Y = original_Y[lengths(original_Y) != 0]

  original_EST = lapply(1:length(EST), function(i){
    if((i %in% outliers) == FALSE){
      return(EST[[i]])
    }
  })
  original_EST = original_EST[lengths(original_EST) != 0]

  return(list(original_Y=original_Y, cluster_Y=cluster_Y, 
              original_EST=original_EST))
}