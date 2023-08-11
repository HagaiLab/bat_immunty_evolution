############################## Creates CSV files with x axis and y axis corresponding to specie cell types #################################
# Make a loop that make all possible combinations from de_species_list, and compares all species one to each other
corrCellTypesPWT_Intersection = function(de_sp1, de_sp2, sp_pair = c("species1", "species2"),
                            ct_sets = NULL, sim_thr = 0.2,
                            gene_filter = NULL, norm_all = T, filter_ct = T){
  # which cells to use
  if(is.null(ct_sets)){
    ct_sets = list(unname(colnames(de_sp1$avg_exp)),
                   unname(colnames(de_sp2$avg_exp)))
    names(ct_sets) = sp_pair
  }
  # using marker intersection
  # Get a list of all available de_species in de_species_list
  
  # Get a list of all intersected markers for each de_species and combine them
  all_intersected_markers <- lapply(de_species_list, function(de_species) {
    unlist(de_species$intersected_markers)
  })
  combined_intersected_markers <- unique(unlist(all_intersected_markers))

  # Intersect the combined markers with each de_species and save the results in s1 and s2
  s1 <- Reduce(intersect, lapply(de_species_list, function(de_species) {
    intersected_markers <- unlist(de_species$intersected_markers)
    intersect(intersected_markers, combined_intersected_markers)
  }))

  s2 <- Reduce(intersect, lapply(de_species_list, function(de_species) {
    intersected_markers <- unlist(de_species$intersected_markers)
    intersect(intersected_markers, combined_intersected_markers)
  }))
  #s1 = intersect(unique(unlist(de_species1$intersected_markers)),c(unique(unlist(de_species2$intersected_markers)),unique(unlist(de_species3$intersected_markers))
                                                         #unique(unlist(de_species4$intersected_markers))#,unique(unlist(de_species5$intersected_markers)),
                                                         #unique(unlist(de_species6$intersected_markers)),unique(unlist(de_species7$intersected_markers))
  #))
  #s2 = intersect(unique(unlist(de_species1$intersected_markers)),c(unique(unlist(de_species2$intersected_markers)),unique(unlist(de_species3$intersected_markers))
                                                         #unique(unlist(de_species4$intersected_markers))#,unique(unlist(de_species5$intersected_markers)),
                                                         #unique(unlist(de_species6$intersected_markers)),unique(unlist(de_species7$intersected_markers))
  #))
  ort_match = data.frame(s1, s2)
  print(dim(ort_match))
  # only use group of genes
  if(!is.null(gene_filter)){
    ort_match = ort_match[ort_match$s1 %in% gene_filter[sp_pair[1]] &
                            ort_match$s2 %in% gene_filter[sp_pair[2]],]
  }
  # check if enough genes are being used
  if(nrow(ort_match)>0){
    message(paste0(sp_pair[1], " vs ", sp_pair[2], ": ", nrow(ort_match), " genes used."))
  } else{
    stop("No genes left to perform comparison. Revise gene filters or ortholog matching.")
  }
  
  if(filter_ct){ # use only ct present (will print those absent)
    if(any(!(ct_sets[[1]] %in% colnames(de_sp1$avg_exp)))){
      noct = ct_sets[[1]][!(ct_sets[[1]] %in% colnames(de_sp1$avg_exp))]
      message(paste0(c("Cell types absent in", sp_pair[1], ":", noct), collapse = " "))
    }
    ct_sets[[1]] = ct_sets[[1]][ct_sets[[1]] %in% colnames(de_sp1$avg_exp)]
    if(any(!(ct_sets[[2]] %in% colnames(de_sp2$avg_exp)))){
      noct = ct_sets[[2]][!(ct_sets[[2]] %in% colnames(de_sp2$avg_exp))]
      message(paste0(c("Cell types absent in", sp_pair[2], ":", noct), collapse = " "))
    }
    ct_sets[[2]] = ct_sets[[2]][ct_sets[[2]] %in% colnames(de_sp2$avg_exp)]
  }
  
  de_sp1$avg_exp <- de_sp1$avg_exp[intersect(rownames(de_sp1$avg_exp), ort_match$s1),]
  de_sp2$avg_exp <- de_sp2$avg_exp[intersect(rownames(de_sp2$avg_exp), ort_match$s2),]
  # normalise means
  m_list_sub = list(de_sp1$avg_exp, de_sp2$avg_exp)
  names(m_list_sub) = sp_pair
  if(norm_all){
    m_list_sub = lapply(m_list_sub[sp_pair],
                        function(x) t(apply(x, 1, function(y) y/mean(y))))
  } else{
    for(spi in 1:length(sp_pair)){
      sp = sp_pair[spi]
      if(length(ct_sets[[spi]])>1){
        m_list_sub[[sp]] = t(apply(m_list_sub[[sp]][,ct_sets[[sp]]],
                                   1, function(y) y/mean(y)))
      } else{
        m_list_sub[[sp]] = data.frame(m_list_sub[[sp]][,ct_sets[[sp]]])
        colnames(m_list_sub[[sp]]) = ct_sets[[sp]]
      }
      m_list_sub[[sp]][is.na(m_list_sub[[sp]])] = 0
    }
  }
  # get correlations, determine max per row and column
  sp1 = sp_pair[1]
  sp2 = sp_pair[2]
  sp1dat = m_list_sub[[sp1]][ort_match$s1,ct_sets[[sp1]],drop=FALSE]
  sp2dat = m_list_sub[[sp2]][ort_match$s2,ct_sets[[sp2]],drop=FALSE]
  if(any(c(dim(sp1dat), dim(sp2dat))==0)) stop("At least one species had no cell types to compare.")
  cort = psych::corr.test(sp1dat, sp2dat, method = "sp",
                          adjust = "fdr", alpha = 0.01, ci = F)
  # define most similar cell type
  if(is.null(sim_thr)){
    cort$maxrow = apply(cort$r, 1, which.max)
    cort$maxcol = apply(cort$r, 2, which.max)
  } else{
    mod_mat = cort$r
    mod_mat[mod_mat<sim_thr] = -3
    cort$maxrow = apply(mod_mat, 1, which.max)
    cort$maxcol = apply(mod_mat, 2, which.max)
  }
  
  cort[[paste0("sp1_", sp_pair[1])]] = sp1dat
  cort[[paste0("sp2_", sp_pair[2])]] = sp2dat
  return(cort)
}


############################### or if union go to --> ##################################
corrCellTypesPWT_Union = function(de_sp1, de_sp2, sp_pair = c("species1", "species2"),
                                 ort_scope = NULL, ct_sets = NULL, sim_thr = 0.2,
                                 gene_filter = NULL, norm_all = T, filter_ct = T){
  
  
  # which cells to use
  if(is.null(ct_sets)){
    ct_sets = list(unname(colnames(de_sp1$avg_exp)),
                   unname(colnames(de_sp2$avg_exp)))
    names(ct_sets) = sp_pair
  }
  # Intersect the combined markers with each de_species and save the results in s1 and s2
  s1 <- Reduce(union, lapply(de_species_list, function(de_species) {
    union_markers <- unlist(de_species$union_markers)
  }))

  s2 <- Reduce(union, lapply(de_species_list, function(de_species) {
    union_markers <- unlist(de_species$union_markers)
  }))
  # using marker intersection	
  #s1 = union(unique(unlist(de_species1$union_markers)),c(unique(unlist(de_species2$union_markers)),unique(unlist(de_species4$union_markers))
  #unique(unlist(de_species4$union_markers))#,unique(unlist(de_species5$union_markers)),
  #unique(unlist(de_species6$union_markers)),unique(unlist(de_species7$union_markers))
  #                                                      ))
  #s2 = union(unique(unlist(de_species1$union_markers)),c(unique(unlist(de_species2$union_markers)),unique(unlist(de_species4$union_markers))
  #unique(unlist(de_species4$union_markers))#,unique(unlist(de_species5$union_markers)),
  #unique(unlist(de_species6$union_markers)),unique(unlist(de_species7$union_markers))
  #                                                      ))
  ort_match = data.frame(s1, s2)	
  print(dim(ort_match))
 
  # only use group of genes
  if(!is.null(gene_filter)){
    ort_match = ort_match[ort_match$s1 %in% gene_filter[sp_pair[1]] &
                            ort_match$s2 %in% gene_filter[sp_pair[2]],]
  }
  
  # check if enough genes are being used
  if(nrow(ort_match)>0){
    message(paste0(sp_pair[1], " vs ", sp_pair[2], ": ", nrow(ort_match), " genes used."))
  } else{
    stop("No genes left to perform comparison. Revise gene filters or ortholog matching.")
  }
  
  if(filter_ct){ # use only ct present (will print those absent)
    if(any(!(ct_sets[[1]] %in% colnames(de_sp1$avg_exp)))){
      noct = ct_sets[[1]][!(ct_sets[[1]] %in% colnames(de_sp1$avg_exp))]
      message(paste0(c("Cell types absent in", sp_pair[1], ":", noct), collapse = " "))
    }
    ct_sets[[1]] = ct_sets[[1]][ct_sets[[1]] %in% colnames(de_sp1$avg_exp)]
    
    if(any(!(ct_sets[[2]] %in% colnames(de_sp2$avg_exp)))){
      noct = ct_sets[[2]][!(ct_sets[[2]] %in% colnames(de_sp2$avg_exp))]
      message(paste0(c("Cell types absent in", sp_pair[2], ":", noct), collapse = " "))
    }
    ct_sets[[2]] = ct_sets[[2]][ct_sets[[2]] %in% colnames(de_sp2$avg_exp)]
  }
  de_sp1$avg_exp <- de_sp1$avg_exp[intersect(rownames(de_sp1$avg_exp), ort_match$s1),]
  de_sp2$avg_exp <- de_sp2$avg_exp[intersect(rownames(de_sp2$avg_exp), ort_match$s2),]
  # normalise means
  m_list_sub = list(de_sp1$avg_exp, de_sp2$avg_exp)
  names(m_list_sub) = sp_pair
  if(norm_all){
    m_list_sub = lapply(m_list_sub[sp_pair],
                        function(x) t(apply(x, 1, function(y) y/mean(y))))
  } else{
    for(spi in 1:length(sp_pair)){
      sp = sp_pair[spi]
      if(length(ct_sets[[spi]])>1){
        m_list_sub[[sp]] = t(apply(m_list_sub[[sp]][,ct_sets[[sp]]],
                                   1, function(y) y/mean(y)))
      } else{
        m_list_sub[[sp]] = data.frame(m_list_sub[[sp]][,ct_sets[[sp]]])
        colnames(m_list_sub[[sp]]) = ct_sets[[sp]]
      }
      m_list_sub[[sp]][is.na(m_list_sub[[sp]])] = 0
    }
  }
  
  # get correlations, determine max per row and column
  sp1 = sp_pair[1]
  sp2 = sp_pair[2]
  sp1dat = m_list_sub[[sp1]][ort_match$s1,ct_sets[[sp1]],drop=FALSE]
  sp1dat[is.na(sp1dat)] <- 0
  sp2dat = m_list_sub[[sp2]][ort_match$s2,ct_sets[[sp2]],drop=FALSE]
  sp2dat[is.na(sp2dat)] <- 0
  missing_in_1 = rownames(sp2dat)[!rownames(sp2dat) %in% rownames(sp1dat)]
  if (length(missing_in_1) > 0){
    x = matrix(0.0, length(missing_in_1), ncol(sp1dat))
    rownames(x) = c(missing_in_1)
    sp1dat <- rbind(x, sp1dat)
  }
  missing_in_2 = rownames(sp1dat)[!rownames(sp1dat) %in% rownames(sp2dat)]
  if (length(missing_in_2) > 0){
    x = matrix(0.0, length(missing_in_2), ncol(sp2dat))
    rownames(x) = c(missing_in_2)
    sp2dat <- rbind(x, sp2dat)
  }
  sp1dat = sp1dat[order(row.names(sp1dat)), ]
  sp2dat = sp2dat[order(row.names(sp2dat)), ]
  if(any(c(dim(sp1dat), dim(sp2dat))==0)) stop("At least one species had no cell types to compare.")
  cort = psych::corr.test(sp1dat, sp2dat, method = "sp",
                          adjust = "fdr", alpha = 0.01, ci = F)
  if(is.null(sim_thr)){
    cort$maxrow = apply(cort$r, 1, which.max)
    cort$maxcol = apply(cort$r, 2, which.max)
  } else{
    mod_mat = cort$r
    mod_mat[mod_mat<sim_thr] = -3
    cort$maxrow = apply(mod_mat, 1, which.max)
    cort$maxcol = apply(mod_mat, 2, which.max)
  }
  
  cort[[paste0("sp1_", sp_pair[1])]] = sp1dat
  cort[[paste0("sp2_", sp_pair[2])]] = sp2dat
  return(cort)
}


######################## Find markers CORR ################################ THIS PART WILL WORK ! ###########################
calculatePWDET = function(objects, group.by, fc.thr = 0.2, pv.thr = 0.05){
  require(Seurat)
  require(parallel)
  library(viridis)

  # function to calculate the markers (in parallel)
  PWDET = function(mat, ds){
    # Find markers for current group versus all other groups
    mk = lapply(mat, function(x) FindMarkers(ds, only.pos = TRUE, ident.1 = x, verbose = T))
    names(mk) <- mat
    return(mk)
  }

  # Initialize list to store results
  de_list <- list()

  # Loop through all objects
  for (i in 1:length(objects)) {
    obj <- objects[[i]]

    # Get unique group names
    group_names = unique(obj@meta.data[,group.by])

    # calculate pairwise DE (parallelised)
    obj <- SetIdent(obj, value = group.by)
    mk_list = PWDET(group_names, obj)

    avg_exp = AverageExpression(obj, assays = "RNA", group.by = group.by)$RNA

    # Save results to list with consistent naming convention
    name <- paste0("de_species", i)
    de_list[[name]] <- list(avg_exp = avg_exp, markers = mk_list)
  }

  # Return list of results
  return(de_list)
}


intersection_union <- function(objects, group.by, fc.thr = 0.2, pv.thr = 0.05) {
  # calculate pairwise DE using calculatePWDET
  de_species_list <- calculatePWDET(objects, group.by, fc.thr, pv.thr)
  
  common_genes <- list()
  all_genes <- list()
  cell_types <- Reduce(intersect, lapply(de_species_list, function(x) names(x$markers)))
  
  for (ct in cell_types) {
    common_genes[[ct]] <- Reduce(intersect, lapply(de_species_list, function(x) rownames(x$markers[[ct]])))
    all_genes[[ct]] <- Reduce(union, lapply(de_species_list, function(x) rownames(x$markers[[ct]])))
  }
  
  for (i in 1:length(de_species_list)) {
    de_species_list[[i]]$intersected_markers <- common_genes
    de_species_list[[i]]$union_markers <- all_genes
  }
  
  return(de_species_list)
}
de_species_list <- intersection_union(list(specie1, specie2, specie3), "inter_species_ann", 0.2, 0.05)
sp_pair = c("human", "mouse")
corr = corrCellTypesPWT_Intersection(de_species_list$de_species1, de_species_list$de_species2, sp_pair = sp_pair)
plotCorr(corr, sp_pair, "intersect")
corr = corrCellTypesPWT_Union(de_species_list$de_species1, de_species_list$de_species2, sp_pair = sp_pair)
plotCorr(corr, sp_pair, "union")
