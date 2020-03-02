zigzag$methods(

  initializeOutputFiles = function(prefix, post_pred_boolean = FALSE){

    model_paramlist <<- c("gen", "lnl", "s0", "s1", "tau")

    for(param in seq(num_libraries)){

      model_paramlist[length(model_paramlist)+1] <<- paste0("alpha_r_library_", param)

    }

    model_paramlist <<- c(model_paramlist, "spike_probability", "inactive_mean", "inactive_variance")

    for(param in c("active_mean","active_variance") ){

      for(k in 1:num_active_components){
        model_paramlist[length(model_paramlist)+1] <<- paste0(param, "_component",k)
      }

    }

    model_paramlist[length(model_paramlist)+1] <<- "weight_active"

    for(k in 1:num_active_components ){
      model_paramlist[length(model_paramlist)+1] <<- paste0("weight_within_active_component_",k)
    }

    write.table(matrix(model_paramlist,nrow=1),
                file=paste0(output_directory, "/", prefix,"_model_parameters.log"), col.names = F, row.names = F, sep="\t", quote = F)

    ## posterior for Yg of candidate genes

    write.table(matrix(c("gen", candidate_gene_list), nrow = 1),
                file=paste0(output_directory, "/", prefix, "_yg_candidate_genes.log"), col.names = F, row.names = F, sep = "\t", quote = F)

    ## posterior for sigma_g of candidate genes

    write.table(matrix(c("gen", candidate_gene_list), nrow = 1),
                file=paste0(output_directory, "/", prefix, "_sigmag_candidate_genes.log"), col.names = F, row.names = F, sep = "\t", quote = F)

    ## posterior predicitve simulation file
    if(post_pred_boolean){

      write.table(matrix(c("gen", "Wass_Lower",  "Wass_Upper", "Rums"), nrow = 1),
                  file = paste0(output_directory, "/", prefix, ".post_predictive_output.log"),
                  col.names = F, row.names = F, sep = "\t", quote = F)

    }

  },

  writeToOutputFiles = function(prefix, gen){

    # Write parameter posterior samples to file
    Ygrow=c(lnl_trace[[length(lnl_trace)]], s0, s1, tau, alpha_r,
            spike_probability, inactive_means, inactive_variances)

    for(k in 1:num_active_components ){
      Ygrow[length(Ygrow)+1] <- active_means[k]
    }


    for(k in 1:num_active_components ){
      Ygrow[length(Ygrow)+1] <- active_variances[k]
    }

    Ygrow[length(Ygrow)+1] <- weight_active

    for(k in 1:num_active_components ){
      Ygrow[length(Ygrow)+1] <- weight_within_active[k]
    }

    Ygrow=round(as.numeric(Ygrow),digits=4)
    write.table(matrix(c(gen, Ygrow),nrow=1),
                file=paste0(output_directory, "/", prefix,"_model_parameters.log"),
                append=T, sep="\t",row.names=F, col.names=F)

  },

  writeToYgSigmagOutputFiles = function(prefix, gen){

    ### write candidate gene Yg posterior samples to file
    geneExp <- Yg[match(candidate_gene_list, gene_names)] *
      (1-inactive_spike_allocation[match(candidate_gene_list, gene_names)]) +
      -10 * inactive_spike_allocation[match(candidate_gene_list, gene_names)]

    write.table(matrix(c(gen, round(geneExp, digits = 6)), nrow=1),
                file=paste0(output_directory, "/", prefix, "_yg_candidate_genes.log"),
                append=T, sep="\t",row.names=F, col.names=F)


    ### write candidate gene sigma_g posterior samples to file
    geneSigma_g <- sigma_g[match(candidate_gene_list, gene_names)] *
      (1-inactive_spike_allocation[match(candidate_gene_list, gene_names)]) +
      0 * inactive_spike_allocation[match(candidate_gene_list, gene_names)]

    write.table(matrix(c(gen, round(geneSigma_g, digits = 6)), nrow=1),
                file=paste0(output_directory, "/", prefix, "_sigmag_candidate_genes.log"),
                append=T, sep="\t",row.names=F, col.names=F)

  },

  computeGeneExpressionProbs_writeToFile = function(prefix){

    write.table("prob_active", file=paste0(output_directory, "/", prefix, "_probability_active.tab"), quote = F, row.names = "gene", col.names = F, sep="\t")


    probs = round(allocation_trace/sum(allocation_trace[1,]), digits = 3)

    write.table(format(1 - probs[,1], digits = 3),file=paste0(output_directory, "/", prefix, "_probability_active.tab"), append = T, quote = F,
                row.names=gene_names, col.names = F, sep="\t")

    write.table(matrix(c(sapply(0:num_active_components,function(k){
      return(paste0("prob_",k))})), nrow =1 ,byrow=T), file=paste0(output_directory, "/", prefix, "_probability_in_component.tab"),
      row.names = "gene", quote = F, col.names = F, sep="\t")

    write.table(format(probs, digits = 3),file=paste0(output_directory, "/", prefix, "_probability_in_component.tab"),
                append= T, quote = F, row.names=gene_names, col.names = F, sep="\t")

  }

)
