#' @title zigzag reference class.
#'
#' @description  This class contains all the tools for computing the posterior probability
#' each gene in a dataset is active given a set of libraries. Creating a zigzag object sets all the
#' hyperparameters and initializes all model parameters.
#'
#' @usage zigzag$new(data, gene_length = NULL, candidate_gene_list = "random",
#' num_active_components = 1,
#' weight_active_shape_1 = 2,
#' weight_active_shape_2 = 2,
#' inactive_means_prior_shape = 1,
#' inactive_means_prior_rate = 1/3,
#' inactive_variances_prior_min = 0.01,
#' inactive_variances_prior_max = 5,
#' spike_prior_shape_1 = 1,
#' spike_prior_shape_2 = 1,
#' active_means_dif_prior_shape = 1,
#' active_means_dif_prior_rate = 1/3,
#' active_variances_prior_min = 0.01,
#' active_variances_prior_max = 5,
#' shared_active_variances = TRUE,
#' output_directory = "output",
#' threshold_a = 0,
#' threshold_i = threshold_a[1],
#' beta = 1,
#' tau_shape = 1,
#' tau_rate = 1,
#' s0_mu = -1,
#' s0_sigma = 2,
#' s1_shape = 1,
#' s1_rate = 2,
#' alpha_r_shape = 1,
#' alpha_r_rate = 1/10, ...)
#'
#' @field Xg A G x R matrix containing the log-relative-expression of G genes in R libraries.
#' @field Yg A vector of length G containing the true latent log-relative-expression for G, genes.
#' @field output_directory A directory where burnin and and mcmc log files and plots are created and stored.
#' @field gene_names A character vector of length G contianing the gene names.
#' @field num_libraries The number of libraries in the data.
#' @field num_active_components The number of Normally-distributed subcomponents of the active distribution.
#' @field num_transcripts The number of genes or transcripts in the data depending on whether the data is gene or transcript level expression.
#' @field gene_lengths A vector of length G, containing either the mean lengths of gene transcripts or the lengths of transcripts.
#' @field threshold_a A vector containing the lower-bound of the prior distribution for the active component means.
#' If the length of the threshold_a vector is < than num_active_components then threshold_a[1] is used as the lower bound for the first active component
#' and subsequent active means are parameterized with orderd offsets rather than component-specific
#' lower boundaries. For example, if the user specifies 3 active subcomponents and sets threshold_a = c(1,3), then
#' the lower boundary for subcomponent 1 will be set 1, and the prior distributions for the difference between mean 1 and 2,
#' and the difference between mean 2 and 3 will be gamma distributed.
#' @field threshold_i The upper boundary for the inactive distribution.
#' @field candidate_gene_list A character vector of genes of interest for which mcmc log files for Yg and sigma_g should be output.
#' @field shared_active_variances A boolean indicating if all active subcomponents share a variance component. If FALSE, a variance parameter is estimated for each component.
#' @field sigma_g_upper_bound The upper boundary for the truncated log-normal prior for gene-specific variances. Default = Inf
#' @field weight_active_shape_1 .
#' @field weight_active_shape_2 .
#' @field inactive_means_prior_shape .
#' @field inactive_means_prior_rate .
#' @field inactive_variances_prior_min .
#' @field inactive_variances_prior_max .
#' @field inactive_variances_prior_shape .
#' @field inactive_variances_prior_rate .
#' @field inactive_variances_prior_log_min .
#' @field inactive_variances_prior_log_max .
#' @field spike_prior_shape_1 .
#' @field spike_prior_shape_2 .
#' @field active_means_dif_prior_shape .
#' @field active_means_dif_prior_rate .
#' @field active_variances_prior_min .
#' @field active_variances_prior_max .
#' @field active_variances_prior_shape .
#' @field active_variances_prior_rate .
#' @field active_variances_prior_log_min .
#' @field active_variances_prior_log_max .
#' @field weight_within_active_alpha .
#' @field s0_mu .
#' @field s0_sigma .
#' @field s1_rate .
#' @field s1_shape .
#' @field Sg .
#' @field s0 .
#' @field s1 .
#' @field tau_shape .
#' @field tau_rate .
#' @field tau .
#' @field sigma_g .
#' @field s0tau_trace .
#' @field alpha_r_shape .
#' @field alpha_r_rate .
#' @field p_x .
#' @field alpha_r .
#' @field inactive_spike_allocation .
#' @field candidate_gene_list .
#' @field allocation_active_inactive .
#' @field allocation_active_inactive_prob .
#' @field inactive_means .
#' @field inactive_variances .
#' @field spike_probability .
#' @field active_means_dif .
#' @field active_means .
#' @field active_variances .
#' @field weight_active .
#' @field weight_within_active .
#' @field allocation_within_active .
#' @field beta A power exponent for the likelihood function. If set to 0, then the mcmc samples from the joint prior.


zigzag <- setRefClass(

  Class = "zigzag",

  field = c(

    ####################
    # misc. properties #
    ####################

    gen = "numeric",

    field_list_values = "list",
    field_list_names = "list",

    sqrt2pi = "numeric",

    active_gene_set = "character",
    active_gene_set_idx = "numeric",

    Xg = "matrix",

    rwm = "numeric",
    rwv = "numeric",

    Yg = "numeric",
    Yg_proposed = "numeric",
    Yg_trace = "matrix",

    output_directory = "character",
    gene_names = "character",
    num_libraries = "numeric",
    num_active_components = "integer",
    component_matrix = "matrix",
    num_transcripts = "integer",
    gene_lengths = "numeric",

    multi_thresh_a = "logical",
    threshold_a = "numeric",
    threshold_i = "numeric",

    inf_tol = "numeric",
    model_paramlist = "character",
    lnl_trace = "list",

    inactive_mean_tuningParam = "numeric",
    inactive_variance_tuningParam = "numeric",
    spike_probability_tuningParam = "numeric",
    active_mean_tuningParam = "numeric",
    active_variance_tuningParam = "numeric",
    mixture_weight_tuningParam = "numeric",
    tuningParam_s0 = "numeric",
    tuningParam_s1 = "numeric",
    tuningParam_tau = "numeric",
    tuningParam_s0tau = "numeric",
    tuningParam_alpha_r ="numeric",
    tuningParam_yg = "numeric",
    tuningParam_sigma_g = "numeric",
    tuningParam_multi_sigma = "numeric",
    tuningParam_sigma_mu = "numeric",

    multi_sigma_trace = "list",
    sigma_mu_trace = "list",

    active_idx = "numeric",
    inactive_idx = "numeric",

    beta = "numeric",
    beta_oneLib = "numeric",

    temperature = "numeric",

    proposal_vector = "character",
    proposal_probs = "numeric",
    level_1_probs = "numeric",
    hieararchical_probs_oneLibrary = "numeric",
    combined_probs = "numeric",
    proposal_list = "list",

    ##########
    # priors #
    ##########

    weight_active_shape_1 = "numeric",
    weight_active_shape_2 = "numeric",

    inactive_means_prior_shape = "numeric",
    inactive_means_prior_rate = "numeric",
    inactive_variances_prior_min = "numeric",
    inactive_variances_prior_max = "numeric",
    inactive_variances_prior_shape = "numeric",
    inactive_variances_prior_rate = "numeric",
    inactive_variances_prior_log_min = "numeric",
    inactive_variances_prior_log_max = "numeric",


    spike_prior_shape_1 = "numeric",
    spike_prior_shape_2 = "numeric",

    active_means_dif_prior_shape = "numeric",
    active_means_dif_prior_rate = "numeric",
    active_variances_prior_min = "numeric",
    active_variances_prior_max = "numeric",
    active_variances_prior_shape = "numeric",
    active_variances_prior_rate = "numeric",
    active_variances_prior_log_min = "numeric",
    active_variances_prior_log_max = "numeric",


    weight_within_active_alpha = "numeric",

    s0_mu = "numeric",
    s0_sigma = "numeric",
    s1_rate = "numeric",
    s1_shape = "numeric",

    ##############
    # parameters #
    ##############



    #### hierarchical #####
    Sg = "numeric",
    s0 = "numeric",
    s0_trace = "list",

    s1 = "numeric",
    s1_trace = "list",

    tau_shape = "numeric",
    tau_rate = "numeric",
    tau = "numeric",
    tau_trace = "list",
    sigma_g = "numeric",
    sigma_g_trace = "matrix",
    sigma_g_upper_bound = "numeric",

    s0tau_trace = "list",

    alpha_r_shape = "numeric",
    alpha_r_rate = "numeric",
    p_x = "matrix",
    alpha_r = "numeric",
    alpha_r_trace = "list",

    # whether Yg is allocated to the spike value inf_tol
    no_detect_spike = "numeric",
    inactive_spike_allocation = "numeric",
    in_spike_idx = "numeric",
    out_spike_idx = "numeric",
    all_zero_idx = "numeric",

    candidate_gene_list = "character",

    XgLikelihood = "numeric",
    YgLikelihood = "numeric",
    sigma_g_probability = "numeric",

    #######################



    allocation_active_inactive = "integer",
    allocation_active_inactive_proposed = "integer",
    allocation_trace = "matrix",

    allocation_active_inactive_prob = "numeric",
    allocation_active_inactive_prob_proposed = "numeric",

    inactive_means = "numeric",
    inactive_means_proposed = "numeric",
    inactive_means_trace = "list",
    inactive_means_prob = "numeric",
    inactive_means_prob_proposed = "numeric",

    inactive_variances = "numeric",
    inactive_variances_proposed = "numeric",
    inactive_variances_trace = "list",
    inactive_variances_prob = "numeric",
    inactive_variances_prob_proposed = "numeric",

    spike_probability = "numeric",
    spike_probability_proposed = "numeric",
    spike_probability_trace = "list",
    spike_probability_prob = "numeric",
    spike_probability_prob_proposed = "numeric",

    active_means_dif = "numeric",
    active_means_dif_proposed = "numeric",
    active_means_trace = "list",
    active_means_dif_prob = "numeric",
    active_means_dif_prob_proposed = "numeric",
    active_means = "numeric",

    shared_active_variances = "logical",
    active_variances = "numeric",
    active_variances_proposed = "numeric",
    active_variances_trace = "list",
    active_variances_prob = "numeric",
    active_variances_prob_proposed = "numeric",

    # probability that a transcript is active
    weight_active = "numeric",
    weight_active_proposed = "numeric",
    weight_active_trace = "list",
    weight_active_prob = "numeric",
    weight_active_prob_proposed = "numeric",

    # probablity of an active transcript being in component k
    weight_within_active = "numeric",
    weight_within_active_proposed = "numeric",
    weight_within_active_trace = "list",
    weight_within_active_prob = "numeric",
    weight_within_active_prob_proposed = "numeric",
    mixture_weight_trace = "list",

    allocation_within_active = "list",
    allocation_within_active_proposed = "list"

  ) # end fields

)






