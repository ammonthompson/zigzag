#' zigzag reference class.
#'
#' @usage zigzag$new(data,
#' gene_length = NULL,
#' candidate_gene_list = "random",
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
#' multi_ta = FALSE,
#' threshold_a = 0,
#' threshold_i = threshold_a[1],
#' beta = 1,
#' temperature = 1,
#' tau_shape = 1,
#' tau_rate = 1,
#' s0_mu = -1,
#' s0_sigma = 2,
#' s1_shape = 1,
#' s1_rate = 2,
#' alpha_r_shape = 1,
#' alpha_r_rate = 1/10,
#' active_gene_set = NULL, ...)'
#' @field gen .
#' @field field_list_values .
#' @field field_list_names .
#' @field sqrt2pi .
#' @field active_gene_set .
#' @field active_gene_set_idx .
#' @field Xg .
#' @field Yg .
#' @field Yg_proposed .
#' @field Yg_trace .
#' @field output_directory .
#' @field gene_names .
#' @field num_libraries .
#' @field num_active_components .
#' @field component_matrix .
#' @field num_transcripts .
#' @field gene_lengths .
#' @field multi_thresh_a .
#' @field threshold_a .
#' @field threshold_i .
#' @field inf_tol .
#' @field model_paramlist .
#' @field lnl_trace .
#' @field inactive_mean_tuningParam .
#' @field inactive_variance_tuningParam .
#' @field spike_probability_tuningParam .
#' @field active_mean_tuningParam .
#' @field active_variance_tuningParam .
#' @field mixture_weight_tuningParam .
#' @field tuningParam_s0 .
#' @field tuningParam_s1 .
#' @field tuningParam_tau .
#' @field tuningParam_s0tau .
#' @field tuningParam_alpha_r .
#' @field tuningParam_yg .
#' @field tuningParam_sigma_g .
#' @field tuningParam_multi_sigma .
#' @field tuningParam_sigma_mu .
#' @field multi_sigma_trace .
#' @field sigma_mu_trace .
#' @field active_idx .
#' @field inactive_idx .
#' @field beta .
#' @field beta_oneLib .
#' @field temperature .
#' @field proposal_vector .
#' @field proposal_probs .
#' @field level_1_probs .
#' @field hieararchical_probs_oneLibrary .
#' @field combined_probs .
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
#' @field s0_trace .
#' @field s1 .
#' @field s1_trace .
#' @field tau_shape .
#' @field tau_rate .
#' @field tau .
#' @field tau_trace .
#' @field sigma_g .
#' @field sigma_g_trace .
#' @field s0tau_trace .
#' @field alpha_r_shape .
#' @field alpha_r_rate .
#' @field p_x .
#' @field alpha_r .
#' @field alpha_r_trace .
#' @field all_means .
#' @field all_variances .
#' @field no_detect_spike .
#' @field inactive_spike_allocation .
#' @field in_spike_idx .
#' @field out_spike_idx .
#' @field all_zero_idx .
#' @field candidate_gene_list .
#' @field XgLikelihood .
#' @field YgLikelihood .
#' @field sigma_g_probability .
#' @field allocation_active_inactive .
#' @field allocation_active_inactive_proposed .
#' @field allocation_trace .
#' @field allocation_active_inactive_prob .
#' @field allocation_active_inactive_prob_proposed .
#' @field inactive_means .
#' @field inactive_means_proposed .
#' @field inactive_means_trace .
#' @field inactive_means_prob .
#' @field inactive_means_prob_proposed .
#' @field inactive_variances .
#' @field inactive_variances_proposed .
#' @field inactive_variances_trace .
#' @field inactive_variances_prob .
#' @field inactive_variances_prob_proposed .
#' @field spike_probability .
#' @field spike_probability_proposed .
#' @field spike_probability_trace .
#' @field spike_probability_prob .
#' @field spike_probability_prob_proposed .
#' @field active_means_dif .
#' @field active_means_dif_proposed .
#' @field active_means_trace .
#' @field active_means_dif_prob .
#' @field active_means_dif_prob_proposed .
#' @field active_means .
#' @field shared_active_variances .
#' @field active_variances .
#' @field active_variances_proposed .
#' @field active_variances_trace .
#' @field active_variances_prob .
#' @field active_variances_prob_proposed .
#' @field weight_active .
#' @field weight_active_proposed .
#' @field weight_active_trace .
#' @field weight_active_prob .
#' @field weight_active_prob_proposed .
#' @field weight_within_active .
#' @field weight_within_active_proposed .
#' @field weight_within_active_trace .
#' @field weight_within_active_prob .
#' @field weight_within_active_prob_proposed .
#' @field mixture_weight_trace .
#' @field allocation_within_active .
#' @field allocation_within_active_proposed .

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
    #sim_xg = "list",

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

    s0tau_trace = "list",

    alpha_r_shape = "numeric",
    alpha_r_rate = "numeric",
    p_x = "matrix",
    alpha_r = "numeric",
    alpha_r_trace = "list",

    all_means = "list",
    all_variances = "list",

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






