######################################################################################################################
##  This is the code for the simulation study of the manuscript 'Power prior models for treatment estimation        ##
##  in a small n, sequential, multiple assignment, randomized trial'                                                ##
##                                                                                                                  ##
##  Section 1: All functions for the simulation study                                                               ##
##  Section 2: The assignment of the parameters used in the simulation study                                        ##
##  Section 3: Running the simulation study and obtaining the output table                                          ##                                                                                          ##
##  Section 4: Giving outputs of estimation of delta and pi and a histogram of delta                                ##
##                                                                                                                  ##
##  The code was last run on December 12, 2021, and the version of R was 4.0.3.                                    ##
######################################################################################################################


library(MCMCpack)
library(coda)
library(tidyverse)
library(rjags)
library(optimx)

############################################################################
############################################################################
###                                                                      ###
###                         SECTION 1: FUNCTIONS                         ###
###                                                                      ###
############################################################################
############################################################################

# Generate dataset for SMART
standard_SMART <- function(pi_1A, pi_1B, pi_1C, discount_y, discount_n1, discount_n2, n_A, n_B, n_C){
  pi_2A_y <- pi_1A * discount_y[1]        # Second stage response rate of responders to A 
  pi_2B.A_n <- pi_1A * discount_n1[2]     # Second stage response rate of non-responders to B who receive A in the second stage
  pi_2C.A_n <- pi_1A * discount_n1[3]     # Second stage response rate of non-responders to C who receive A in the second stage
  pi_2B_y <- pi_1B * discount_y[2]        # Second stage response rate of responders to B
  pi_2A.B_n <- pi_1B * discount_n1[1]     # Second stage response rate of non-responders to A who receive B in the second stage
  pi_2C.B_n <- pi_1B * discount_n2[3]     # Second stage response rate of non-responders to C who receive B in the second stage
  pi_2C_y <- pi_1C * discount_y[3]        # Second stage response rate of responders to C
  pi_2A.C_n <- pi_1C * discount_n2[1]     # Second stage response rate of non-responders to A who receive C in the second stage
  pi_2B.C_n <- pi_1C * discount_n2[2]     # Second stage response rate of non-responders to B who receive C in the second stage
  ### stage 1 ###
  # trt A
  n_resp_A <- rbinom(n=1,size=n_A,prob=pi_1A)
  n_nonresp_A <- n_A - n_resp_A
  n_nonresp_A_stage2_B <- ifelse(n_nonresp_A%%2==1,(n_nonresp_A+1)/2,n_nonresp_A/2)
  n_nonresp_A_stage2_C <- n_nonresp_A - n_nonresp_A_stage2_B
  
  # trt B
  n_resp_B <- rbinom(n=1,size=n_B,prob=pi_1B)
  n_nonresp_B <- n_B - n_resp_B
  n_nonresp_B_stage2_A <<- ifelse(n_nonresp_B%%2==1,(n_nonresp_B+1)/2,n_nonresp_B/2)
  n_nonresp_B_stage2_C <<- n_nonresp_B - n_nonresp_B_stage2_A
  
  # trt C
  n_resp_C <- rbinom(n=1,size=n_C,prob=pi_1C)
  n_nonresp_C <- n_C - n_resp_C
  n_nonresp_C_stage2_A <<- ifelse(n_nonresp_C%%2==1,(n_nonresp_C+1)/2,n_nonresp_C/2)
  n_nonresp_C_stage2_B <<- n_nonresp_C - n_nonresp_C_stage2_A
  
  ### stage 2 ###
  # stage I = trt A
  n_resp_A.A.Y <- rbinom(n=1,size=n_resp_A,prob=pi_2A_y)
  data_A.A.Y <- data.frame(treatment_stageI = rep(1,n_resp_A),
                           response_stageI = rep(1,n_resp_A),
                           treatment_stageII = rep(1,n_resp_A),
                           response_stageII = c(rep(1,n_resp_A.A.Y),rep(0,n_resp_A-n_resp_A.A.Y)))
  
  n_resp_A.B.Y <- rbinom(n=1,size=n_nonresp_A_stage2_B,prob=pi_2A.B_n)
  data_A.B.Y <- data.frame(treatment_stageI = rep(1,n_nonresp_A_stage2_B),
                           response_stageI = rep(0,n_nonresp_A_stage2_B),
                           treatment_stageII = rep(2,n_nonresp_A_stage2_B),
                           response_stageII = c(rep(1,n_resp_A.B.Y),rep(0,n_nonresp_A_stage2_B-n_resp_A.B.Y)))
  
  n_resp_A.C.Y <- rbinom(n=1,size=n_nonresp_A_stage2_C,prob=pi_2A.C_n)
  data_A.C.Y <- data.frame(treatment_stageI = rep(1,n_nonresp_A_stage2_C),
                           response_stageI = rep(0,n_nonresp_A_stage2_C),
                           treatment_stageII = rep(3,n_nonresp_A_stage2_C),
                           response_stageII = c(rep(1,n_resp_A.C.Y),rep(0,n_nonresp_A_stage2_C-n_resp_A.C.Y)))
  
  # stage I = trt B
  n_resp_B.B.Y <- rbinom(n=1,size=n_resp_B,prob=pi_2B_y)
  data_B.B.Y <- data.frame(treatment_stageI = rep(2,n_resp_B),
                           response_stageI = rep(1,n_resp_B),
                           treatment_stageII = rep(2,n_resp_B),
                           response_stageII = c(rep(1,n_resp_B.B.Y),rep(0,n_resp_B-n_resp_B.B.Y)))
  
  n_resp_B.A.Y <- rbinom(n=1,size=n_nonresp_B_stage2_A,prob=pi_2B.A_n)
  data_B.A.Y <- data.frame(treatment_stageI = rep(2,n_nonresp_B_stage2_A),
                           response_stageI = rep(0,n_nonresp_B_stage2_A),
                           treatment_stageII = rep(1,n_nonresp_B_stage2_A),
                           response_stageII = c(rep(1,n_resp_B.A.Y),rep(0,n_nonresp_B_stage2_A-n_resp_B.A.Y)))
  
  n_resp_B.C.Y <- rbinom(n=1,size=n_nonresp_B_stage2_C,prob=pi_2B.C_n)
  data_B.C.Y <- data.frame(treatment_stageI = rep(2,n_nonresp_B_stage2_C),
                           response_stageI = rep(0,n_nonresp_B_stage2_C),
                           treatment_stageII = rep(3,n_nonresp_B_stage2_C),
                           response_stageII = c(rep(1,n_resp_B.C.Y),rep(0,n_nonresp_B_stage2_C-n_resp_B.C.Y)))
  
  # stage I = trt C
  n_resp_C.C.Y <- rbinom(n=1,size=n_resp_C,prob=pi_2C_y)
  data_C.C.Y <- data.frame(treatment_stageI = rep(3,n_resp_C),
                           response_stageI = rep(1,n_resp_C),
                           treatment_stageII = rep(3,n_resp_C),
                           response_stageII = c(rep(1,n_resp_C.C.Y),rep(0,n_resp_C-n_resp_C.C.Y)))
  
  n_resp_C.A.Y <- rbinom(n=1,size=n_nonresp_C_stage2_A,prob=pi_2C.A_n)
  data_C.A.Y <- data.frame(treatment_stageI = rep(3,n_nonresp_C_stage2_A),
                           response_stageI = rep(0,n_nonresp_C_stage2_A),
                           treatment_stageII = rep(1,n_nonresp_C_stage2_A),
                           response_stageII = c(rep(1,n_resp_C.A.Y),rep(0,n_nonresp_C_stage2_A-n_resp_C.A.Y)))
  
  n_resp_C.B.Y <- rbinom(n=1,size=n_nonresp_C_stage2_B,prob=pi_2C.B_n)
  data_C.B.Y <- data.frame(treatment_stageI = rep(3,n_nonresp_C_stage2_B),
                           response_stageI = rep(0,n_nonresp_C_stage2_B),
                           treatment_stageII = rep(2,n_nonresp_C_stage2_B),
                           response_stageII = c(rep(1,n_resp_C.B.Y),rep(0,n_nonresp_C_stage2_B-n_resp_C.B.Y)))
  
  data_stageI.II <- rbind(data_A.A.Y,data_A.B.Y,data_A.C.Y,
                          data_B.B.Y,data_B.A.Y,data_B.C.Y,
                          data_C.C.Y,data_C.A.Y,data_C.B.Y)
  return(list(dataset = data_stageI.II, n_resp_I = c(n_resp_A, n_resp_B, n_resp_C),
              n_II = c(n_resp_A, n_nonresp_A_stage2_B, n_nonresp_A_stage2_C,
                       n_resp_B, n_nonresp_B_stage2_A, n_nonresp_B_stage2_C,
                       n_resp_C, n_nonresp_C_stage2_A, n_nonresp_C_stage2_B),
              n_resp_II = c(n_resp_A.A.Y, n_resp_A.B.Y, n_resp_A.C.Y,
                            n_resp_B.B.Y, n_resp_B.A.Y, n_resp_B.C.Y,
                            n_resp_C.C.Y, n_resp_C.A.Y, n_resp_C.B.Y)))
} 

# log likelihood + prior for modified power prior model (weights as random variables)
complete_likelihood_modified_power_prior_2_delta <-
  function(parameters, s2_n_group, s2_n_resp, s1_n_group, s1_n_resp,
           prior_delta_a, prior_delta_b, prior_pi_a, prior_pi_b) {
    
    # Assign the indices of subgroups of trts
    stage2_A <- c(1,5,8)
    stage2_B <- c(4,2,9)
    stage2_C <- c(7,3,6)
    
    # Specify pis and deltas
    pii <- parameters[1:3]
    deltas <- parameters[4:5]
    if (all(pii < 1 & pii > 0) & all(deltas > 0 & deltas < 1)) {
      log_prior0 <- 
        (sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_A]) + prior_pi_a - 1) * log(pii[1]) + 
        (sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_A] - s2_n_resp[stage2_A])) + 
           prior_pi_b - 1) * log(1 - pii[1]) +
        (sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_B]) + prior_pi_a - 1) * log(pii[2]) + 
        (sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_B] - s2_n_resp[stage2_B])) + 
           prior_pi_b - 1) * log(1 - pii[2]) +
        (sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_C]) + prior_pi_a - 1) * log(pii[3]) + 
        (sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_C] - s2_n_resp[stage2_C])) + 
           prior_pi_b - 1) * log(1 - pii[3]) + 
        sum(log(dbeta(deltas, prior_delta_a, prior_delta_b))) -
        lbeta(sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_A]) + prior_pi_a, 
              sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_A] - s2_n_resp[stage2_A])) + prior_pi_b) - 
        lbeta(sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_B]) + prior_pi_a, 
              sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_B] - s2_n_resp[stage2_B])) + prior_pi_b) - 
        lbeta(sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_C]) + prior_pi_a, 
              sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_C] - s2_n_resp[stage2_C])) + prior_pi_b)
      
      log_lik <- 
        s1_n_resp[1] * log(pii[1]) + (s1_n_group[1] - s1_n_resp[1]) * log(1 - pii[1]) +
        s1_n_resp[2] * log(pii[2]) + (s1_n_group[2] - s1_n_resp[2]) * log(1 - pii[2]) + 
        s1_n_resp[3] * log(pii[3]) + (s1_n_group[3] - s1_n_resp[3]) * log(1 - pii[3])
      
    } else {
      log_prior0 <- -10^200
      log_lik <- -10^200
    }
    log_lik_prior0 <- log_lik + log_prior0
    return(log_lik_prior0)
  }

# normalized marginal log likelihood for the power prior model with the marginal likelihood criteria  
complete_marginal_log_likelihood_2_delta <- 
  function(delta, s2_n_group, s2_n_resp, s1_n_group, s1_n_resp,
           prior_pi_a, prior_pi_b) {
    
    # Assign the indices of subgroups of trts
    stage2_A <- c(1,5,8)
    stage2_B <- c(4,2,9)
    stage2_C <- c(7,3,6)  
    
    loglik <-
      lbeta(s1_n_resp[1] + sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_A]) + prior_pi_a, 
            s1_n_group[1] - s1_n_resp[1] + sum(delta[c(1, 2, 2)] * 
                                                 (s2_n_group[stage2_A] - s2_n_resp[stage2_A])) + prior_pi_b) +
      lbeta(s1_n_resp[2] + sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_B]) + prior_pi_a, 
            s1_n_group[2] - s1_n_resp[2] + sum(delta[c(1, 2, 2)] * 
                                                 (s2_n_group[stage2_B] - s2_n_resp[stage2_B])) + prior_pi_b) +
      lbeta(s1_n_resp[3] + sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_C]) + prior_pi_a, 
            s1_n_group[3] - s1_n_resp[3] + sum(delta[c(1, 2, 2)] * 
                                                 (s2_n_group[stage2_C] - s2_n_resp[stage2_C])) + prior_pi_b) -
      lbeta(sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_A]) + prior_pi_a, 
            sum(delta[c(1, 2, 2)] * (s2_n_group[stage2_A] - s2_n_resp[stage2_A])) + prior_pi_b) -
      lbeta(sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_B]) + prior_pi_a, 
            sum(delta[c(1, 2, 2)] * (s2_n_group[stage2_B] - s2_n_resp[stage2_B])) + prior_pi_b) -
      lbeta(sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_C]) + prior_pi_a, 
            sum(delta[c(1, 2, 2)] * (s2_n_group[stage2_C] - s2_n_resp[stage2_C])) + prior_pi_b)
    return(-loglik)
  }

# penalized marginal log likelihood (unnormalized)
complete_penalized_marginal_log_likelihood_2_delta <- 
  function(delta, s2_n_group, s2_n_resp, s1_n_group, s1_n_resp,
           prior_pi_a, prior_pi_b) {
    
    # Assign the indices of subgroups of trts
    stage2_A <- c(1,5,8)
    stage2_B <- c(4,2,9)
    stage2_C <- c(7,3,6)  
    
    loglik <-
      lbeta(s1_n_resp[1] + sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_A]) + prior_pi_a, 
            s1_n_group[1] - s1_n_resp[1] + sum(delta[c(1, 2, 2)] * 
                                                 (s2_n_group[stage2_A] - s2_n_resp[stage2_A])) + prior_pi_b) +
      lbeta(s1_n_resp[2] + sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_B]) + prior_pi_a, 
            s1_n_group[2] - s1_n_resp[2] + sum(delta[c(1, 2, 2)] * 
                                                 (s2_n_group[stage2_B] - s2_n_resp[stage2_B])) + prior_pi_b) +
      lbeta(s1_n_resp[3] + sum(delta[c(1, 2, 2)] * s2_n_resp[stage2_C]) + prior_pi_a, 
            s1_n_group[3] - s1_n_resp[3] + sum(delta[c(1, 2, 2)] * 
                                                 (s2_n_group[stage2_C] - s2_n_resp[stage2_C])) + prior_pi_b)
    # number of subjects assigned to each subgroup 
    n_subgroup <- c(sum(s2_n_group[c(1, 4, 7)]), sum(s2_n_group[c(2:3, 5:6, 8:9)]))
    
    penalized_loglik <- -2 * loglik + sum(log(n_subgroup) / delta)
    
    return(penalized_loglik)
  }

# log likelihood + prior for power prior model (weights as fixed numbers)

complete_likelihood_power_prior_2_delta <- 
  function(pii, deltas, s2_n_group, s2_n_resp, s1_n_group, s1_n_resp,
           prior_pi_a, prior_pi_b) {
    
    # Assign the indices of subgroups of trts
    stage2_A <- c(1,5,8)
    stage2_B <- c(4,2,9)
    stage2_C <- c(7,3,6)  
    
    if (all(pii < 1 & pii > 0)) {
      log_prior0 <- 
        (sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_A]) + prior_pi_a - 1) * log(pii[1]) +
        (sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_A] - s2_n_resp[stage2_A])) + 
           prior_pi_b - 1) * log(1 - pii[1]) + 
        (sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_B]) + prior_pi_a - 1) * log(pii[2]) +
        (sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_B] - s2_n_resp[stage2_B])) + 
           prior_pi_b - 1) * log(1 - pii[2]) + 
        (sum(deltas[c(1, 2, 2)] * s2_n_resp[stage2_C]) + prior_pi_a - 1) * log(pii[3]) +
        (sum(deltas[c(1, 2, 2)] * (s2_n_group[stage2_C] - s2_n_resp[stage2_C])) + 
           prior_pi_b - 1) * log(1 - pii[3])
      log_lik <- 
        s1_n_resp[1] * log(pii[1]) + (s1_n_group[1] - s1_n_resp[1]) * log(1 - pii[1]) +
        s1_n_resp[2] * log(pii[2]) + (s1_n_group[2] - s1_n_resp[2]) * log(1 - pii[2]) +
        s1_n_resp[3] * log(pii[3]) + (s1_n_group[3] - s1_n_resp[3]) * log(1 - pii[3])
    } else {
      log_prior0 <- -10^200
      log_lik <- -10^200
    }
    log_lik_prior0 <- log_lik + log_prior0
    return(log_lik_prior0)
  }

# FET with 2 deltas in total
complete_FET_2_delta <- 
  function(s2_n_group, s2_n_resp, s1_n_group, s1_n_resp) {
    
    # Assign the indices of subgroups of trts
    stage2_A <- c(1,5,8)
    stage2_B <- c(4,2,9)
    stage2_C <- c(7,3,6)  
    
    FT_matrix1_1 <- matrix(c(s1_n_group[1] - s1_n_resp[1], s1_n_resp[1], 
                             s1_n_resp[1] - s2_n_resp[stage2_A[1]], s2_n_resp[stage2_A[1]]),
                           nrow = 2, ncol = 2)
    FT_matrix1_2 <- matrix(c(s1_n_group[2] - s1_n_resp[2], s1_n_resp[2], 
                             s1_n_resp[2] - s2_n_resp[stage2_B[1]], s2_n_resp[stage2_B[1]]),
                           nrow = 2, ncol = 2)
    FT_matrix1_3 <- matrix(c(s1_n_group[3] - s1_n_resp[3], s1_n_resp[3], 
                             s1_n_resp[3] - s2_n_resp[stage2_C[1]], s2_n_resp[stage2_C[1]]),
                           nrow = 2, ncol = 2)
    FT_matrix2_1 <- matrix(c(s1_n_group[1] - s1_n_resp[1], s1_n_resp[1], 
                             sum(s2_n_group[stage2_A[2:3]]) - sum(s2_n_resp[stage2_A[2:3]]), sum(s2_n_resp[stage2_A[2:3]])),
                           nrow = 2, ncol = 2)
    FT_matrix2_2 <- matrix(c(s1_n_group[2] - s1_n_resp[2], s1_n_resp[2], 
                             sum(s2_n_group[stage2_B[2:3]]) - sum(s2_n_resp[stage2_B[2:3]]), sum(s2_n_resp[stage2_B[2:3]])),
                           nrow = 2, ncol = 2)
    FT_matrix2_3 <- matrix(c(s1_n_group[3] - s1_n_resp[3], s1_n_resp[3], 
                             sum(s2_n_group[stage2_C[2:3]]) - sum(s2_n_resp[stage2_C[2:3]]), sum(s2_n_resp[stage2_C[2:3]])),
                           nrow = 2, ncol = 2)
    delta_1 <- c(fisher.test(FT_matrix1_1)$p.value,
                 fisher.test(FT_matrix1_2)$p.value,
                 fisher.test(FT_matrix1_3)$p.value) %>% 
      mean
    delta_2 <- c(fisher.test(FT_matrix2_1)$p.value,
                 fisher.test(FT_matrix2_2)$p.value,
                 fisher.test(FT_matrix2_3)$p.value) %>% 
      mean
    return(c(delta_1, delta_2))
  }

# Bhattacharyya's overlap measure with 2 deltas in total
complete_BM_2_delta <- 
  function(s2_n_group, s2_n_resp, s1_n_group, s1_n_resp, prior_pi_a, prior_pi_b) {
    
    # Assign the indices of subgroups of trts
    stage2_A <- c(1,5,8)
    stage2_B <- c(4,2,9)
    stage2_C <- c(7,3,6)  
    
    a1_A <- s1_n_resp[1] + prior_pi_a
    a1_B <- s1_n_resp[2] + prior_pi_a
    a1_C <- s1_n_resp[3] + prior_pi_a
    b1_A <- s1_n_group[1] - s1_n_resp[1] +  prior_pi_b
    b1_B <- s1_n_group[2] - s1_n_resp[2] +  prior_pi_b
    b1_C <- s1_n_group[3] - s1_n_resp[3] +  prior_pi_b
    a2_A_r <- s2_n_resp[stage2_A[1]] + prior_pi_a
    a2_B_r <- s2_n_resp[stage2_B[1]] + prior_pi_a
    a2_C_r <- s2_n_resp[stage2_C[1]] + prior_pi_a
    b2_A_r <- s2_n_group[stage2_A[1]] - s2_n_resp[stage2_A[1]] + prior_pi_b
    b2_B_r <- s2_n_group[stage2_B[1]] - s2_n_resp[stage2_B[1]] + prior_pi_b
    b2_C_r <- s2_n_group[stage2_C[1]] - s2_n_resp[stage2_C[1]] + prior_pi_b
    a2_A_nr <- sum(s2_n_resp[stage2_A[2:3]]) + prior_pi_a
    a2_B_nr <- sum(s2_n_resp[stage2_B[2:3]]) + prior_pi_a
    a2_C_nr <- sum(s2_n_resp[stage2_C[2:3]]) + prior_pi_a
    b2_A_nr <- sum(s2_n_group[stage2_A[2:3]] - s2_n_resp[stage2_A[2:3]]) + prior_pi_b
    b2_B_nr <- sum(s2_n_group[stage2_B[2:3]] - s2_n_resp[stage2_B[2:3]]) + prior_pi_b
    b2_C_nr <- sum(s2_n_group[stage2_C[2:3]] - s2_n_resp[stage2_C[2:3]]) + prior_pi_b
    
    BM_A_r <- beta((a1_A + a2_A_r) / 2, (b1_A + b2_A_r) / 2) / sqrt(beta(a1_A, b1_A) * beta(a2_A_r, b2_A_r))
    BM_B_r <- beta((a1_B + a2_B_r) / 2, (b1_B + b2_B_r) / 2) / sqrt(beta(a1_B, b1_B) * beta(a2_B_r, b2_B_r))
    BM_C_r <- beta((a1_C + a2_C_r) / 2, (b1_C + b2_C_r) / 2) / sqrt(beta(a1_C, b1_C) * beta(a2_C_r, b2_C_r))
    BM_A_nr <- beta((a1_A + a2_A_nr) / 2, (b1_A + b2_A_nr) / 2) / sqrt(beta(a1_A, b1_A) * beta(a2_A_nr, b2_A_nr))
    BM_B_nr <- beta((a1_B + a2_B_nr) / 2, (b1_B + b2_B_nr) / 2) / sqrt(beta(a1_B, b1_B) * beta(a2_B_nr, b2_B_nr))
    BM_C_nr <- beta((a1_C + a2_C_nr) / 2, (b1_C + b2_C_nr) / 2) / sqrt(beta(a1_C, b1_C) * beta(a2_C_nr, b2_C_nr))
    delta_1 <- mean(BM_A_r, BM_B_r, BM_C_r)
    delta_2 <- mean(BM_A_nr, BM_B_nr, BM_C_nr)
    
    return(c(delta_1, delta_2))
  }

# Help generate tables to plot Figure 4 in the manuscript
bias_rMSE_organizer <- function(result_table){
  order_final <- c("MPP", "pen", "EB", "BM", "FET", "BJSM", "d1", "d0")
  bias_rMSE <- result_table %>% 
    dplyr::filter(grepl("pi_.*_MPP|pi_.*_pen|pi_.*_EB|pi_.*_d1|pi_.*_d0|pi_.*_BJSM|pi_.*_BM|pi_.*_FET", type)) %>% 
    separate(type, into = c("pi", "arm", "method"), sep = "_") %>% 
    arrange(match(method, order_new), desc(method), arm) %>% 
    dplyr::select(arm, method, bias, rMSE) %>% 
    pivot_wider(names_from = method, values_from = c(bias, rMSE)) %>% 
    dplyr::select(matches(order_final)) %>% 
    summarise_all(mean)
  return(bias_rMSE)
}
###########################################################################
###########################################################################
###                                                                     ###
###                     SECTION 2: PARAMETERS SETUP                     ###
###                                                                     ###
###########################################################################
###########################################################################
jags_path <- "Your_path_to_bugs_code"

##### The true response rates of treatments for scenarios 1 - 7 

### Scenario 1
# The pre-specified response rate of treatment A in the first stage
pi_1A <- 0.2        
# The pre-specified response rate of treatment B in the first stage
pi_1B <- 0.3          
# The pre-specified response rate of treatment C in the first stage
pi_1C <- 0.4    
# The pre-specified values of linkage parameters for responders to treatments A, B, C
discount_y <- c(1, 1, 1)  # 0.2*1=0.2; 0.3*1=0.3; 0.4*1=0.4  
# The pre-specified values of linkage parameters for non-responders 
# who receive first second-stage treatment
discount_n1 <- c(1, 1, 1)   # 0.3*1=0.3 (A->B); 0.2*1=0.2 (B->A); 0.2*1=0.2 (C->A)  
discount_n2 <- c(1, 1, 1)   # 0.4*1=0.4 (A->C); 0.4*1=0.4 (B->C); 0.3*1=0.3 (C->B)

### Scenario 2
pi_1A <- 0.2        
pi_1B <- 0.3          
pi_1C <- 0.4    
discount_y <- c(2, 2, 2)   # 0.2*2=0.4; 0.3*2=0.6; 0.4*2=0.8  
discount_n1 <- c(1, 1, 1)  # 0.3*1=0.3 (A->B); 0.2*1=0.2 (B->A); 0.2*1=0.2 (C->A)   
discount_n2 <- c(1, 1, 1)  # 0.4*1=0.4 (A->C); 0.4*1=0.4 (B->C); 0.3*1=0.3 (C->B)

### Scenario 3
pi_1A <- 0.2        
pi_1B <- 0.3          
pi_1C <- 0.4    
discount_y <- c(1, 1, 1)  # 0.2*2=0.4; 0.3*2=0.6; 0.4*2=0.8  
discount_n1 <- c(1/2, 1/2, 1/2)  # 0.3*1/2=0.15 (A->B); 0.2*1/2=0.1 (B->A); 0.2*1/2=0.1 (C->A)   
discount_n2 <- c(1/2, 1/2, 1/2)  # 0.4*1/2=0.2 (A->C); 0.4*1/2=0.2 (B->C); 0.3*1/2=0.15 (C->B)

### Scenario 4
pi_1A <- 0.2        
pi_1B <- 0.3          
pi_1C <- 0.4    
discount_y <- c(2, 2, 2)  # 0.2*2=0.4; 0.3*2=0.6; 0.4*2=0.8  
discount_n1 <- c(3/2, 3/2, 3/2)  # 0.3*3/2=0.45 (A->B); 0.2*3/2=0.3 (B->A); 0.2*3/2=0.3 (C->A)  
discount_n2 <- c(3/2, 3/2, 3/2)  # 0.4*3/2=0.6 (A->C); 0.4*3/2=0.6 (B->C); 0.3*3/2=0.45 (C->B)

### Scenario 5
pi_1A <- 0.2        
pi_1B <- 0.3          
pi_1C <- 0.4    
discount_y <- c(3, 2, 3/2)  # 0.2*3=0.6; 0.3*2=0.6; 0.4*3/2=0.6  
discount_n1 <- c(2, 2, 2)  # 0.3*2=0.6 (A->B); 0.2*2=0.4 (B->A); 0.2*2=0.4 (C->A)  
discount_n2 <- c(1/2, 1/2, 1/2)  # 0.4*1/2=0.2 (A->C); 0.4*1/2=0.2 (B->C); 0.3*1/2=0.15 (C->B) 

### Scenario 6
pi_1A <- 0.3        
pi_1B <- 0.3          
pi_1C <- 0.3    
discount_y <- c(1, 1, 1)  # 0.3*1=0.3; 0.3*1=0.3; 0.3*1=0.3  
discount_n1 <- c(1, 1, 1)   # 0.3*1=0.3 (A->B); 0.3*1=0.3 (B->A); 0.3*1=0.3 (C->A)  
discount_n2 <- c(1, 1, 1)   # 0.3*1=0.3 (A->C); 0.3*1=0.3 (B->C); 0.3*1=0.3 (C->B)

### Scenario 7
pi_1A <- 0.3        
pi_1B <- 0.3          
pi_1C <- 0.3    
discount_y <- c(2/3, 1, 3/4)  # 0.3*2/3=0.2; 0.3*1=0.3; 0.3*4/3=0.4  
discount_n1 <- c(1, 2/3, 2/3)   # 0.3*1=0.3 (A->B); 0.3*2/3=0.2 (B->A); 0.3*2/3=0.2 (C->A)  
discount_n2 <- c(4/3, 4/3, 1)   # 0.3*4/3=0.4 (A->C); 0.3*4/3=0.4 (B->C); 0.3*1=0.3 (C->B)



# Number of iterations for MCMC chain
niter <- 50000

# Sample sizes for each arm in stage 1
n_A <- n_B <- n_C <- 30

# Number of simulation runs
# In the paper, the number is 10000
n_sim <- 5


#################################################################
##           Section 2.1: Bayesian Joint Stage Model           ##
#################################################################

# Mean of pi of treatment = 1/(1+1) = 0.5
pi_prior.a <- c(1,1,1)
pi_prior.b <- c(1,1,1)

# mean of beta0 = r/mu = 2/2 = 1, beta0 ~ gamma(r,mu)
beta0_prior.r <- 2
beta0_prior.mu <- 2     
#mean of beta1 = r/mu = 2/2 = 1, beta1 ~ gamma(r,mu)
beta1_prior.r <- 2           
beta1_prior.mu <- 2
# Number of arms in SMART
NUM_ARMS <- 3
# Number of MCMC samples in BJSM
MCMC_SAMPLE <- 6000
# Number of burn-in iterations in BJSM
BURN.IN <- 1000
# Number of MCMC chains
n_MCMC_chain <- 1

#################################################################
##           Section 2.2: Modified Power Prior Model           ##
#################################################################

# Prior hyperparameters of delta (weight) (prior mean of delta = 0.5)
prior_delta_a <- 1
prior_delta_b <- 1
# # (prior mean of delta = 0.2)
# prior_delta_a <- 0.4
# prior_delta_b <- 1.6
# # (prior mean of delta = 0.8)
# prior_delta_a <- 1.6
# prior_delta_b <- 0.4

# Prior hyperparameters of pi (response rate)
prior_pi_a <- 1
prior_pi_b <- 1

###########################################################################
###########################################################################
###                                                                     ###
###               SECTION 3: MODEL FITTING AND ESTIMATION               ###
###                                                                     ###
###########################################################################
###########################################################################


trt_result_MPP <- tibble()
trt_result_EB <- tibble()
trt_result_pen <- tibble()
trt_result_d1 <- tibble()
trt_result_d0 <- tibble()
trt_result_BJSM <- tibble()
trt_result_FET <- tibble()
trt_result_BM <- tibble()
pi_hat_bjsm <- c()

for (j in 1:n_sim) {
  mydata <- standard_SMART(pi_1A, pi_1B, pi_1C, discount_y, 
                           discount_n1, discount_n2, n_A, n_B, n_C)
  dataset <- mydata[[1]]
  # number of first stage responders to each treatment
  s1_n_resp <- mydata[[2]]
  # number of subjects to second-stage treatments following the first stage treatments
  s2_n_group <- mydata[[3]] 
  # number of second stage responders 
  s2_n_resp <- mydata[[4]]   
  
  #################################################################
  ##                      Section 3.1: BJSM                      ##
  #################################################################
  
  jag.model.name <- "BJSM_2beta_gamma_gamma.bug"  
  tryCatch({
    jag <- jags.model(file.path(jags_path,jag.model.name),
                      data=list(n = nrow(dataset),
                                num_arms = NUM_ARMS,
                                Y1 = dataset$response_stageI,
                                Y2 = dataset$response_stageII,
                                treatment_stageI = dataset$treatment_stageI,
                                treatment_stageII = dataset$treatment_stageII,
                                response1 = dataset$response_stageI + 1,
                                
                                #prior
                                pi_prior.a = pi_prior.a,
                                pi_prior.b = pi_prior.b,
                                beta0_prior.a = beta0_prior.r,    # gamma
                                beta0_prior.b = beta0_prior.mu,    # gamma
                                beta1_prior.a = beta1_prior.r,    # gamma
                                beta1_prior.c = beta1_prior.mu     # gamma
                      ),
                      n.chains=n_MCMC_chain,n.adapt = BURN.IN)   
    posterior_sample <- coda.samples(jag,
                                     c('pi','beta'),
                                     MCMC_SAMPLE)
  },
  warning = function(war){},
  error = function(err){},
  finally = {
    print(j)     # show the number of iterations run 
  }
  )
  out_post <- posterior_sample[[1]]
  pi_hat_bjsm <- rbind(pi_hat_bjsm, apply(out_post[,3:5], 2, mean))
  
  trt_result_BJSM <- 
    summary(out_post)[[1]] %>% 
    as_tibble(rownames = "type") %>% 
    dplyr::filter(grepl("pi", type)) %>% 
    mutate(type = c("pi_A_BJSM",
                    "pi_B_BJSM",
                    "pi_C_BJSM")) %>% 
    dplyr::select(type,
                  Mean,
                  SD) %>% 
    bind_rows(trt_result_BJSM, .)
  
  #################################################################
  ##              Section 3.2: MPP model estimation              ##
  #################################################################
  
  # The indices of the subgroups of patients involved in estimation of each trt effect
  subgroup_A <- c(1, 4, 6)
  subgroup_B <- c(3, 2, 6)
  subgroup_C <- c(5, 2, 4)
  
  try({
    init_param <- rep(0.1, 5)
    result <- MCMCmetrop1R(complete_likelihood_modified_power_prior_2_delta, theta.init = init_param,
                           burnin = 1 * niter, mcmc = 1 * niter, verbose = 0,
                           thin = 1 ,V = 0.0008 * diag(length(init_param)),
                           s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                           s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                           prior_delta_a = prior_delta_a, prior_delta_b = prior_delta_b, 
                           prior_pi_a = prior_pi_a,prior_pi_b = prior_pi_b)
    
    result <- MCMCmetrop1R(complete_likelihood_modified_power_prior_2_delta, theta.init = init_param,
                           burnin = 1 * niter, mcmc = 0.8 * niter, verbose = 0,
                           thin = 1 ,V = 0.1 * cov(result[trunc(0.8 * niter) : trunc(1 * niter),]),
                           s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                           s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                           prior_delta_a = prior_delta_a, prior_delta_b = prior_delta_b, 
                           prior_pi_a = prior_pi_a,prior_pi_b = prior_pi_b)
    
    result <- MCMCmetrop1R(complete_likelihood_modified_power_prior_2_delta, theta.init = init_param,
                           burnin = 1 * niter, mcmc = 0.8 * niter, verbose = 0,
                           thin = 1 ,V = 0.1 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                           s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                           s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                           prior_delta_a = prior_delta_a, prior_delta_b = prior_delta_b, 
                           prior_pi_a = prior_pi_a,prior_pi_b = prior_pi_b)
    
    result <- MCMCmetrop1R(complete_likelihood_modified_power_prior_2_delta, theta.init = init_param,
                           burnin = 1 * niter, mcmc = 0.8 * niter, verbose = 0,
                           thin = 5 ,V = 0.2 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                           s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                           s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                           prior_delta_a = prior_delta_a, prior_delta_b = prior_delta_b, 
                           prior_pi_a = prior_pi_a,prior_pi_b = prior_pi_b)
    
    
    result <- as.mcmc(as.data.frame((result)))
    varnames(result)[1:3] <- c("pi_A_MPP", "pi_B_MPP", "pi_C_MPP")
    for (i in 1:2) {
      varnames(result)[i + 3] <- paste0("delta_", i, "_MPP")
    }
    trt_result_MPP <- summary(result)[[1]] %>%
      as_tibble(rownames = "type") %>%
      dplyr::select(type, Mean, SD) %>%
      bind_rows(trt_result_MPP, .)
  })
  
  ##################################################################
  ##      Section 3.3: PP model estimation using MLC approach     ##
  ##################################################################
  
  
  
  # Estimate pi using the estimated delta
  
  # Initial value of delta for MLC
  delta_init <- rep(0.1, 2)

  
  test2 <- optimx(
    par = delta_init, 
    fn = complete_marginal_log_likelihood_2_delta, 
    gr = NULL,
    s2_n_group = s2_n_group, 
    s2_n_resp = s2_n_resp,
    s1_n_group = c(n_A, n_B, n_C), 
    s1_n_resp = s1_n_resp,
    prior_pi_a = prior_pi_a, 
    prior_pi_b = prior_pi_b,
    method = "Rvmmin", 
    lower = rep(0, length(delta_init)),
    upper = rep(1, length(delta_init)),
    control = list(maxit = 1000, trace = 1)
  )
  
  delta_hat <- c(test2$p1, test2$p2)
  
  # result table of delta hat
  delta_result <- 
    tibble(type = paste0("delta_", seq(1, 2, 1), "_EB"),
           Mean = delta_hat)
  
  init_param <- rep(0.1, 3)
  
  # Estimate pi using the estimated delta
  try({
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = init_param, burnin = 1 * niter, mcmc = 1 * niter, verbose = 0,
                   thin = 1, V = 0.0008 * diag(length(init_param)), deltas = delta_hat,
                   s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.2 * cov(result[trunc(0.8 * niter):trunc(1 * niter), ]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.4 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result),], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 5, 
                   V = 0.3 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- as.mcmc(as.data.frame((result))) 
    varnames(result)[1:3] <- c("pi_A_EB", "pi_B_EB", "pi_C_EB")
    
    trt_result_EB <- 
      summary(result)[[1]] %>% 
      as_tibble(rownames = "type") %>% 
      bind_rows(delta_result) %>% 
      dplyr::select(type,
                    Mean,
                    SD) %>% 
      bind_rows(trt_result_EB, .)
    
    
  })
  
  #############################################################################
  ##  Section 3.4: PP model estimation using penalized likelihood criterion  ##
  #############################################################################
  
  try({
    # Initial value of delta using penalized-likelihood criterion
    delta_init <- rep(0.1, 2)
    
    test2 <- optimx(
      par = delta_init,
      fn = complete_penalized_marginal_log_likelihood_2_delta,
      gr = NULL,
      s2_n_group = s2_n_group,
      s2_n_resp = s2_n_resp,
      s1_n_group = c(n_A, n_B, n_C),
      s1_n_resp = s1_n_resp,
      prior_pi_a = prior_pi_a,
      prior_pi_b = prior_pi_b,
      method = "Rvmmin",
      lower = rep(0, length(delta_init)),
      upper = rep(1, length(delta_init)),
      control = list(maxit = 1000, trace = 1)
    )
    
    delta_hat <- c(test2$p1, test2$p2)
    
    # result table of delta hat
    delta_result <- 
      tibble(type = paste0("delta_", seq(1, 2, 1), "_pen"),
             Mean = delta_hat)
    
    init_param <- rep(0.1, 3)
    # Estimate pi using the estimated delta
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = init_param, burnin = 1 * niter, mcmc = 1 * niter, verbose = 0,
                   thin = 1, V = 0.0008 * diag(length(init_param)), deltas = delta_hat,
                   s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.2 * cov(result[trunc(0.8 * niter):trunc(1 * niter), ]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.4 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result),], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 5, 
                   V = 0.3 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- as.mcmc(as.data.frame((result))) 
    varnames(result)[1:3] <- c("pi_A_pen", "pi_B_pen", "pi_C_pen")
    
    trt_result_pen <- 
      summary(result)[[1]] %>% 
      as_tibble(rownames = "type") %>% 
      bind_rows(delta_result) %>% 
      dplyr::select(type,
                    Mean,
                    SD) %>% 
      bind_rows(trt_result_pen, .)
    
  })
  
  
  ##################################################################
  ##             Section 3.4: PP model with delta = 1             ##
  ##################################################################
  delta_hat <- rep(1, 2)
  
  # result table of delta hat
  delta_result <- 
    tibble(type = paste0("delta_", seq(1, 2, 1),"_d1"),
           Mean = delta_hat)
  
  init_param <- rep(0.3, 3)
  
  # Estimate pi using the estimated delta
  try({
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = init_param, burnin = 1 * niter, mcmc = 1 * niter, verbose = 0,
                   thin = 1, V = 0.0008 * diag(length(init_param)), deltas = delta_hat,
                   s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.2 * cov(result[trunc(0.8 * niter):trunc(1 * niter), ]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.4 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result),], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 5, 
                   V = 0.3 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- as.mcmc(as.data.frame((result))) 
    varnames(result)[1:3] <- c("pi_A_d1", "pi_B_d1", "pi_C_d1")
    
    trt_result_d1 <- 
      summary(result)[[1]] %>% 
      as_tibble(rownames = "type") %>% 
      bind_rows(delta_result) %>% 
      dplyr::select(type,
                    Mean,
                    SD) %>% 
      bind_rows(trt_result_d1, .)
  })
  
  ##################################################################
  ##             Section 3.6: PP model with delta = 0             ##
  ##################################################################
  delta_hat <- rep(0, 2)
  
  # result table of delta hat
  delta_result <- 
    tibble(type = paste0("delta_", seq(1, 2, 1), "_d0"),
           Mean = delta_hat)
  
  init_param <- rep(0.3, 3)
  
  # Estimate pi using the estimated delta
  try({
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = init_param, burnin = 1 * niter, mcmc = 1 * niter, verbose = 0,
                   thin = 1, V = 0.0008 * diag(length(init_param)), deltas = delta_hat,
                   s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.2 * cov(result[trunc(0.8 * niter):trunc(1 * niter), ]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.4 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result),], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 5, 
                   V = 0.3 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- as.mcmc(as.data.frame((result))) 
    varnames(result)[1:3] <- c("pi_A_d0", "pi_B_d0", "pi_C_d0")
    
    trt_result_d0 <- 
      summary(result)[[1]] %>% 
      as_tibble(rownames = "type") %>% 
      bind_rows(delta_result) %>% 
      dplyr::select(type,
                    Mean,
                    SD) %>% 
      bind_rows(trt_result_d0, .)
  })
  
  ##################################################################
  ##                Section 3.7: Fisher Exact Test                ##
  ##################################################################
  
  ##### FET with 2 deltas
  
  delta_hat <- complete_FET_2_delta(s2_n_group, s2_n_resp, c(n_A, n_B, n_C), s1_n_resp)
  
  # result table of delta hat
  delta_result <- 
    tibble(type = paste0("delta_", seq(1, 2, 1), "_FET"),
           Mean = delta_hat)
  
  init_param <- rep(0.3, 3)
  
  # Estimate pi using the estimated delta
  try({
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = init_param, burnin = 1 * niter, mcmc = 1 * niter, verbose = 0,
                   thin = 1, V = 0.0008 * diag(length(init_param)), deltas = delta_hat,
                   s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.2 * cov(result[trunc(0.8 * niter):trunc(1 * niter), ]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.4 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result),], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 5, 
                   V = 0.3 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- as.mcmc(as.data.frame((result))) 
    varnames(result)[1:3] <- c("pi_A_FET", "pi_B_FET", "pi_C_FET")
    
    trt_result_FET <- 
      summary(result)[[1]] %>% 
      as_tibble(rownames = "type") %>% 
      bind_rows(delta_result) %>% 
      dplyr::select(type,
                    Mean,
                    SD) %>% 
      bind_rows(trt_result_FET, .)
    
    
  })
  
  ##################################################################
  ##              Section 3.8: Bhattacharyya Measure              ##
  ##################################################################
  
  ##### BM with 2 deltas
  
  delta_hat <- complete_BM_2_delta(s2_n_group, s2_n_resp, c(n_A, n_B, n_C), s1_n_resp, prior_pi_a, prior_pi_b)
  
  # result table of delta hat
  delta_result <- 
    tibble(type = paste0("delta_", seq(1, 2, 1), "_BM"),
           Mean = delta_hat)
  
  init_param <- rep(0.3, 3)
  
  # Estimate pi using the estimated delta
  try({
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = init_param, burnin = 1 * niter, mcmc = 1 * niter, verbose = 0,
                   thin = 1, V = 0.0008 * diag(length(init_param)), deltas = delta_hat,
                   s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta,
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.2 * cov(result[trunc(0.8 * niter):trunc(1 * niter), ]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result), ], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 1, 
                   V = 0.4 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- 
      MCMCmetrop1R(complete_likelihood_power_prior_2_delta, 
                   theta.init = result[nrow(result),], burnin = 1 * niter, mcmc = 0.8 * niter, 
                   verbose = 0, thin = 5, 
                   V = 0.3 * cov(result[trunc(0.6 * niter) : trunc(0.8 * niter),]),
                   deltas = delta_hat, s2_n_group = s2_n_group, s2_n_resp = s2_n_resp,
                   s1_n_group = c(n_A, n_B, n_C), s1_n_resp = s1_n_resp, 
                   prior_pi_a = prior_pi_a, prior_pi_b = prior_pi_b)
    
    result <- as.mcmc(as.data.frame((result))) 
    varnames(result)[1:3] <- c("pi_A_BM", "pi_B_BM", "pi_C_BM")
    
    trt_result_BM <- 
      summary(result)[[1]] %>% 
      as_tibble(rownames = "type") %>% 
      bind_rows(delta_result) %>% 
      dplyr::select(type,
                    Mean,
                    SD) %>% 
      bind_rows(trt_result_BM, .)
    
    
  })
}

result_table <- bind_rows(trt_result_MPP,
                         trt_result_EB,
                         trt_result_pen,
                         trt_result_d1,
                         trt_result_d0,
                         trt_result_BJSM,
                         trt_result_FET,
                         trt_result_BM)

############################################################################
############################################################################
###                                                                      ###
###                  SECTION 4: PRESENTATION OF RESULTS                  ###
###                                                                      ###
############################################################################
############################################################################

##########################################################################
##  Section 4.1: Estimation of pi and delta (Tables 2, 3 and Figure 4)  ##
##########################################################################

## The estimated delta and its standard deviation can be seen in the following table for the specific scenario
final_table <- result_table %>% 
  group_by(type) %>% 
  dplyr::summarize(estimate = mean(Mean),
                   emp_sd = sd(Mean)) %>%
  mutate(true = c(rep(NA, 14), rep(c(pi_1A, pi_1B, pi_1C), each = 8)),
         bias = estimate - true,
         rMSE = sqrt(bias^2 + emp_sd^2)) 


## To plot Figure 4, need to generate "final_table" using scenarios 4-7
## Save only the bias and rMSE for each scenario as follows

bias_rMSE <- bias_rMSE_organizer(final_table)

## The four tables below are the stored bias and rMSE information from scenarios 4-7, respectively

# bias_rMSE4 <- bias_rMSE
# bias_rMSE5 <- bias_rMSE
# bias_rMSE6 <- bias_rMSE
# bias_rMSE7 <- bias_rMSE
# bias_rMSE_all <- bind_rows(
#   bias_rMSE4,
#   bias_rMSE5,
#   bias_rMSE6,
#   bias_rMSE7
# )

## Uncomment the following code (lines 1098-1144) to plot Figure 4

# bias_barplot <- bias_rMSE_all %>% 
#   dplyr::select(matches("bias")) %>% 
#   rownames_to_column(var = "scenario") %>% 
#   pivot_longer(cols = starts_with("bias"),
#                names_to = "method", 
#                values_to = "bias")
# 
# rMSE_barplot <- bias_rMSE_all %>% 
#   dplyr::select(matches("rMSE")) %>% 
#   rownames_to_column(var = "scenario") %>% 
#   pivot_longer(cols = starts_with("rMSE"),
#                names_to = "method", 
#                values_to = "rMSE")
# 
# bias_barplot$scenario <- factor(bias_barplot$scenario, levels = c("4", "5", "6", "7"), 
#                                 labels = c("Scenario 4", "Scenario 5", "Scenario 6", "Scenario 7"))
# 
# rMSE_barplot$scenario <- factor(rMSE_barplot$scenario, levels = c("4", "5", "6", "7"), 
#                                 labels = c("Scenario 4", "Scenario 5", "Scenario 6", "Scenario 7"))
# 
# bias_barplot$method <- factor(bias_barplot$method, levels = c("bias_MPP", "bias_pen", "bias_EB", "bias_BM",
#                                                               "bias_FET", "bias_BJSM", "bias_d1", "bias_d0"), 
#                               labels = c("MPP", "PLC", "MLC", "BOM", 
#                                          "FET", "BJSM", expression(delta=1), expression(delta=0)))
# 
# rMSE_barplot$method <- factor(rMSE_barplot$method, levels = c("rMSE_MPP", "rMSE_pen", "rMSE_EB", "rMSE_BM",
#                                                               "rMSE_FET", "rMSE_BJSM", "rMSE_d1", "rMSE_d0"), 
#                               labels = c("MPP", "PLC", "MLC", "BOM", 
#                                          "FET", "BJSM", expression(delta=1), expression(delta=0)))
# 
# bias_4567 <- ggplot(bias_barplot %>% filter(scenario %in% c("Scenario 4", "Scenario 5", "Scenario 6", "Scenario 7")), 
#        aes(method, abs(bias))) + 
#   geom_bar(aes(method), stat="identity", position="dodge") +
#   facet_grid(. ~ scenario) + 
#   scale_x_discrete(labels = c("MPP", "PLC", "MLC", "BOM", 
#                               "FET", "BJSM", expression(delta~"=1"), expression(delta~"=0"))) + 
#   labs(y = "mean absolute bias of three treatment estimates", title = "bias") + 
#   theme(axis.text.x = element_text(size = 7))
# 
# rMSE_4567 <- ggplot(rMSE_barplot %>% filter(scenario %in% c("Scenario 4", "Scenario 5", "Scenario 6", "Scenario 7")), 
#                     aes(method, rMSE)) + 
#   geom_bar(aes(method), stat="identity", position="dodge") +
#   facet_grid(. ~ scenario) + 
#   scale_x_discrete(labels = c("MPP", "PLC", "MLC", "BOM", 
#                               "FET", "BJSM", expression(delta~"=1"), expression(delta~"=0"))) + 
#   labs(y = "mean rMSE of three treatment estimates", title = "rMSE") + 
#   theme(axis.text.x = element_text(size = 7))

#################################################################
##      Section 4.2: Histogram of delta (Figures 2 and 3)      ##
#################################################################

hist_table <- result_table %>% 
  dplyr::filter(grepl("delta", type)) %>% 
  dplyr::filter(!grepl("d0|d1", type)) %>% 
  separate(type, into = c("delta", "subgroup", "method"), sep = "_") %>% 
  dplyr::select(subgroup, method, Mean)

hist_table$method <- factor(hist_table$method, levels = c("MPP", "pen", "EB", "BM", "FET"), 
                        labels = c("MPP", "PLC", "MLC", "BOM", "FET"))

levels(hist_table$method)[levels(hist_table$method) == "BM"] <- "BOM"

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))
hist_delta <- ggplot(hist_table, aes(Mean, fill = subgroup)) + 
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.5) + 
  facet_grid(. ~ method) +
  scale_fill_manual(values = c("blue", "red"), labels = c(expression(delta[1]), expression(delta[2]))) + 
  ggtitle("Scenario X") + 
  labs(x = expression(delta))
hist_delta