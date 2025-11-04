library(mvnfast)
require(ezECM)

iters <- 250
Bcat_accuracy <- Bb_accuracy <- accuracies <- accuracy_detection <- rep(NA, times = iters)



## Global Arguments
alpha <- 0.95
alphatilde <- 0.05
vic <- "EX"

Cdetect <- matrix(1, ncol = 2, nrow = 2) - diag(2) # order of action columns: EX, not EX
Ccat <- matrix(1, ncol = 3, nrow = 3) - diag(3) # order of action columns:  c("EX","DEQ", "SEQ")

### Loading Data

Y <- read.csv(system.file("extdata", "ECM_validation_data_scrubbed.csv",
                          package="ezECM"))[,c("Source.Type", "pLP", "pCF_EKA", "pCF_GBA", "pCF_WRA", "pCF_YKA")]

all_missing <- apply(Y[,-1], 1, function(X){all(is.na(X))})

Y <- Y[-which(all_missing),]
names(Y)[1] <- "event"

transform <- "logit"#"arcsin"

Y_all <- Y

#iters <- 250 <- this is what the original MC experiment used
BT <- c(500, 50500)
thinning <- 5

########################
###### Begin Function Arguments
########################

for(ii in 1:iters){

Y <- Y_all


########### ########### ########### ########### ###########
########### Training and Testing Data Split
########### ########### ########### ########### ###########

Ytilde_ind <- sample(1:(nrow(Y)), size = 100)
Ytilde <- Y[Ytilde_ind,]

Y <- Y[-Ytilde_ind,]

Ytilde_truth <- Ytilde$event
Ytilde$event <- NULL

### naming vars as previous B-ECM model testing
train_full <- Y
testing <- Ytilde
test_truth <- Ytilde_truth

########################
###### End Function Arguments
########################

Ylist <- list()
cats <- unique(Y$event)

if(length(cats) <= 1){
  stop("Number of event categories must be 2 or greater.")
}

discriminants <- names(Y)[-which(names(Y) == "event")]

for(k in cats){
  Ylist[[k]] <- unname(as.matrix(Y[Y$event == k, -which(names(Y) == "event")]))
  if(transform == "logit"){
    Ylist[[k]] <- log(Ylist[[k]]) - log(1-Ylist[[k]])
  }else if(transform == "arcsin"){
    Ylist[[k]] <- (2/pi) * asin(sqrt(Ylist[[k]]))
  }else{
    stop("Set the variable transform <- 'logit', or transform <- 'asin'")
  }
}

cats <- sort(cats)

if(transform == "logit"){
  Ytilde <- log(Ytilde) - log(1-Ytilde)
  #Psi <- 0.5
  eta <- 0
}else{
  Ytilde <- (2/pi) * asin(sqrt(Ytilde))
  Psi <- 0.1# * diag(p)
  eta <- 0.5
}



Ylist <- Ylist[cats]

## Record indices of data for specific missingness patterns

missing_list <- list()
N_list <- list()

for(k in cats){

  N_list[[k]] <- list()

  na_to_bi <- apply(Ylist[[k]], 1, function(X){paste0(as.numeric(!is.na(X)), collapse = "")})
  u_missing <- unique(na_to_bi)
  #u_missing <- paste("m", u_missing, sep = "") ## doesn't work for indexing below, but may make sense to do this later
  missing_list[[k]] <- list()

  for(m in u_missing){

    missing_list[[k]][[m]] <- which(na_to_bi == m)

    Nmk <- length(missing_list[[k]][[m]])
    N_list[[k]][[m]] <- Nmk
  }
}

####################################################################################################################################################################################
########### Prior allocation
########### Prior mean will be the empirical mean over all obeservations
########### Prior covariance will be diagonal with all elements having the same value,
########### making implementation for now more simple compared to differing values or full covariance
########### Prior covariance elements will be fairly large relative to the [0,1] space containing the data
########### These choices are with the intent that inference will still be ok when no observations were made for a specific missingness patterns (global mean for k),
########### but low density evaluation for such a case, and for missingness patterns with observations, there is a smaller shift away from a global mean.
########### Further model development could introduce this global vs. local relationship in a more statistically consistent manner
####################################################################################################################################################################################

p <- ncol(Y) - 1
K <- length(cats)
N <- nrow(Y)

emp_mean <- data.frame(matrix(NA, nrow = p, ncol = length(cats)))
names(emp_mean) <- cats

#Psi <- 0.1# * diag(p)
#eta <- 0.5
atilde <- 0.5 # also only considering jeffreys prior over all combinations of missingness and events
nu <- 5

Psi <- list()

for(k in cats){

  emp_mean_temp <- apply(Ylist[[k]],2, mean, na.rm = TRUE)
  if(any(is.nan(emp_mean_temp))){
    emp_mean_temp[which(is.nan(emp_mean_temp))] <- eta
  }
  emp_mean[[k]] <- emp_mean_temp

  Psi[[k]] <- max(apply(Ylist[[k]],2,var,na.rm = TRUE), na.rm =TRUE)
}

# uncomment below to set the prior mean for all event categories to 0
#emp_mean[] <- 0

########### ########### ########### ########### ###########
########### Training Parameter Computation
########### ########### ########### ########### ###########

mu_tilde <- list()
Sig_tilde <- list()
nu_tilde <- list()
lpmztilde <- list()
lpmztilde_k <- list() ## missingness probabilities given \tilde{z}_k for use with typicality index
Nvec <-unname(sapply(Ylist, nrow))

missing_combos <- (2^p - 1) * K
missing_combos_k <- 2^p - 1
atilde_sum <- missing_combos * atilde
pmztilde_denom <- atilde_sum + N
lpmztilde_denom <- log(pmztilde_denom)

for(k in cats){

  mu_tilde[[k]] <- list()
  Sig_tilde[[k]] <- list()
  lpmztilde[[k]] <- list()

  kk <- which(cats %in% k)


  u_missing <- names(missing_list[[k]])
  for(m in u_missing){

    ptilde <- sum(as.numeric(unlist(strsplit(m, split = ""))))
    Nmk <- length(missing_list[[k]][[m]])

    lpmztilde[[k]][[m]] <- log(Nmk + atilde) - lpmztilde_denom
    lpmztilde_k[[k]][[m]] <- log(Nmk + atilde) - log(Nvec[kk] + missing_combos_k*atilde)


    nu_tilde[[k]][[m]] <- Nmk + 1  #Nmk + nu + 1 - ptilde

    Ymk <- Ylist[[k]][missing_list[[k]][[m]], ,drop = FALSE]

    ## Here special care will need to be taken if code allows different priors
    if(ptilde < p){
      Ymk <- Ymk[, -which(apply(Ymk,2,function(X){all(is.na(X))})), drop = FALSE]
    }

    ybar <- (1/Nmk)*t(Ymk ) %*% rep(1, times = Nmk) #apply(Ymk,2,mean, simplify = FALSE)

    etamk <- matrix(emp_mean[[k]][as.numeric(unlist(strsplit(m, split = ""))) == 1], ncol = ptilde, nrow = 1)

    mu_tilde[[k]][[m]] <- (Nmk * ybar + t(etamk))/(Nmk + 1)

    Eta <- rep(1, times = Nmk) %*% etamk
    IN1J <- diag(Nmk) - 1/(Nmk + 1) * matrix(1, ncol = Nmk, nrow = Nmk)
    Sig_tilde[[k]][[m]] <- ((Nmk + 2)/((Nmk + 1)*(Nmk + ptilde + 1 - ptilde))) * (Psi[[k]] * diag(ptilde) + t(Ymk - Eta) %*% IN1J %*% (Ymk - Eta))#((Nmk + 2)/((Nmk + 1)*(Nmk + nu + 1 - ptilde))) * (Psi * diag(ptilde) + t(Ymk - Eta) %*% IN1J %*% (Ymk - Eta))

  }

}


########### ########### ########### ########### ###########
########### Prediction
########### ########### ########### ########### ###########



lpztildek_numerator <- matrix(NA, ncol = K, nrow = nrow(Ytilde))

ltyp_ind <- matrix(NA, ncol = K, nrow = nrow(Ytilde))

for(i in 1:(nrow(Ytilde))){

  ytilde <- Ytilde[i, ]
  ytilde <- unlist(unname(ytilde))

  mtilde <- paste0(as.numeric(!is.na(ytilde)), collapse = "")
  ptilde <- sum(as.numeric(unlist(strsplit(mtilde, split = ""))))

  if(any(is.na(ytilde))){
    ytilde <- ytilde[-which(is.na(ytilde))]
  }

  ldens <- rep(NA, times = K)

  for(k in cats){

    kk <- which(cats %in% k)

    if(mtilde %in% names(missing_list[[k]])){


      mutilde_use <- mu_tilde[[k]][[mtilde]]
      Sigtilde_use <- Sig_tilde[[k]][[mtilde]]
      nutilde_use <- nu_tilde[[k]][[mtilde]]
      lpmztilde_use <- lpmztilde[[k]][[mtilde]]

      shift <- (ytilde - mutilde_use)
      btilde <- as.numeric((t(shift) %*% solve(Sigtilde_use, shift))/ptilde)

      mutilde_use <- as.vector(mutilde_use)

      ldens[kk] <- mvnfast::dmvt(ytilde, mu = mutilde_use, sigma = Sigtilde_use, df = nutilde_use, log = TRUE)

      ltyp_ind[i, kk] <- pf(btilde, ptilde, nutilde_use, log.p = TRUE) #+ lpmztilde_k[[k]][[mtilde]]


    }else{


      mutilde_use <- as.vector(emp_mean[[k]][as.numeric(unlist(strsplit(mtilde, split = ""))) == 1])
      Sigtilde_use <- 2 * (1/( 1)) * Psi[[k]] * diag(ptilde) #2 * (1/(nu + 1 - ptilde)) * Psi * diag(ptilde)
      nutilde_use <-  1 #(nu + 1 - ptilde)
      lpmztilde_use <- log(atilde) - lpmztilde_denom

      shift <- (ytilde - mutilde_use)
      btilde <- as.numeric((t(shift) %*% solve(Sigtilde_use, shift))/ptilde)

      mutilde_use <- as.vector(mutilde_use)

      ldens[kk] <- mvnfast::dmvt(ytilde, mu = mutilde_use, sigma = Sigtilde_use, df = nutilde_use, log = TRUE)

      ltyp_ind[i, kk] <- pf(btilde, ptilde, nutilde_use, log.p = TRUE) #+ log(atilde) - log(Nvec[kk] + missing_combos_k * atilde)


    }

    lpztildek_numerator[i,kk] <- ldens[kk] + lpmztilde_use


  }

}


#####################################################################
######### Typicality index, over all missingness patterns for each k
########################################################################

# typ_ind_temp <- rep(NA, times = 2^p -1)
# typ_ind <- matrix(NA, ncol = K, nrow = nrow(Ytilde))
#
# all_missing_patterns <- do.call(expand.grid, replicate(p, c(0, 1), simplify = FALSE))
# all_missing_patterns <- apply(all_missing_patterns, 1, function(X){paste(X, collapse = "")})
#
# lpmztilde_k_denom <- log(Nvec + missing_combos_k*atilde)
#
#
#
# for(i in 1:(nrow(Ytilde))){
#
#   ytilde <- Ytilde[i, ]
#   ytilde <- unlist(unname(ytilde))
#
#   mtilde <- paste0(as.numeric(!is.na(ytilde)), collapse = "")
#   ptilde <- sum(as.numeric(unlist(strsplit(mtilde, split = ""))))
#
#   if(any(is.na(ytilde))){
#     ytilde <- ytilde[-which(is.na(ytilde))]
#   }
#
#   for(k in cats){
#
#     for(m in all_missing_patterns){
#
#       mu_use <- emp_mean[[k]]
#
#
#       m_ind <- which(unlist(strsplit(m, split ="")) == "1")
#
#       mm <- which(all_missing_patterns %in% m)
#       ptilde_m <- sum(as.numeric(unlist(strsplit(m, split = ""))))
#
#       if(m %in% names(missing_list[[k]])){
#
#
#         mutilde_use <- as.vector(mu_tilde[[k]][[m]])
#         Sigtilde_use <- Sig_tilde[[k]][[m]]
#         nutilde_use <- nu_tilde[[k]][[m]]
#         lpmztilde_use <- lpmztilde_k[[k]][[m]]
#
#         shift <-
#
#
#
#         typ_ind_temp[mm] <- mvnfast::dmvt(ytilde, mu = mutilde_use, sigma = Sigtilde_use, df = nutilde_use, log = TRUE)
#
#       }else{
#
#
#
#         mutilde_use <- as.vector(emp_mean[[k]][as.numeric(unlist(strsplit(mtilde, split = ""))) == 1])
#         Sigtilde_use <- 2 * (1/( 1)) * Psi * diag(ptilde) #2 * (1/(nu + 1 - ptilde)) * Psi * diag(ptilde)
#         nutilde_use <-  1 #(nu + 1 - ptilde)
#         lpmztilde_use <- log(atilde) - lpmztilde_k_denom
#
#         typ_ind_temp[mm] <- mvnfast::dmvt(ytilde, mu = mutilde_use, sigma = Sigtilde_use, df = nutilde_use, log = TRUE)
#
#
#       }
#
#       lpztildek_numerator[i,kk] <- ldens[kk] + lpmztilde_use
#
#
#
#     }
#
#     }
#
#   }
#
# }



pztildek_numerator <- exp(lpztildek_numerator)
pztildek_denom <- pztildek_numerator %*% rep(1, times = K)
lpztildek_denom <- log(pztildek_denom)

lpztildek <-  (lpztildek_numerator - as.vector(lpztildek_denom))

pztildek <- exp(lpztildek)

min_exp_loss <- apply(pztildek, 1, which.max)

ztilde <- cats[min_exp_loss]

## Full categorization accuracy?
accuracies[ii] <- mean(ztilde == Ytilde_truth)

p_binary <- data.frame(pztildek %*% matrix(c(0,1,0,1,0,1), ncol = 2, nrow = 3))
names(p_binary) <- c("EX","EQ")

min_bin_loss <- apply(p_binary, 1, which.max)
ztilde_bin <- c("EX","EQ")[min_bin_loss]

accuracy_detect_i <- rep(NA, times = nrow(Ytilde))

typ_ind <- exp(ltyp_ind)

for(i in 1:(nrow(Ytilde))){

  ztrue <- Ytilde_truth[i]

  if(ztrue %in% c("SEQ", "DEQ")){
    if(ztilde_bin[i] == "EQ"){
      accuracy_detect_i[i] <- TRUE
    }else{
      accuracy_detect_i[i] <- FALSE
    }
  }else if(ztrue == "EX"){
    if(ztilde_bin[i] == "EX"){
      accuracy_detect_i[i] <- TRUE
    }else{
      accuracy_detect_i[i] <- FALSE
    }
  }

  }

accuracy_detection[ii] <- mean(accuracy_detect_i)


### Fitting the previous models on this subset of training and testing data

# cECM_fit <- cECM(x = train_full, transform = TRUE, newdata = testing)
#
# cECM_out <-  apply(cECM_decision(pval = cECM_fit, alphatilde = alphatilde,
#                                  vic = as.character(vic),
#                                  cat_truth = test_truth)$events[,1:3] ,2,
#                    sum, na.rm = TRUE)/100


bayes_fit <- BayesECM(Y = train_full, BT = BT, verb = TRUE)

bayes_pred <- predict(bayes_fit, Ytilde = testing)

p_detect <- matrix(bayes_pred$epz$EX, ncol = 2, nrow = length(bayes_pred$epz$EX))

p_detect[,2] <- 1- p_detect[,2]

# Column order:  EX, not EX
loss_detect <- p_detect %*% Cdetect

min_loss_bin <- c("EX", "EQ")[apply(loss_detect, 1, which.min)]
truth_binary <- rep("EQ", times = length(Ytilde_truth))
truth_binary[Ytilde_truth == "EX"] <- "EX"

Bb_accuracy[ii] <- mean(min_loss_bin == truth_binary)


p_cat <- bayes_pred$epz[c("EX","DEQ", "SEQ")]

# Column order: c("EX","DEQ", "SEQ")
loss_cat <- as.matrix(p_cat) %*% Ccat

min_loss_cat <- c("EX","DEQ", "SEQ")[apply(loss_cat,1,which.min)]

Bcat_accuracy[ii] <- mean(min_loss_cat == Ytilde_truth)

}



####################################################
############# Scratch work for coding up experiment
####################################################

# Important functions
cECM <- cECM(x = train_full, transform = TRUE, newdata = testing)

cECM_out <-  apply(cECM_decision(pval = cECM, alphatilde = alphatilde,
                                 vic = as.character(vic),
                                 cat_truth = test_truth)$events[,1:3] ,2,
                   sum, na.rm = TRUE)


bayes_fit <- BayesECM(Y = train_full)

bayes_pred <- predict(bayes_fit, Ytilde = testing,
                      mixture_weights = mixture_weights)















####################################################
############# Scratch work for coding up experiment
####################################################


### Current version of the plot
fctr <- 1

pdf(file = "missingness_bxplt.pdf", width = 8*fctr, height = 4.52*fctr)
group_labels <- c("Detection", "Event Categorization")   # Labels for the x-axis ticks
pltcols <- hcl.colors(44, palette = "viridis", rev = TRUE)[c(1.5,5,10, 15, 18.5)*2 ]
mcols <- hcl.colors(100, palette = "Zissou1")[c(82,70)]
box_colors   <- c(pltcols[3], mcols[1], pltcols[5], mcols[2])

df <- data.frame(detect = Bb_accuracy, m_detect = accuracy_detection, cat = Bcat_accuracy, m_cat =accuracies)

positions <- c(1,2,4, 5)  # space between 2 and 4 is the group gap

# Plot

par(mar = c(3,4,0.5,0.5), cex.axis = 1.1, cex.lab = 1.1)
boxplot(df,
        at = positions,                  # Position boxes with gap
        col = box_colors,                 # Colors from vector
        xaxt = "n",                        # Suppress default x-axis
        ylab = "Model Accuracy",
        main = "",
        ylim = c(0.72,1))

# Add group labels centered on groups
axis(1, at = c(mean(c(1, 2)), mean(c(4, 5))), labels = group_labels)

# Add legend
legend("topright",
       legend = c("M-B-ECM", "M-B-ECM and Missingness", "M-B-ECM Cat", "M-B-ECM Cat and Missingness"),
       fill = box_colors,
       bty = "n",
       ncol = 2, cex = 1.1)
       #horiz = TRUE)
dev.off()

######## Plot of Current Results

cols <- hcl.colors(10, palette = "Zissou1")[c(1,8)]

png("boxplot.png", width = 2800, height = 1200, res = 300)

par(mar = c(3,4,2,0.5))
boxplot(cbind(accuracy_detection, accuracies)*100, pars =list(boxwex =0.5, staplewex = 0.3), names = c("", ""), main = "B-ECM With Missingness Informing Categorization", col = cols, ylim = c(0.65,1)*100,
        ylab = "Accuracy (%)"#,
        #cex.main = 0.9,
        #cex.lab = 0.9,
        #cex.axis = 0.9
        )
axis(side = 1, at = 1:2, tick = FALSE, labels = c("Binary\nDetection", "Full Event\nCategorization"), line = 1, cex.axis = 0.9)
dev.off()



### Taking a look at some properties of the total of the data
Y <- read.csv(system.file("extdata", "ECM_validation_data_scrubbed.csv",
                          package="ezECM"))[,c("Source.Type", "pLP", "pCF_EKA", "pCF_GBA", "pCF_WRA", "pCF_YKA")]

all_missing <- apply(Y[,-1], 1, function(X){all(is.na(X))})

Y <- Y[-which(all_missing),]
names(Y)[1] <- "event"

transform <- "logit"#"arcsin"


Ytilde_ind <- sample(1:(nrow(Y)), size = 100)
Ytilde <- Y[Ytilde_ind,]

Y <- Y[-Ytilde_ind,]

Ytilde_truth <- Ytilde$event
Ytilde$event <- NULL


Ylist <- list()
cats <- unique(Y$event)

if(length(cats) <= 1){
  stop("Number of event categories must be 2 or greater.")
}

discriminants <- names(Y)[-which(names(Y) == "event")]

for(k in cats){
  Ylist[[k]] <- unname(as.matrix(Y[Y$event == k, -which(names(Y) == "event")]))
  if(transform == "logit"){
    Ylist[[k]] <- log(Ylist[[k]]) - log(1-Ylist[[k]])
  }else{
    Ylist[[k]] <- (2/pi) * asin(sqrt(Ylist[[k]]))
  }
}

lapply(Ylist, function(X){apply(X,2,function(x){mean(!is.na(x))})})

lapply(Ylist, function(X){apply(X,2,function(x){mean(x, na.rm = TRUE)})})

xgrid <- seq(from = -5, to = 5, length.out =250)

Xgrid <- matrix(xgrid, ncol = 1, nrow = length(xgrid))

exdens <- dmvt(Xgrid, mu = mu_tilde[["EX"]][["10000"]], sigma = Sig_tilde[["EX"]][["10000"]], df = nu_tilde[["EX"]][["10000"]])
seqdens <- dmvt(Xgrid, mu = mu_tilde[["SEQ"]][["10000"]], sigma = Sig_tilde[["SEQ"]][["10000"]], df = nu_tilde[["SEQ"]][["10000"]])

plot(xgrid, exdens, type = "l", col = "red", main = "Overlap of predictive densities for Explosions (red) and Shallow Earthquakes (blue)")
lines(xgrid, seqdens, col = "blue")

lapply(lpmztilde[["EX"]], exp)
lapply(lpmztilde[["SEQ"]], exp)
