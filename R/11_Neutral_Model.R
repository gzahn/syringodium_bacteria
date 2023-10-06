

# Define functions ####


## mle ####
#  File src/library/stats4/R/mle.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/
setClass("mle", representation(call = "language",
                               coef = "numeric",
                               fullcoef = "numeric",
                               vcov = "matrix",
                               min = "numeric",
                               details = "list",
                               minuslogl = "function",
                               nobs = "integer",
                               method = "character"))

setClass("summary.mle", representation(call = "language",
                                       coef = "matrix",
                                       m2logL = "numeric"))

setClass("profile.mle", representation(profile="list",
                                       summary="summary.mle"))

mle <- function(minuslogl, start = formals(minuslogl), method = "BFGS",
                fixed = list(), nobs, ...)
{
  # Insert sanity checks here...
  call <- match.call()
  n <- names(fixed)
  fullcoef <- formals(minuslogl)
  if(any(! n %in% names(fullcoef)))
    stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
  fullcoef[n] <- fixed
  if(!missing(start) && (!is.list(start) || is.null(names(start))))
    stop("'start' must be a named list")
  start[n] <- NULL
  start <- sapply(start, eval.parent) # expressions are allowed
  nm <- names(start)
  oo <- match(nm, names(fullcoef))
  if (anyNA(oo))
    stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
  start <- start[order(oo)]
  nm <- names(start) ## reordered names needed
  f <- function(p){
    l <- as.list(p)
    names(l) <- nm
    l[n] <- fixed
    do.call("minuslogl", l)
  }
  oout <- if (length(start))
    optim(start, f, method = method, hessian = TRUE, ...)
  else list(par = numeric(), value = f(start))
  coef <- oout$par
  vcov <- if(length(coef)) solve(oout$hessian) else matrix(numeric(), 0L, 0L)
  min <- oout$value
  fullcoef[nm] <- coef
  new("mle", call = call, coef = coef, fullcoef = unlist(fullcoef),
      vcov = vcov, min = min, details = oout, minuslogl = minuslogl,
      nobs = if(missing(nobs)) NA_integer_ else nobs,
      method = method)
}

setGeneric("coef")
setMethod("coef", "mle", function(object) object@fullcoef )
setMethod("coef", "summary.mle", function(object) object@coef )

setMethod("show", "mle", function(object){
  cat("\nCall:\n")
  print(object@call)
  cat("\nCoefficients:\n")
  print(coef(object))
})

setMethod("show", "summary.mle", function(object){
  cat("Maximum likelihood estimation\n\nCall:\n")
  print(object@call)
  cat("\nCoefficients:\n")
  print(coef(object))
  cat("\n-2 log L:", object@m2logL, "\n")
})

setGeneric("summary")
setMethod("summary", "mle", function(object, ...){
  cmat <- cbind(Estimate = object@coef,
                `Std. Error` = sqrt(diag(object@vcov)))
  m2logL <- 2*object@min
  new("summary.mle", call = object@call, coef = cmat, m2logL = m2logL)
})

setGeneric("profile")
setMethod("profile", "mle",
          function (fitted, which = 1L:p, maxsteps = 100,
                    alpha = 0.01, zmax = sqrt(qchisq(1 - alpha, 1L)),
                    del = zmax/5, trace = FALSE, ...)
          {
            onestep <- function(step)
            {
              bi <- B0[i] + sgn * step * del * std.err[i]
              fix <- list(bi)
              names(fix) <- pi
              call$fixed <- c(fix, fix0)
              pfit <- tryCatch(eval.parent(call, 2L), error = identity)
              if(inherits(pfit, "error")) return(NA)
              else {
                zz <- 2*(pfit@min - fitted@min)
                ri <- pv0
                ri[, names(pfit@coef)] <- pfit@coef
                ri[, pi] <- bi
                
                if (zz > -0.001)
                  zz <- max(zz, 0)
                else stop("profiling has found a better solution, so original fit had not converged")
                z <- sgn * sqrt(zz)
                pvi <<- rbind(pvi, ri)
                zi <<- c(zi, z)
              }
              if (trace) cat(bi, z, "\n")
              z
            }
            ## Profile the likelihood around its maximum
            ## Based on profile.glm in MASS
            summ <- summary(fitted)
            std.err <- summ@coef[, "Std. Error"]
            Pnames <- names(B0 <- fitted@coef)
            pv0 <- t(as.matrix(B0))
            p <- length(Pnames)
            prof <- vector("list", length = length(which))
            names(prof) <- Pnames[which]
            call <- fitted@call
            call$minuslogl <- fitted@minuslogl
            ndeps <- eval.parent(call$control$ndeps)
            parscale <- eval.parent(call$control$parscale)
            fix0 <- eval.parent(call$fixed)
            for (i in which) {
              zi <- 0
              pvi <- pv0
              pi <- Pnames[i]
              if (!is.null(ndeps)) call$control$ndeps <- ndeps[-i]
              if (!is.null(parscale)) call$control$parscale <- parscale[-i]
              for (sgn in c(-1, 1)) {
                if (trace)
                  cat("\nParameter:", pi, c("down", "up")[(sgn + 1)/2 + 1], "\n")
                step <- 0
                z <- 0
                
                ## This logic was a bit frail in some cases with
                ## high parameter curvature. We should probably at least
                ## do something about cases where the mle call fails
                ## because the parameter gets stepped outside the domain.
                ## (We now have.)
                
                call$start <- as.list(B0)
                lastz <- 0
                while ((step <- step + 1) < maxsteps && abs(z) < zmax) {
                  z <- onestep(step)
                  if(is.na(z)) break
                  lastz <- z
                }
                if(abs(lastz) < zmax) {
                  ## now let's try a bit harder if we came up short
                  for(dstep in c(0.2, 0.4, 0.6, 0.8, 0.9)) {
                    z <- onestep(step - 1 + dstep)
                    if(is.na(z) || abs(z) > zmax) break
                  }
                } else if(length(zi) < 5L) { # try smaller steps
                  mxstep <- step - 1L
                  step <- 0.5
                  while ((step <- step + 1L) < mxstep) onestep(step)
                }
              }
              si <- order(pvi[, i])
              prof[[pi]] <- data.frame(z = zi[si])
              prof[[pi]]$par.vals <- pvi[si,, drop=FALSE]
            }
            new("profile.mle", profile = prof, summary = summ)
          })

setGeneric("plot")
setMethod("plot", signature(x="profile.mle", y="missing"),
          function (x, levels, conf = c(99, 95, 90, 80, 50)/100, nseg = 50,
                    absVal = TRUE, ...)
          {
            ## Plot profiled likelihood
            ## Based on profile.nls (package stats)
            obj <- x@profile
            
            confstr <- NULL
            if (missing(levels)) {
              levels <- sqrt(qchisq(pmax(0, pmin(1, conf)), 1))
              confstr <- paste0(format(100 * conf), "%")
            }
            if (any(levels <= 0)) {
              levels <- levels[levels > 0]
              warning("levels truncated to positive values only")
            }
            if (is.null(confstr)) {
              confstr <- paste0(format(100 * pchisq(levels^2, 1)), "%")
            }
            mlev <- max(levels) * 1.05
            nm <- names(obj)
            opar <- par(mar = c(5, 4, 1, 1) + 0.1)
            if (absVal) {
              ## OBS: it is important to use the actual names for indexing,
              ## in case profile(....,which=....) was used
              for (i in nm) {
                ## <FIXME> This does not need to be monotonic
                sp <- splines::interpSpline(obj[[i]]$par.vals[, i], obj[[i]]$z,
                                            na.action=na.omit)
                bsp <- splines:: backSpline(sp)
                ## </FIXME>
                xlim <- predict(bsp, c(-mlev, mlev))$y
                if (is.na(xlim[1L]))
                  xlim[1L] <- min(obj[[i]]$par.vals[, i])
                if (is.na(xlim[2L]))
                  xlim[2L] <- max(obj[[i]]$par.vals[, i])
                dev.hold()
                plot(abs(z) ~ par.vals[, i], data = obj[[i]], xlab = i,
                     ylim = c(0, mlev), xlim = xlim, ylab = expression(abs(z)),
                     type = "n")
                avals <- rbind(as.data.frame(predict(sp)),
                               data.frame(x = obj[[i]]$par.vals[, i],
                                          y = obj[[i]]$z))
                avals$y <- abs(avals$y)
                lines(avals[order(avals$x), ], col = 4)
                abline(v = predict(bsp, 0)$y, h=0, col = 3, lty = 2)
                for (lev in levels) {
                  ## Note: predict may return NA if we didn't profile
                  ## far enough in either direction. That's OK for the
                  ## "h" part of the plot, but the horizontal line made
                  ## with "l" disappears.
                  pred <- predict(bsp, c(-lev, lev))$y
                  lines(pred, rep(lev, 2), type = "h", col = 6, lty = 2)
                  pred <- ifelse(is.na(pred), xlim, pred)
                  lines(pred, rep(lev, 2), type = "l", col = 6, lty = 2)
                }
                dev.flush()
              }
            }
            else {
              for (i in nm) {
                ## <FIXME> This does not need to be monotonic
                sp <- splines::interpSpline(obj[[i]]$par.vals[, i], obj[[i]]$z,
                                            na.action=na.omit)
                bsp <- splines::backSpline(sp)
                ## </FIXME>
                xlim <- predict(bsp, c(-mlev, mlev))$y
                x0 <- predict(bsp, 0)$y
                if (is.na(xlim[1L]))
                  xlim[1L] <- min(obj[[i]]$par.vals[, i])
                if (is.na(xlim[2L]))
                  xlim[2L] <- max(obj[[i]]$par.vals[, i])
                dev.hold()
                plot(z ~ par.vals[, i], data = obj[[i]], xlab = i,
                     ylim = c(-mlev, mlev), xlim = xlim, ylab = expression(z),
                     type = "n")
                lines(predict(sp), col = 4)
                abline(h = 0, v=x0, col = 3, lty = 2)
                for (lev in levels) {
                  pred <- predict(bsp, c(-lev, lev))$y
                  lines(pred, c(-lev, lev), type = "h", col = 6, lty = 2)
                  pred <- ifelse(is.na(pred), xlim, pred)
                  lines(c(x0,pred[2L]), rep(lev, 2), type = "l", col = 6, lty = 2)
                  lines(c(pred[1L],x0), rep(-lev, 2), type = "l", col = 6, lty = 2)
                }
                dev.flush()
              }
            }
            par(opar)
          })

setGeneric("confint")
setMethod("confint", "profile.mle",
          function (object, parm, level = 0.95, ...)
          {
            ## Calculate confidence intervals based on likelihood
            ## profiles
            of <- object@summary
            pnames <- rownames(of@coef)
            if (missing(parm))
              parm <- seq_along(pnames)
            if (is.character(parm))
              parm <- match(parm, pnames, nomatch = 0L)
            a <- (1 - level)/2
            a <- c(a, 1 - a)
            pct <- paste(round(100 * a, 1), "%")
            ci <- array(NA, dim = c(length(parm), 2L),
                        dimnames = list(pnames[parm], pct))
            cutoff <- qnorm(a)
            for (pm in parm) {
              pro <- object@profile[[pnames[pm]]]
              sp <- if (length(pnames) > 1L)
                spline(x = pro[, "par.vals"][, pm], y = pro[, 1L])
              else spline(x = pro[, "par.vals"], y = pro[, 1L])
              ci[pnames[pm], ] <- approx(sp$y, sp$x, xout = cutoff)$y
            }
            drop(ci)
          })

setMethod("confint", "mle",
          function (object, parm, level = 0.95, ...)
          {
            cat("Profiling...\n")
            confint(profile(object), alpha = (1 - level)/4, parm, level, ...)
          })

setGeneric("nobs")
setMethod("nobs", "mle", function (object, ...)
  if("nobs" %in% slotNames(object)) object@nobs else NA_integer_)

setGeneric("logLik")
setMethod("logLik", "mle", function(object, ...) {
  if(!missing(...))
    warning("extra arguments discarded")
  val <- -object@min
  if ("nobs" %in% slotNames(object) && # introduced in 2.13.0
      !is.na(no <- object@nobs)) attr(val, "nobs") <- no
  attr(val, "df") <- length(object@coef)
  class(val) <- "logLik"
  val
})

setGeneric("vcov")
setMethod("vcov", "mle", function (object, ...) object@vcov)

setGeneric("update")
setMethod("update", "mle", function (object, ..., evaluate = TRUE)
{
  call <- object@call
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) ) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate) eval(call, parent.frame()) else call
})

## sncm.fit ####

#Adam Burns - 2/10/2015
#aburns2@uoregon.edu
#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
#spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#pool: A community table for defining source community (optional; Default=NULL).
#taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.

spp <- spp.rare


sncm.fit <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  # require(stats4)
  
  options(warn=-1)
  
  #Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(pool)){
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  
  #Combine
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  
  #Calculate the limit of detection
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle2(sncm.LL, start=list(m=0.1, sigma=0.1),method="Nelder-Mead")
  # mle2()
  # ?mle()
  # ##Calculate Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, N, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle2(bino.LL, start=list(mu=0, sigma=0.1),method="Nelder-Mead")
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)
  
  ##Goodness of fit for Poisson model
  pois.pred <- ppois(d, N*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##Results
  if(stats==TRUE){
    fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    return(B)
  }
}




# Setup ####
library(tidyverse)
library(phyloseq)
library(bbmle)

set.seed(666)

#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
#spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#pool: A community table for defining source community (optional; Default=NULL).
#taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.


ps <- readRDS("./output/clean_phyloseq_object.RDS")

# Neutral model (Sloan, 2006) ####
# pull out components for neutral model
spp <- ps@otu_table %>% as("matrix") %>% as.data.frame()
spp.rare <- vegan::rrarefy(spp,min(rowSums(spp))) # rarefy to minumum depth

taxon <- ps@tax_table %>% as('matrix') %>% as.data.frame()


## Fit neutral model ####
sncm_stats <- sncm.fit(spp.rare,taxon=taxon)
sncm_fit <- sncm.fit(spp.rare,taxon=taxon,stats = FALSE) %>% 
  mutate(mean_realbund = spp.rare %>% apply(1,function(x){x/sum(x)}) %>% apply(1,mean),
         neutral_group = case_when(freq < pred.lwr ~ "Below",
                                   freq > pred.upr ~ "Above",
                                   freq >= pred.lwr & freq <= pred.upr ~ "Within"))



sncm_fit %>% 
  filter(mean_realbund > 0) %>% 
  ggplot(aes(y=freq,x=log(mean_realbund),color=neutral_group)) +
  geom_point()

outside_taxa <- 
row.names(sncm_fit)[which(sncm_fit$neutral_group != "Within")]

outside_taxa_df <- 
sncm_fit[outside_taxa,]
sig_taxa <- readRDS("./output/final_significant_taxa.RDS")

outside_taxa_df <- 
outside_taxa_df[outside_taxa_df$Genus %in% sig_taxa,]
outside_taxa_df <- outside_taxa_df[-1,]
outside_taxa_df <- outside_taxa_df[-6,]



outside_taxa_df %>% 
  ggplot(aes(y=Genus)) +
  geom_errorbarh(aes(xmin=pred.lwr,xmax=pred.upr),height=.2) +
  geom_point(aes(x=freq),color='red',size=3) +
  theme_minimal() +
  theme(axis.text.y = element_text(face='bold.italic')) +
  labs(x="Frequency")
ggsave("./output/figs/neutral_model_predictions_and_sig_taxa.png",dpi=300,width = 3,height = 3)



## Do repeated rarefaction ####

sncm.list <- list()
# for(i in 1:100){
#   spp.rare <- vegan::rrarefy(spp,min(rowSums(spp))) # rarefy to minumum depth
#   sncm <- sncm.fit(spp.rare,taxon=taxon)
#   sncm.list[[i]] <- sncm
# }


# Interpretation ####

## migration rate ####
# This estimated migration rate is the probability that a random loss 
# of an individual in a local community will be replaced by dispersal from 
# the metacommunity, as opposed to reproduction within the local community, 
# and can thus be interpreted as a measure of dispersal limitation.

migration_rate <- sncm.list %>% map_dbl("m")

### confidence intervals ####
# Binomial proportion 95% confidence intervals around the model predictions 
# were calculated using the Wilson score interval in the HMisc package in R

migration_rate_ci <- sncm.list %>% map_dbl("m.ci")

### set neutral boundaries ####
neutral_boundaries <- data.frame(mean = mean(migration_rate)) %>% 
  mutate(upper=mean + mean(migration_rate_ci),
         lower=mean - mean(migration_rate_ci))


# Divide OTU table ####

# find taxa that occur lower, within, and higher than neutral boundary

# divide up by side of WL
west <- ps %>% 
  subset_samples(east_west == "West")
east <- ps %>% 
  subset_samples(east_west == "East")



# run neutral model on samples (east vs west)
west_df <- west %>% otu_table() %>% as('matrix') %>% as.data.frame()
west_rare <- vegan::rrarefy(west_df,sample = min(rowSums(west_df)))
west_sncm_stats <- sncm.fit(west_rare,taxon=taxon)
west_sncm_fit <- sncm.fit(west_rare,taxon=taxon,stats = FALSE)

east_df <- east %>% otu_table() %>% as('matrix') %>% as.data.frame()
east_rare <- vegan::rrarefy(east_df,sample = min(rowSums(east_df)))
east_sncm_stats <- sncm.fit(east_rare,taxon=taxon)
east_sncm_fit <- sncm.fit(east_rare,taxon=taxon,stats = FALSE)

west_mean_ra <- 
west %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% as("matrix") %>% as.data.frame()
west_mean_ra[west_mean_ra == 0] <- NA
west_mean_ra <- apply(west_mean_ra,2,function(x){mean(x,na.rm=TRUE)})


east_mean_ra <- 
  east %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  otu_table() %>% as("matrix") %>% as.data.frame()
east_mean_ra[east_mean_ra == 0] <- NA
east_mean_ra <- apply(east_mean_ra,2,function(x){mean(x,na.rm=TRUE)})

west_sncm_fit <- 
west_sncm_fit %>% 
  mutate(mean_relabund = (west_mean_ra),
         neutral_group = case_when(freq < pred.lwr ~ "Below",
                                   freq > pred.upr ~ "Above",
                                   freq >= pred.lwr & freq <= pred.upr ~ "Within"))

east_sncm_fit <- 
  east_sncm_fit %>% 
  mutate(mean_relabund = (east_mean_ra),
         neutral_group = case_when(mean_relabund < pred.lwr ~ "Below",
                                   mean_relabund > pred.upr ~ "Above",
                                   mean_relabund >= pred.lwr & mean_relabund <= pred.upr ~ "Within"))


west_sncm_fit %>% 
  ggplot(aes(x=log10(mean_relabund),y=freq,color=neutral_group)) +
  geom_point()

east_sncm_fit %>% 
  ggplot(aes(x=log10(mean_relabund),y=freq,color=neutral_group)) +
  geom_point(alpha=.5)

# get relative abundance values for merged samples (east vs west)
merged <- east %>% transform_sample_counts(function(x){x/sum(x)})

east <- merged %>% subset_samples(sample_names(merged) == "East") %>% otu_table() %>% as('matrix') %>% t() %>% as.data.frame()
west <- merged %>% subset_samples(sample_names(merged) == "West") %>% otu_table() %>% as('matrix') %>% t() %>% as.data.frame()

which(east > neutral_boundaries$upper)
which(west > neutral_boundaries$upper)

plot(east)

mean(migration_rate)

# build relative abundance table for east and west
relabund_table <- data.frame(west,east,asv=row.names(west)) %>% 
  pivot_longer(-asv,names_to = "location",values_to = "rel_abund")

# plot
relabund_table %>% 
  filter(rel_abund > 0) %>% # remove empty taxa
  ggplot(aes(x=rel_abund,color=location)) +
  geom_density() +
  lims(x=c(0,.0002))
  

neutral_boundaries

# rarefy to minimum samplesum
west_spp <- otu_table(west) %>% as("matrix") %>% as.data.frame()
west_spp_rare <- vegan::rrarefy(west_spp,min(rowSums(west_spp)))

east_spp <- otu_table(east) %>% as("matrix") %>% as.data.frame()
east_spp_rare <- vegan::rrarefy(east_spp,min(rowSums(east_spp)))

# get relative abundance values
west_spp_ra <- vegan::decostand(west_spp_rare,'total',MARGIN = 1) %>% as.data.frame()
east_spp_ra <- vegan::decostand(east_spp_rare,'total',MARGIN = 1) %>% as.data.frame()



