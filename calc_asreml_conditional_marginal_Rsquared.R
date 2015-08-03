# Get design matrix of fixed effects from model
Fmat <- model.matrix(eval(asreml.fit$fixed.formula)[-2], glm.spl)

# Get variance of fixed effects by multiplying coefficients by design matrix
VarF <- sum(var(as.vector(rev(asreml.fit$coefficients$fixed) %*% t(Fmat))), summary(asreml.fit)$varcomp[1:5, 2])

# Get variance of random effects by extracting variance components
VarRand <- sum(summary(asreml.fit)$varcomp[6, 2])

# Get residual variance
VarResid <- summary(asreml.fit)$varcomp[7, 2]

varTotal <- VarF + VarRand + VarResid

#calculate marginal R-squared
marR2 <- VarF/varTotal

#calculate conditional R-squared
condR2 <- (VarF + VarRand)/varTotal

cbind(VarF, VarRand, VarResid)



