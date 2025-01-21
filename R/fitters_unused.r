# # https://stats.stackexchange.com/questions/645614/reimplementation-of-massglm-nb-with-fastglm
# speedglm.fit.nb <- function(y, X, .formula, .offset, sparse, maxit, tol = 1e-07) {

#     cat("\tInitializing alternation...\n")
#     fit0 <- speedglm::speedglm.wfit(
#         y,
#         X,
#         family = poisson(),
#         sparse = sparse,
#         offset = .offset,
#         maxit = maxit
#     )
#     mu0 <- predictSparse(fit0, X, .offset)
#     theta0 <- MASS::theta.ml(y, mu0, limit = maxit)
#     fam0 <- MASS::negative.binomial(theta = as.vector(theta0))
    
#     iter <- 0
#     dev <- 2
#     while(dev > tol | iter <= maxit){

#         iter <- iter + 1
#         cat("\tAlternation ", iter, ", current theta ", signif(theta0, 5), "...\n", sep = "")
#         rm(fit0)

#         fit0 <- speedglm::speedglm.wfit(
#             y,
#             X,
#             family = fam0,
#             sparse = sparse,
#             offset = .offset,
#             maxit = maxit,
#             etastart = log(mu0),
#             trace = TRUE
#         )
#         mu0 <- predictSparse(fit0, X, .offset)
#         theta1 <- MASS::theta.ml(y, mu0, limit = theta.maxit)
#         fam0 <- MASS::negative.binomial(theta = theta1)
        
#         dev <- abs(theta1 - theta0)
#         theta0 <- theta1

#     }

#     l <- list(model = fit0, theta = theta0)
#     return(l)

# }

# #' @export
# predictSparse <- function(x, newdata, ...) {

#    UseMethod("predictSparse", x)

# }

# #' @export
# predictSparse.speedglm <- function(model, newdata, .offset) {

#     b <- coef(model)
#     b[is.na(b)] <- 0.0
#     mu <- model$family$linkinv((newdata %*% b)[, 1] + .offset)

#     return(mu)

# }
