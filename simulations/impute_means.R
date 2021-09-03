rm(list = ls())
library(data.table)
require(doSNOW)
library(mice)
source(file.path("simulations", "utils.R"))

data_wide <- read_wide_data()

data_wide_split <- split(data_wide, data_wide$time)

# average the scores for different items across raters and within
colnms <- paste0("IFBE_", 0:22)

data_wide_mean  <- data_wide[, lapply(.SD, mean), by = list(patient, time), .SDcols=colnms]
other_nms       <- setdiff(colnames(data_wide), colnames(data_wide_mean))
other_nms       <- setdiff(other_nms, c("rater", "rater_group", "patient", "time"))

data_wide_other <- data_wide[, lapply(.SD, function(x) x[!is.na(x)][1]), by = list(patient, time), .SDcols=other_nms]

setorderv(data_wide_other, c("patient", "time"))
setorderv(data_wide_mean,  c("patient", "time"))

data_wide$Aantal_Patienten
data_wide_other$Aantal_Patienten

# returns TRUE, which indicates data frames can be merged by binding the columns
all(
  data_wide_other$patient == data_wide_mean$patient,
  data_wide_other$time    == data_wide_mean$time
)

data_wider <- cbind(
  data_wide_other[seq(1, 210, 2), 1:10],
  data_wide_mean [seq(1, 210, 2), 3:25],
  data_wide_mean [seq(2, 210, 2), 3:25]
)
colnames(data_wider)[11:56] <- paste0(colnms, "_t_", rep(0:1, each = length(colnms)))

sum(is.na(data_wider)) # 21 missing values

hasNA   <- sapply(data_wider, anyNA)
classes <- sapply(data_wider, class)
classes[hasNA]

# impute missing values with mice -- there are no missings in the target variable violent_after
predictorMatrix <- matrix(1L, ncol(data_wider), ncol(data_wider), dimnames = list(colnames(data_wider), colnames(data_wider)))
diag(predictorMatrix)                 <- 0L
predictorMatrix[, "time"]             <- 0L
predictorMatrix[, "Aantal_Patienten"] <- 0L


cores <- parallel::detectCores()
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

pb <- txtProgressBar(min=1, max=8000, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

result <- foreach(i=1:100, .options.snow=opts, .combine='rbind') %dopar% {
  Sys.sleep(1)
}

close(pb)
stopCluster(cl)
doParallel::registerDoParallel(cores = parallel::detectCores())
set.seed(500)

## Parallelized execution
M <- 5
miceout <- foreach(i = seq_len(M), .combine = ibind) %dorng% {
    cat("### Started iteration", i, "\n")
    miceout <- mice(data = df_before, m = 1, print = TRUE,
                    predictorMatrix = predictorMatrix, method = dryMice$method)
    cat("### Completed iteration", i, "\n")
    ## Make sure to return the output
    miceout
}
data_wider_imputed <- mice::mice(data_wider, m = 1, maxit = 1, method = 'pmm', seed = 500, predictorMatrix = predictorMatrix)
data_wider_imputed$loggedEvents
data_wider_imputed

saveRDS(data_wider,         file.path("data", "data_wider.rds"))
saveRDS(data_wider_imputed, file.path("data", "data_wider_imputed.rds"))
