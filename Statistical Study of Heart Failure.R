# Data Preparation
# Load the heart dataset
heart.df <- read.csv('heart3.csv')

# Convert categorical variables to factors
heart.df$Sex <- as.factor(heart.df$Sex)
heart.df$ChestPainType <- as.factor(heart.df$ChestPainType)
heart.df$RestingECG <- as.factor(heart.df$RestingECG)
heart.df$ExerciseAngina <- as.factor(heart.df$ExerciseAngina)
heart.df$ST_Slope <- as.factor(heart.df$ST_Slope)

# Summary of the dataset to check for missing or unusual values
summary(heart.df)

# Replace 0 values with NA for variables where 0 is not meaningful
heart.df$Cholesterol[heart.df$Cholesterol == 0] <- NA
heart.df$RestingBP[heart.df$RestingBP == 0] <- NA

# Imputation of missing values using MICE (Predictive Mean Matching)
# Ensure the mice package is installed before running this step
# install.packages("mice")
library(mice)
heart <- complete(mice(heart.df, method = "pmm"))

# Replace negative values in "Oldpeak" with 0
heart$Oldpeak[heart$Oldpeak < 0] <- 0

# Summary of the cleaned dataset
summary(heart)

# Backward elimination to optimize the logistic regression model
# Fit the full model with all predictors
fullbinmodel <- glm(HeartDisease ~ ., data = heart, family = binomial)
summary(fullbinmodel)

# Iteratively remove predictors with the least significance
mod2 <- update(fullbinmodel, . ~ . - RestingECG)
summary(mod2)

mod3 <- update(mod2, . ~ . - RestingBP)
summary(mod3)

mod4 <- update(mod3, . ~ . - Cholesterol)
summary(mod4)

mod5 <- update(mod4, . ~ . - MaxHR)
summary(mod5)

# Final logistic regression model with interaction term
heart.binmodel <- glm(HeartDisease ~ Age + Sex + ChestPainType + FastingBS +
                        ExerciseAngina + RestingECG + ChestPainType:RestingECG + Oldpeak + ST_Slope,
                      data = heart, family = binomial)
summary(heart.binmodel)

# Plotting diagnostics for the final model
plot(heart.binmodel, which = 4)

# Reduced model without interaction term
redmodel <- glm(HeartDisease ~ Age + Sex + ChestPainType + FastingBS +
                  ExerciseAngina + RestingECG + Oldpeak + ST_Slope,
                data = heart, family = binomial)

# Likelihood Ratio Test (LRT) to compare models
LRT <- deviance(redmodel) - deviance(heart.binmodel)
p_value <- 1 - pchisq(LRT, df.residual(redmodel) - df.residual(heart.binmodel))
p_value  # P-value for the LRT

# Train-test split (Assuming train_part and test_part are pre-defined)
# Fit the model on training data
heart_train.binmodel <- glm(HeartDisease ~ Age + Sex + ChestPainType + FastingBS +
                              ExerciseAngina + RestingECG + ChestPainType:RestingECG + Oldpeak + ST_Slope,
                            data = train_part, family = binomial)
summary(heart_train.binmodel)

# Evaluate model performance using ROC curve
library(pROC)
predictedprob <- predict(heart_train.binmodel, test_part, type = "response")
roc_heart <- roc(test_part$HeartDisease, predictedprob, plot = TRUE)
roc_heart  # Displays AUC and plot

# Determine the best cutoff point
best_cutoff_heart <- coords(roc_heart, "best", ret = "threshold")
best_cutoff_heart

# Convert probabilities to binary predictions based on the cutoff
yhat <- rep(0, nrow(test_part))
yhat[predictedprob >= best_cutoff_heart] <- 1

# Confusion matrix to evaluate predictions
table(test_part$HeartDisease, yhat)

# Hosmer-Lemeshow test to assess goodness-of-fit
hosmerlem <- function(y, yhat, g = 10) {
  cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0, 1, 1/g)), include.lowest = TRUE)
  obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq <- sum((obs - expect)^2 / expect)
  P <- 1 - pchisq(chisq, g - 2)
  c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
}

hosmerlem(train_part$HeartDisease, heart_train.binmodel$fitted.values)

