library(caret)
library(data.table)
library(ggplot2)
source("C:/Users/Administrator/Documents/malaria project/load_gen.r")
load("C:/Users/Administrator/Documents/malaria project/malariadata.RData")

Malariadata <- inputs$data_all_wa_ea_drc
data_all <- na.aggregate(Malariadata)
DRC_data <- Malariadata[DRC_indices, ]

# Set a seed for reproducibility
set.seed(90)

covariate.names <- inputs$covariate.names
length(covariate.names)

# Empirical logit and IHS transform on labels
data_all[,"pcent_mortality"] <- emplogit2(data_all["no_dead"],data_all["no_tested"])
theta2 <- optimise(IHS.loglik, lower = 0.001, upper = 50, x = data_all[,"pcent_mortality"],maximum = TRUE)

print(paste("theta2 maximum", theta2$maximum,sep=" "))
data_all[,"pcent_mortality"] <- IHS(data_all[,"pcent_mortality"],theta2$maximum)

model_lhs <- paste(c(inputs$covariate.names),collapse='+')

#Format for RF
form <- as.formula(paste("pcent_mortality ~", model_lhs))

# Specify hyperparameters
rfTuneGrid <- data.frame(mtry = 150)
rfPredJ <- list()
actual_values <- list()

colnames(rfTuneGrid) <- "mtry"
paramsTuneGrid <- data.frame(col_sample_rate_per_tree = 0.7, nodesize = 20)

# Create 10 folds of the data
folds <- createFolds(data_all$no_tested, k = 10)

# Run the entire process 10 times
for(i in 1:10) {
  set.seed(i)  # Set a different seed each time
  
  # Create indices for the test set
  test_inds <- folds[[i]]
  
  # Subset your data to create the test set
  test_set <- data_all[test_inds, ]
  
  # Store actual values
  actual_values[[i]] <- test_set$pcent_mortality
  
  # Create indices for the training set
  train_inds <- setdiff(seq_len(nrow(data_all)), test_inds)
  
  # Subset your data to create the training set
  train_set <- data_all[train_inds, ]
  
  # Replace missing values with the mean value of the respective feature
  preProcValues <- preProcess(train_set, method = 'medianImpute')
  train_set <- predict(preProcValues, train_set)
  
  #Train the model
  rfFitJ <- train(form, data = train_set, method = "rf", verbose = TRUE, tuneGrid = rfTuneGrid, ntree = 100, replace = TRUE, metric = "RMSE", col_sample_rate_per_tree = paramsTuneGrid$col_sample_rate_per_tree, nodesize = paramsTuneGrid$nodesize)
  
  #Get out of sample predictions
  rfPredJ[[i]] <- predict(rfFitJ,newdata=test_set)
  
  # Calculate RMSE for each fold
  rmse <- sqrt(mean((rfPredJ[[i]] - actual_values[[i]])^2))
  print(paste("Fold", i, "RMSE:", rmse))
}

# Calculate overall RMSE
all_preds <- do.call(c, rfPredJ)
all_actuals <- do.call(c, actual_values)
overall_rmse <- sqrt(mean((all_preds - all_actuals)^2))
print(paste("Overall RMSE:", overall_rmse))

# Inverse of the empirical logit transformation
inv_emplogit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

# Inverse of the IHS transformation
inv_ihs <- function(x, theta) {
  return((1/theta)*sinh(theta * x))
}

# Apply the inverse transformations to the predictions and actual values
all_preds_inv <- inv_emplogit(inv_ihs(all_preds, theta2$maximum))
all_actuals_inv <- inv_emplogit(inv_ihs(all_actuals, theta2$maximum))
DRC_predictions_transformed <- all_preds_inv[DRC_indices]
DRC_actuals_transformed<-all_actuals_inv[DRC_indices]


# Plot the actual values vs predicted values 
library(gridExtra)  # for arranging plots
plots <- list()  # to store the plots

#Plot the RF out-of-sample prediction for DRC
df_DRC <- data.frame(Predicted = DRC_predictions_transformed, Actual = DRC_actuals_transformed)

p_DRC <- ggplot(df_DRC, aes(x = Actual, y = Predicted)) + 
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "solid") + # Set linetype to solid
  xlab("Observed Mortality") +
  ylab("Predicted Mortality") +
  theme(
    axis.text.x = element_text(face = "bold", size = 30), 
    axis.text.y = element_text(face = "bold", size = 30),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    panel.background = element_rect(fill = "white"),  # This line sets the plot background to white
    panel.border = element_rect(colour = "black", fill = NA, size = 1)  # Adds a border to the plot
  )

print(p_DRC)

#Get the R2 Value for RF out of sample prediction
# Fit a linear regression model
lm_DRC <- lm(Predicted ~ Actual, data = df_DRC)

# Extract the R-squared value from the model summary
R2_DRC <- summary(lm_DRC)$r.squared

print(R2_DRC)


#Get the linear regression line for RF out of sample prediction
# Fit a linear model
lm_fit <- lm(Predicted ~ Actual, data = df_DRC)

# Extract coefficients
coefficients <- coef(lm_fit)
intercept <- coefficients[1]
slope <- coefficients[2]

# Print the equation of the line
equation <- sprintf("y = %.2fx + %.2f", slope, intercept)
print(equation)

#Get the RMSE value for RF out of sample prediction
# Fit a linear model
lm_fit <- lm(Predicted ~ Actual, data = df_DRC)

# Get predicted values
predictions <- predict(lm_fit, df_DRC)

# Calculate residuals
residuals <- df_DRC$Predicted - predictions

# Calculate RMSE
rmse <- sqrt(mean(residuals^2))
print(rmse)


