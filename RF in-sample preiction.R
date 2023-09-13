library(caret)
library(data.table)
library(ggplot2)
# Load datasets and source functions
source("C:/Users/Administrator/Documents/malaria project/load_gen.r")
load("inputs_data_all_wa_tree5d.r")
load("stk_val_inds5.r")
load("non_spatial_val_sets6.r")
load("~/malaria project/malariadata.RData")

# Extract the data
Malariadata <- inputs$data_all_wa_ea_drc
data_all <- na.aggregate(Malariadata)
DRC_data <- Malariadata[DRC_indices, ]

# Initial seed
initial_seed <- 90

# Set seeds for reproducibility
set.seed(initial_seed)

covariate.names <- inputs$covariate.names
original_pcent_mortality <- data_all[,"pcent_mortality"]

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
colnames(rfTuneGrid) <- "mtry"
paramsTuneGrid <- data.frame(col_sample_rate_per_tree = 0.7, nodesize = 30)

# Train the model on the entire dataset
rfFitJ_all <- train(form, data = data_all, method = "rf", verbose = TRUE, tuneGrid = rfTuneGrid, ntree = 100, replace = TRUE, metric = "RMSE", col_sample_rate_per_tree = paramsTuneGrid$col_sample_rate_per_tree, nodesize = paramsTuneGrid$nodesize)

# Predict on the same dataset
rfPred_all <- predict(rfFitJ_all, newdata = data_all)

# Calculate RMSE for in-sample prediction
rmse_all <- sqrt(mean((rfPred_all - data_all$pcent_mortality)^2))
print(paste("In-sample RMSE:", rmse_all))

# New transformation function
transform_predictions <- function(predictions, theta2) {
  inv_ihs <- function(x, theta) {
    return((1/theta)*sinh(theta * x))
  }
  inv_emplogit <- function(x) {
    return(exp(x) / (1 + exp(x)))
  }
  return(inv_emplogit(inv_ihs(predictions, theta2)))
}

# Apply the new transformation function
rfPred_all_inv <- transform_predictions(rfPred_all, theta2$maximum)
actual_all_inv <- transform_predictions(data_all$pcent_mortality, theta2$maximum)
# Extract transformed predictions for DRC_data using DRC_indices
DRC_predictions_transformed <- rfPred_all_inv[DRC_indices]

# View the transformed predictions for DRC_data
print(DRC_predictions_transformed)
# Extract scaled_original_pcent_mortality values for DRC_data using DRC_indices
DRC_scaled_original_pcent_mortality <- actual_all_inv[DRC_indices]

# View the scaled_original_pcent_mortality values for DRC_data
print(DRC_scaled_original_pcent_mortality)

# Calculate RMSE for DRC_data
rmse_DRC <- sqrt(mean((DRC_predictions_transformed - DRC_scaled_original_pcent_mortality)^2))
print(paste("RMSE for DRC_data:", rmse_DRC))


# Create a dataframe for DRC_data's observed and predicted values
DRC_results_df <- data.frame(
  Observed = DRC_scaled_original_pcent_mortality,
  Predicted = as.vector(DRC_predictions_transformed)
)

# Plot observed vs. predicted values for DRC_data
DRC_plot <- ggplot(DRC_results_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black") + 
  geom_smooth(method = "lm", col = "red") +
  xlab("Observed Mortality") + 
  ylab("Predicted Mortality") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "plain", size = 30, color = "black"),  # Adjusts the x-axis text
    axis.text.y = element_text(face = "plain", size = 30, color = "black"),  # Adjusts the y-axis text
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Adds a border to the plot
    panel.background = element_rect(fill = "white"),  # Sets the plot background to white
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),   # Removes minor grid lines
    axis.line = element_line(colour = "black"),  # Adds the axis lines
    axis.ticks = element_line(colour = "black")   # Adds the axis ticks
  )

print(DRC_plot)

#Get the R2 value for RF DRC in sample prediction
# Fit a linear regression model
lm_DRC <- lm(Predicted ~ Observed, data = DRC_results_df)

# Extract the R-squared value from the model summary
R2_DRC <- summary(lm_DRC)$r.squared

print(R2_DRC)


#Get the RF in sample prediction linear regression line
# Fit a linear model
lm_fit <- lm(Predicted ~ Observed, data = DRC_results_df)

# Extract coefficients
intercept <- coef(lm_fit)[1]
slope <- coef(lm_fit)[2]

# Display the equation
equation <- sprintf("y = %.2fx + %.2f", slope, intercept)
print(equation)

#Get the RMSE value for the RF in sample prediction
# Fit a linear model
lm_fit <- lm(Predicted ~ Observed, data = DRC_results_df)

# Compute residuals
residuals <- resid(lm_fit)

# Calculate RMSE
RMSE <- sqrt(mean(residuals^2))
print(RMSE)


