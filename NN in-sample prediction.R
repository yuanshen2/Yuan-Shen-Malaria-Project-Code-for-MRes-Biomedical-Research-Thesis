# Load required libraries
library(tensorflow)
library(caret)
library(ggplot2)

# Set the Python version for reticulate
use_python("C:/Users/Administrator/AppData/Local/Programs/Python/Python38/python.exe")

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
tf$random$set_seed(initial_seed)

covariate.names <- inputs$covariate.names
original_pcent_mortality <- data_all[,"pcent_mortality"]

# Empirical logit and IHS transform on labels
data_all[,"pcent_mortality"] <- emplogit2(data_all$no_dead, data_all$no_tested)
theta2 <- optimise(IHS.loglik, lower = 0.001, upper = 50, x = data_all[,"pcent_mortality"], maximum = TRUE)$maximum
print(paste("theta2 maximum", theta2, sep=" "))
data_all[,"pcent_mortality"] <- IHS(data_all[,"pcent_mortality"], theta2)

# Splitting data into input and output
cols_to_remove <- c('start_month', 'end_year', 'end_month', 'latitude', 'longitude',
                    'pcent_mortality', 'no_dead', 'no_tested', 'cellnumber')

data_x <- as.matrix(data_all[, !(colnames(data_all) %in% cols_to_remove)])
data_x <- scale(data_x)  # Normalize the data
data_y <- as.numeric(data_all$pcent_mortality)

# Define the Keras model using tf$keras
create_model <- function(neurons=5) {
  model <- tf$keras$models$Sequential()
  model$add(tf$keras$layers$Dense(units = neurons, activation = 'elu', input_shape = list(ncol(data_x))))
  model$add(tf$keras$layers$Dense(units = 1))
  model$compile(
    loss = 'mean_squared_error',
    optimizer = 'adam'
  )
  return(model)
}


# Define custom train function to integrate Keras with caret
early_stopping <- tensorflow::tf$keras$callbacks$EarlyStopping(monitor = "loss", patience = 20)
train_model <- function(x, y, neurons, batch_size) {
  model <- create_model(neurons)
  history <- model$fit(
    x = x, y = y,
    epochs = as.integer(1000),
    batch_size = batch_size,
    validation_split = 0, 
    verbose = 0,
    callbacks = list(early_stopping)
  )
  return(model)
}

# Function to apply inverse transformations to the predictions
inverse_transformations <- function(predictions, theta2) {
  inv_ihs <- function(x, theta) {
    return((1/theta)*sinh(theta * x))
  }
  inv_emplogit <- function(x) {
    return(exp(x) / (1 + exp(x)))
  }
  return(inv_emplogit(inv_ihs(predictions, theta2)))
}

# Convert entire data to numpy arrays
np <- import("numpy")
data_x_np <- np$array(data_x, dtype="float32")
data_y_np <- np$array(data_y, dtype="float32")

# Train model on entire dataset
model <- train_model(data_x_np, data_y_np, as.integer(160), as.integer(254))

# Predict on the entire dataset
predictions <- model$predict(data_x)

# Apply the inverse transformations to predictions
predictions_transformed <- inverse_transformations(predictions, theta2)

# Scale the original observed values
scaled_original_pcent_mortality <- original_pcent_mortality / 100

# Extract transformed predictions for DRC_data using DRC_indices
DRC_predictions_transformed <- predictions_transformed[DRC_indices]

# View the transformed predictions for DRC_data
print(DRC_predictions_transformed)

# Extract scaled_original_pcent_mortality values for DRC_data using DRC_indices
DRC_scaled_original_pcent_mortality <- scaled_original_pcent_mortality[DRC_indices]

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
  geom_point(color = "black") +  # Removed alpha argument
  geom_smooth(method = "lm", col = "red") +
  xlab("Observed Mortality") + 
  ylab("Predicted Mortality") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 30),  # Bold and increase size for x axis ticks
    axis.text.y = element_text(face = "bold", size = 30),  # Bold and increase size for y axis ticks
    axis.title.x = element_text(size = 32),                # Increase size for x axis title
    axis.title.y = element_text(size = 32),                # Increase size for y axis title
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),   # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Adds a border to the plotting area
    axis.line = element_line(color = "black"),  # Adds axis lines
    axis.ticks = element_line(color = "black")  # Adds tick marks to the axes
  )

print(DRC_plot)


#Get the R2 value for NN in-sample prediction
# Fit a linear regression model
lm_DRC <- lm(Predicted ~ Observed, data = DRC_results_df)

# Extract the R-squared value from the model summary
R2_DRC <- summary(lm_DRC)$r.squared

print(R2_DRC)


#Get the NN in-sample prediction linerar regression
# Fit a linear model to the data
lm_fit <- lm(Predicted ~ Observed, data = DRC_results_df)

# Retrieve the coefficients (intercept and slope)
intercept <- coef(lm_fit)[1]
slope <- coef(lm_fit)[2]

# Print the equation
cat("Linear regression equation: y =", intercept, "+", slope, "* x")

#Get the NN in-sample prediction RMSE
# Fit a linear model to the data
lm_fit <- lm(Predicted ~ Observed, data = DRC_results_df)

# Compute the residuals
residuals <- lm_fit$residuals

# Calculate the RMSE
RMSE <- sqrt(mean(residuals^2))

# Print the RMSE
print(RMSE)


