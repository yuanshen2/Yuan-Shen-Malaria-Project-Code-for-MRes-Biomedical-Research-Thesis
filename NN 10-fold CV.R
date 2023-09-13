# Load required libraries
library(keras)
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

# Define the Keras model
create_model <- function(neurons=5) {
  model <- keras_model_sequential()
  model$add(layer_dense(units = neurons, activation = 'elu', input_shape = ncol(data_x)))
  model$add(layer_dense(units = 1))
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

# Manual 10-fold cross-validation
num_folds <- 10
fold_size <- floor(nrow(data_x) / num_folds)

# Create a data frame to accumulate results across folds
results_df <- data.frame(Observed = numeric(), Predicted = numeric(), Fold = factor())

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

# Loop through each fold and store the results in results_df
np <- import("numpy")
for(i in 1:num_folds) {
  # Set a unique seed for each fold
  set.seed(initial_seed + i)
  
  val_indices <- ((i-1) * fold_size + 1):min(i * fold_size, nrow(data_x))
  train_x_fold <- data_x[-val_indices,]
  train_y_fold <- data_y[-val_indices]
  val_x_fold <- data_x[val_indices,]
  val_y_fold <- data_y[val_indices]
  
  # Convert to numpy arrays
  train_x_np <- np$array(train_x_fold, dtype="float32")
  train_y_np <- np$array(train_y_fold, dtype="float32")
  
  model <- train_model(train_x_np, train_y_np, as.integer(160), as.integer(254))
  predictions <- model$predict(val_x_fold)
  
  # Apply the inverse transformations to predictions
  predictions_transformed <- inverse_transformations(predictions, theta2)
  
  # Scale the original observed values
  scaled_original_pcent_mortality <- original_pcent_mortality[val_indices] / 100
  
  # Store the scaled observed values and transformed predictions
  fold_df <- data.frame(Observed = scaled_original_pcent_mortality, 
                        Predicted = as.vector(predictions_transformed), 
                        Fold = factor(i))
  results_df <- rbind(results_df, fold_df)
}

# Now plot using ggplot2 with facet_wrap for the 4x3 layout
peachfold <- ggplot(results_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black", alpha = 0.6) + # Set color to black
  geom_smooth(method = "lm", col = "red") +
  facet_wrap(~Fold, ncol = 4) +
  # ggtitle("NN Observed vs Predicted Mortality Values by Fold") + # This line is commented out
  xlab("Observed Mortality") + ylab("Predicted Mortality") + 
  theme_minimal() +
  guides(color = FALSE) # Remove the color legend

print(peachfold)

# Plotting observed vs predicted values across all folds without title
p_all_folds <- ggplot(results_df, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black") + # Set color to black
  geom_smooth(method = "lm", col = "red") +
  xlab("Observed Mortality") + 
  ylab("Predicted Mortality") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 30), 
    axis.text.y = element_text(face = "bold", size = 30),
    axis.title.x = element_text(size = 32),
    axis.title.y = element_text(size = 32)
  ) +
  guides(color = FALSE) # Remove the color legend

print(p_all_folds)

# Compute RMSE
rmse <- sqrt(mean((results_df$Observed - results_df$Predicted)^2))
print(paste("Overall RMSE:", round(rmse, 4)))

#Calculate R2 values for overall
# Fit a linear model
model <- lm(Predicted ~ Observed, data = results_df)

# Extract R^2
r_squared <- summary(model)$r.squared
print(r_squared)

#Calculate RMSE for all folds
# Compute the squared differences
squared_differences <- (results_df$Observed - results_df$Predicted)^2

# Calculate the mean of the squared differences
mean_squared_error <- mean(squared_differences)

# Compute the RMSE
rmse <- sqrt(mean_squared_error)

# Print the RMSE
print(rmse)


# Fit a linear regression model
model <- lm(Predicted ~ Observed, data = results_df)

# Print the coefficients
print(coef(model))

#get the rank of feature importance
# Get baseline RMSE using validation set (for simplicity, using data_x as a placeholder for validation set)
baseline_predictions <- model$predict(data_x)
compute_rmse <- function(true, predicted) {
  sqrt(mean((true - predicted)^2))
}

baseline_rmse <- compute_rmse(data_y, baseline_predictions)


# Initialize vector to hold importance values
feature_importance <- numeric(ncol(data_x))

# Compute importance for each feature
for (i in 1:ncol(data_x)) {
  set.seed(initial_seed + i)  # Setting a unique seed for each feature shuffle based on the initial_seed.
  
  shuffled_data <- data_x
  shuffled_data[,i] <- sample(data_x[,i])
  
  shuffled_predictions <- model$predict(shuffled_data)
  shuffled_rmse <- compute_rmse(data_y, shuffled_predictions)
  
  feature_importance[i] <- shuffled_rmse - baseline_rmse
}

# Create a dataframe of features and their importance
importance_df <- data.frame(Variable = colnames(data_x), Importance = feature_importance)

# Extract top 20 important features using dplyr
library(dplyr)
top_20_importance <- importance_df %>%
  arrange(-Importance) %>%
  head(20)

print(top_20_importance)
