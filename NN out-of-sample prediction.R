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


results_df$IsDRC <- ifelse(row.names(results_df) %in% row.names(DRC_data), "DRC", "Others")
View(results_df$IsDRC)
# Subset the data for DRC
DRC_subset <- results_df[results_df$IsDRC == "DRC", ]

# Extract the observed and predicted values
DRC_observed <- DRC_subset$Observed
DRC_predicted <- DRC_subset$Predicted
DRC_dataframe <- data.frame(Observed = DRC_observed, Predicted = DRC_predicted)

p_DRC <- ggplot(DRC_dataframe, aes(x = Observed, y = Predicted)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", col = "red") +
  xlab("Observed Mortality") + ylab("Predicted Mortality") + 
  theme_minimal() +
  guides(color = "none") + 
  theme(
    axis.text = element_text(size = 30, face = "plain"),  # Increases size and changes the axis tick labels to plain
    axis.title.x = element_text(size = 32, face = "plain"),  # Increases size and changes the x axis title to plain
    axis.title.y = element_text(size = 32, face = "plain"),   # Increases size and changes the y axis title to plain
    panel.background = element_rect(fill = "white"),  # This line sets the plot background to white
    plot.background = element_rect(fill = "white"),   # This line sets the entire plot area (including margin) background to white
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank(),   # Removes minor grid lines
    panel.border = element_blank(),  # This line removes the outer border
    axis.line = element_line(colour = "black", size = 1),  # Adds axes lines
    axis.ticks = element_line(colour = "black", size = 1),  # Adds the tick marks on both x and y axes
    axis.text.x = element_text(colour = "black", size = 30, face = "plain"),  # Adjusts the x-axis text to plain
    axis.text.y = element_text(colour = "black", size = 30, face = "plain")   # Adjusts the y-axis text to plain
  )

print(p_DRC)


#Get the R2 value for NN Out-of-sample prediction
# Fit a linear regression model
DRC_lm <- lm(Predicted ~ Observed, data = DRC_dataframe)

# Extract the R-squared value from the model summary
R2 <- summary(DRC_lm)$r.squared

print(R2)


#Get the linear regression line of NN out of sample prediction

# Fit a linear model
lm_fit <- lm(Predicted ~ Observed, data = DRC_dataframe)

# Extract coefficients
coefficients <- coef(lm_fit)
intercept <- coefficients[1]
slope <- coefficients[2]

# Print the equation of the line
equation <- sprintf("y = %.2fx + %.2f", slope, intercept)
print(equation)

#Get the RMSE value for NN out-of-sample prediction
# Fit a linear model
lm_fit <- lm(Predicted ~ Observed, data = DRC_dataframe)

# Compute residuals
residuals <- resid(lm_fit)

# Calculate RMSE
RMSE <- sqrt(mean(residuals^2))
print(RMSE)
