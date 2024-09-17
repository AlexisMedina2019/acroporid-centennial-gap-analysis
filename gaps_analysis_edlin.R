# Load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(janitor)
library(readxl)
library(readr)
library(purrr)
library(fitdistrplus)
library(MASS)

# Exploring theoretical distributions for Monte Carlo simulation ----------
# Raw Data
raw_ages <- read_excel("raw_ages.xlsx")

# Acroporid data of Punta Maroma
age <- raw_ages |>
  filter(source == "This paper")

# All years
ages <- age$t_ka

# Fitting multiple distributions
fit_exponential <- fitdist(ages, "exp")
fit_lognormal <- fitdist(ages, "lnorm")
fit_gamma <- fitdist(ages, "gamma")
fit_uniforme <- fitdist(ages, "unif")

# Q-Q plot
par(mfrow = c(1, 1))
qqcomp(
  list(fit_exponential, fit_lognormal, fit_gamma, fit_uniforme),
  legendtext = c("Exponential", "Log Normal", "Gamma", "uniform")
)

# Histogram with fitted density curves
denscomp(
  list(fit_exponential, fit_lognormal, fit_gamma, fit_uniforme),
  legendtext = c("Exponential", "Log Normal", "Gamma", "uniform")
)

# Model selection with AIC and BIC
aic_values <- c(fit_exponential$aic,
                fit_lognormal$aic,
                fit_gamma$aic,
                fit_uniforme$aic)

bic_values <- c(fit_exponential$bic,
                fit_lognormal$bic,
                fit_gamma$bic,
                fit_uniforme$bic)

# Dataframe
Model.full <- data.frame(
  Data = "All data",
  Distribucion = c("Exponential", "Log Normal", "Gamma", "Uniform"),
  AIC = aic_values,
  BIC = bic_values
)

# Exclusion of the most recent 50 years
ages50 <- age$t_ka[age$t_ka > 0.5]

# Fitting multiple distributions
fit50_exponential <- fitdist(ages50, "exp")
fit50_lognormal <- fitdist(ages50, "lnorm")
fit50_gamma <- fitdist(ages50, "gamma")
fit50_uniforme <- fitdist(ages50, "unif")

# Q-Q plot
qqcomp(
  list(fit50_exponential, fit50_lognormal, fit50_gamma, fit50_uniforme),
       legendtext = c("Exponential", "Log Normal", "Gamma", "uniform")
  )
  
# Histogram with fitted density curves
denscomp(
  list(fit50_exponential, fit50_lognormal, fit50_gamma, fit50_uniforme),
    legendtext = c("Exponential", "Log Normal", "Gamma", "uniform")
  )
  
# Comparar AIC y BIC
  aic_values50 <- c(fit50_exponential$aic,
                    fit50_lognormal$aic,
                    fit50_gamma$aic,
                    fit50_uniforme$aic)
  
  bic_values50 <- c(fit50_exponential$bic,
                    fit50_lognormal$bic,
                    fit50_gamma$bic,
                    fit50_uniforme$bic)
  
  # Model selection with AIC and BIC
Model.50 <- data.frame(
    Data = "exclusion of 50 most recent years",
    Distribucion = c("Exponential", "Log Normal", "Gamma", "Uniform"),
    AIC = aic_values50,
    BIC = bic_values50
  )
  
models <- bind_rows(Model.full, Model.50)
#Uniform distribution allows for a better fit in both scenarios  

# Simulation using Uniform distribution ----------------------
## Yucatan
raw_ages <- read_excel("raw_ages.xlsx")
acrop <-raw_ages|>
  filter(source == "This paper")|>
  filter(species_2 == "BRAN")

age_acrop <- acrop$t_ka

age_range <- 1000 * (max(age_acrop) - min(age_acrop))

intervals <- 50 #length of class intervals in years

CI <- round(age_range/intervals, 0) # numbre of classes

# Histogram to visualize the data
hist(age_acrop,
     breaks = CI, 
     main = paste("Yucatan,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Defining class intervals
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Yucatan,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calculate the observed frequencies
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Function to identify sequences of zeros and their length
gaps <- function(vector) {
  secuencias <- data.frame(inicio = integer(), longitud = integer())
  longitud_actual <- 0
  inicio_secuencia <- NULL
  
  for (i in 1:length(vector)) {
    if (vector[i] == 0) {
      if (is.null(inicio_secuencia)) {
        inicio_secuencia <- i
      }
      longitud_actual <- longitud_actual + 1
    } else {
      if (!is.null(inicio_secuencia)) {
        secuencias <- rbind(secuencias, data.frame(inicio = inicio_secuencia, longitud = longitud_actual))
      longitud_actual <- 0
        inicio_secuencia <- NULL
      }
    }
  }
  
  if (!is.null(inicio_secuencia)) {
    secuencias <- rbind(secuencias, data.frame(inicio = inicio_secuencia, longitud = longitud_actual))
  }
  
  secuencias$size <- secuencias$longitud * (age_range/length(observed))
  
  
  return(secuencias)
}

# Monte Carlo Simulation
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  # Parameterized Uniform Distribution
  simulated_data <- runif(length(age_acrop), min = min(age_acrop), max = max(age_acrop))
  
  # Calculate the simulated frequencies
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  # Gap tests
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis
gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

yucatan <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))|>
  mutate(site = "Yucatan")

## Barbados
acrop <-raw_ages|>
  filter(source == "Abdul et al.,16")|>
  filter(species_2 == "BRAN")

age_acrop <- acrop$t_ka

age_range <- 1000 * (max(age_acrop) - min(age_acrop))

intervals <- 50

CI <- round(age_range/intervals, 0)

# Histogram to visualize the data
hist(age_acrop,
     breaks = CI, 
     main = paste("Barbados,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Defining class intervals
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Barbados,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calculate the observed frequencies
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Monte Carlo Simulation
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  simulated_data <- runif(81, min = min(age_acrop), max = max(age_acrop))
  
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis

gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

barbados <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))|>
  mutate(site = "Barbados")

## Belize
acrop <-raw_ages|>
  filter(source == "Gischler et al.,23")|>
  filter(species_2 == "BRAN")

age_acrop <- acrop$t_ka

age_range <- 1000 * (max(age_acrop) - min(age_acrop))

intervals <- 50

CI <- round(age_range/intervals, 0)

# Histogram to visualize the data
hist(age_acrop,
     breaks = CI, 
     main = paste("Barbados,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Defining class intervals
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Barbados,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calculate the observed frequencies
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Monte Carlo simulation
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  simulated_data <- runif(69, min = min(age_acrop), max = max(age_acrop))
  
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts

  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis

gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

belize <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))|>
  mutate(site = "Belize")

# Saint Croix
acrop <-raw_ages|>
  filter(source == "Hubbard et al.,05;13")|>
  filter(species_2 == "BRAN")

age_acrop <- acrop$t_ka

age_range <- 1000 * (max(age_acrop) - min(age_acrop))

intervals <- 50

CI <- round(age_range/intervals, 0)

# Histogram to visualize the data
hist(age_acrop,
     breaks = CI, 
     main = paste("Saint Croix,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Defining class intervals
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Saint Croix,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calculate the observed frequencies
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Monte Carlo simulations
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  simulated_data <- runif(52, min = min(age_acrop), max = max(age_acrop))
  
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis
gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

stc <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))|>
  mutate(site = "Saint Croix")

## All probabilities in a single plot 

bind_rows(yucatan, barbados, belize, stc) |>
  filter(size <= 1250) |>
  mutate(Region = factor(
    site,
    levels = c("Yucatan", "Barbados", "Belize", "Saint Croix"),
    labels = c(
      "Yucatan (n = 64)",
      "Barbados (n = 81)",
      "Belize (n = 69)",
      "Saint Croix (n = 52)"
    )
  )) |>
  ggplot(aes(x = size, y = freq)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_x_continuous(breaks = seq(50, 1250, 50)) +
  geom_hline(yintercept = 0.05,
             colour = "red",
             linewidth = 1) +
  ylab(expression("Probability under " * H[0] * " of uniform distribution") ) +
  xlab("Gap size (years)") +
  geom_line(aes(colour = Region), linewidth = 1) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black")
  )


bind_rows(yucatan, barbados, belize, stc) |>
  write_csv(file = "Table_gap_size_pvalues.csv")


# Large gaps in a composite analyses of all sites -------------------------
raw_ages <- read_excel("raw_ages.xlsx")
acrop <-raw_ages|>
  filter(species_2 == "BRAN")

age_acrop <- acrop$t_ka

age_range <- 1000 * (max(age_acrop) - min(age_acrop))

intervals <- 50

CI <- round(age_range/intervals, 0)

# Histogram to visualize the data
hist(age_acrop,
     breaks = CI, 
     main = paste("All sites,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Defining class intervals
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("All sites,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calculate the observed frequencies
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Monte Carlo simulation
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  simulated_data <- runif(266, min = min(age_acrop), max = max(age_acrop))
  
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis

gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

all_sites <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))|>
  mutate(site = "All sites")


bind_rows(yucatan, barbados, belize, stc, all_sites) |>
  filter(size <= 1250) |>
  mutate(Region = factor(
    site,
    levels = c("Yucatan", "Barbados", "Belize", "Saint Croix", "All sites"),
    labels = c(
      "Yucatan (n = 64)",
      "Barbados (n = 81)",
      "Belize (n = 69)",
      "Saint Croix (n = 52)",
      "All sites (n = 266)"
    )
  )) |>
  ggplot(aes(x = size, y = freq)) +
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_x_continuous(breaks = seq(50, 1250, 50)) +
  geom_hline(yintercept = 0.05,
             colour = "red",
             linewidth = 1) +
  ylab(expression("Probability under " * H[0] * " of uniform distribution") ) +
  xlab("Gap size (years)") +
  geom_line(aes(colour = Region), linewidth = 1) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black")
  )
