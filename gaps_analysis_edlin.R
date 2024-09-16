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

raw_ages <- read_excel("raw_ages.xlsx")
age<-raw_ages|>
  filter(source == "This paper")

edades <- age$t_ka[age$t_ka > 0.5] # Se excluyen los 50 años más recientes

# Ajustar varias distribuciones
ajuste_exponential <- fitdist(edades, "exp")
ajuste_lognormal <- fitdist(edades, "lnorm")
ajuste_gamma <- fitdist(edades, "gamma")
ajuste_uniforme <- fitdist(edades, "unif")

# Gráficos Q-Q para comparar
par(mfrow = c(1, 1))  # Múltiples gráficos en una ventana

qqcomp(list(ajuste_exponential, 
            ajuste_lognormal,
            ajuste_gamma,
            ajuste_uniforme
            #ajuste_poisson,
            #ajuste_geometric,
            #ajuste_nbinom
),
legendtext = c("Exponential",
               "Log Normal",
               "Gamma",
               "uniform"
               #"Poisson",
               #"Geometric",
               #"Negative Binomial"
))

# Histograma con curvas de densidad ajustadas
denscomp(list(ajuste_exponential, 
              ajuste_lognormal,
              ajuste_gamma,
              ajuste_uniforme
              #ajuste_poisson,
              #ajuste_geometric,
              #ajuste_nbinom
),
legendtext = c("Exponential",
               "Log Normal",
               "Gamma",
               "uniform"
               #"Poisson",
               #"Geometric",
               #"Negative Binomial"
))

# Comparar AIC y BIC
aic_values <- c(ajuste_exponential$aic, 
                ajuste_lognormal$aic,
                ajuste_gamma$aic,
                ajuste_uniforme$aic)

bic_values <- c(ajuste_exponential$bic, 
                ajuste_lognormal$bic,
                ajuste_gamma$bic,
                ajuste_uniforme$bic)

# Mostrar AIC y BIC
data.frame(Distribucion = c("Exponential","Log Normal","Gamma","Uniform"),
           AIC = aic_values, BIC = bic_values) # La distribución uniforme presenta el mejor ajuste

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

# Crear el histograma para visualizar los datos
hist(age_acrop,
     breaks = CI, 
     main = paste("Yucatan,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Ajustar el parámetro de la tasa (rate) de la distribución exponencial
#rate <- 1 / mean(age_acrop)  # Tasa estimada asumiendo un aexponencial

# Definir los intervalos de clase
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Yucatan,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calcular las frecuencias observadas
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Función para identificar las secuencias de ceros y su longitud
gaps <- function(vector) {
  # Inicializar variables
  secuencias <- data.frame(inicio = integer(), longitud = integer())
  longitud_actual <- 0
  inicio_secuencia <- NULL
  
  # Recorrer el vector
  for (i in 1:length(vector)) {
    if (vector[i] == 0) {
      # Si encontramos un cero, incrementamos la longitud actual
      if (is.null(inicio_secuencia)) {
        inicio_secuencia <- i
      }
      longitud_actual <- longitud_actual + 1
    } else {
      # Si encontramos un valor diferente de cero y hay una secuencia en curso
      if (!is.null(inicio_secuencia)) {
        secuencias <- rbind(secuencias, data.frame(inicio = inicio_secuencia, longitud = longitud_actual))
        # Reiniciar variables
        longitud_actual <- 0
        inicio_secuencia <- NULL
      }
    }
  }
  
  # Agregar la última secuencia si termina en cero
  if (!is.null(inicio_secuencia)) {
    secuencias <- rbind(secuencias, data.frame(inicio = inicio_secuencia, longitud = longitud_actual))
  }
  
  secuencias$size <- secuencias$longitud * (age_range/length(observed))
  
  
  return(secuencias)
}

# Número de simulaciones
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  # Generamos la simulación asumiendo una distribución uniforme
  simulated_data <- runif(length(age_acrop), min = min(age_acrop), max = max(age_acrop))
  
  # Calcular las frecuencias observadas en los mismos intervalos de clase
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  # Evaluar si hay gaps (frecuencia 0)
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis
gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

temp1 <- gap_size|>
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

# Crear el histograma para visualizar los datos
hist(age_acrop,
     breaks = CI, 
     main = paste("Barbados,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Definir los intervalos de clase
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Barbados,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calcular las frecuencias observadas
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Número de simulaciones
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  # Generamos la simulación asumiendo una distribución uniforme
  simulated_data <- runif(81, min = min(age_acrop), max = max(age_acrop))
  
  # Calcular las frecuencias observadas en los mismos intervalos de clase
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  # Evaluar si hay gaps (frecuencia 0)
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis

gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

temp2 <- gap_size|>
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

# Crear el histograma para visualizar los datos
hist(age_acrop,
     breaks = CI, 
     main = paste("Barbados,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Definir los intervalos de clase
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Barbados,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calcular las frecuencias observadas
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Número de simulaciones
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  # Generamos la simulación asumiendo una distribución uniforme
  simulated_data <- runif(69, min = min(age_acrop), max = max(age_acrop))
  
  # Calcular las frecuencias observadas en los mismos intervalos de clase
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  # Evaluar si hay gaps (frecuencia 0)
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis

gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

temp3 <- gap_size|>
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

# Crear el histograma para visualizar los datos
hist(age_acrop,
     breaks = CI, 
     main = paste("Saint Croix,", intervals, "years intervals"),
     xlab = "Age (ka)",
     ylab = "Frequency")

# Definir los intervalos de clase
breaks <- hist(age_acrop, 
               breaks = CI, 
               main = paste("Saint Croix,", intervals, "years intervals"),
               xlab = "Age (ka)",
               ylab = "Frequency")$breaks

# Calcular las frecuencias observadas
observed <- hist(age_acrop, breaks = breaks, plot = FALSE)$counts

# Número de simulaciones
n_sim <- 10000
simulated_gaps <- numeric(n_sim)
simulated_gaps_size <- vector(mode = "list", length = n_sim)

for (i in 1:n_sim) {
  # Generamos la simulación asumiendo una distribución uniforme
  simulated_data <- runif(52, min = min(age_acrop), max = max(age_acrop))
  
  # Calcular las frecuencias observadas en los mismos intervalos de clase
  sim_observed <- hist(simulated_data, breaks = CI, plot = FALSE)$counts
  
  # Evaluar si hay gaps (frecuencia 0)
  simulated_gaps[i] <- sum(sim_observed == 0)
  simulated_gaps_size[[i]] <- gaps(sim_observed)
}

# Gap size analysis
gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

temp4 <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))|>
  mutate(site = "Saint Croix")

bind_rows(temp1, temp2, temp3, temp4) |>
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
  #ggtitle("Barbados")+
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


bind_rows(temp1, temp2, temp3, temp4) |>
  write_csv(file = "Table_gap_size_pvalues.csv")
