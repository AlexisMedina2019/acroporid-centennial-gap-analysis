# análisis exploratorios

library(ggplot2)
library(dplyr)
library(tidyr)
library(janitor)
library(readxl)
library(purrr)


raw_ages <- read_excel("raw_ages.xlsx")
View(raw_ages)


source <- levels(factor(raw_ages$source))

age<-raw_ages|>
  filter(source == "This paper")|>
  mutate(t_ka2 = round(t_ka, 1))|>
  group_by(species_2)|>
  reframe(tabyl(t_ka2))

|>
  ggplot(aes(x = t_ka2, y = percent))+
  geom_col(aes(colour = species_2))+
  scale_y_continuous(breaks = seq(0,0.2, 0.01))+
  theme_bw()
  #geom_smooth(aes(colour = species_2))+
  #geom_smooth(method = "lm")+
  #geom_smooth()
  
  
# Establecer la semilla para reproducibilidad
set.seed(123)

# Generar 50 observaciones de una distribución exponencial
datos <- rexp(50, rate = 0.5)

# Crear el histograma para visualizar los datos
hist(datos, breaks = "Sturges", main = "Histograma de los datos", xlab = "Antigüedad (años)", ylab = "Frecuencia")

# Ajustar el parámetro de la tasa (rate) de la distribución exponencial
rate <- 1 / mean(datos)  # Tasa estimada

# Definir los intervalos de clase
breaks <- hist(datos, breaks = "Sturges", plot = FALSE)$breaks

# Calcular las frecuencias observadas
observed <- hist(datos, breaks = breaks, plot = FALSE)$counts

# Calcular las frecuencias esperadas bajo la distribución exponencial
expected <- diff(pexp(breaks, rate = rate)) * length(datos)

# Ajustar para evitar problemas con frecuencias muy pequeñas
expected <- ifelse(expected < 1, 1, expected)

# Realizar la prueba de chi-cuadrado
chisq_test <- chisq.test(observed, p = expected/sum(expected))

# Ver el resultado de la prueba
print(chisq_test)

# Calcular los residuos
residuals <- chisq_test$stdres

# Mostrar los residuos
plot(residuals)


## Datos reales
age<-raw_ages|>
  filter(source == "This paper")

edades <- age$t_ka
hist(edades, breaks = 110, main = "Histograma FD", xlab = "Antigüedad (años)", ylab = "Frecuencia")
breaks <- hist(edades, breaks = 110, main = "Histograma FD", xlab = "Antigüedad (años)", ylab = "Frecuencia")$breaks
# Calcular las frecuencias observadas
observed <- hist(s, breaks = breaks, plot = FALSE)$counts

#Datos esperados
rate <- 1 / mean(edades)
# Calcular las frecuencias esperadas bajo la distribución exponencial
expected <- diff(pexp(breaks, rate = rate)) * length(edades)

# Ajustar para evitar problemas con frecuencias muy pequeñas
expected <- ifelse(expected < 1, 1, expected)

# Realizar la prueba de chi-cuadrado
chisq_test <- chisq.test(observed, p = expected/sum(expected))




medias <- age|>
  group_by(species_2)|>
  summarise(mean(t_ka2))



##### identificación de distribución teórica
library(fitdistrplus)
library(MASS)
library(readxl)
library(dplyr)

raw_ages <- read_excel("raw_ages.xlsx")
age<-raw_ages|>
  filter(source == "This paper")

edades <- age$t_ka

# Ajustar varias distribuciones
ajuste_exponential <- fitdist(edades, "exp")
ajuste_lognormal <- fitdist(edades, "lnorm")
ajuste_gamma <- fitdist(edades, "gamma")
ajuste_uniforme <- fitdist(edades, "unif")
#ajuste_poisson <- fitdist(edades, "pois", method = "mme", discrete = TRUE)
#ajuste_geometric <- fitdist(edades, "geom", method = "mme")
#ajuste_nbinom <- fitdist(edades, "nbinom", method = "mse", discrete = TRUE)

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
           AIC = aic_values, BIC = bic_values)

# Bootstraping ------------------------------------------------------------
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
n_sim <- 9999
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

# Comparar con los gaps observados en los datos reales
observed_gaps <- sum(observed == 0)

# Proporción de simulaciones que tienen más o igual número de gaps que los datos reales
p_value_gaps <- mean(simulated_gaps >= observed_gaps)
p_value_gaps

# Gap size analysis

gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())
  
gap_size_p <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))
  
gap_size_p|>
  ggplot(aes(x = size, y = freq))+
  scale_y_continuous(breaks = seq(0,1,0.05))+
  scale_x_continuous(breaks = seq(50,1200,50))+
  geom_hline(yintercept = 0.05, colour = "red", linewidth = 1)+
  ylab("Probability under H0 of unifor distribution")+
  xlab("Gap size (years)")+
  geom_line()+
  theme_bw()


#### Gaps in Barbados data ####
raw_ages <- read_excel("raw_ages.xlsx")
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
n_sim <- 9999
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

# Comparar con los gaps observados en los datos reales
observed_gaps <- sum(observed == 0)

# Proporción de simulaciones que tienen más o igual número de gaps que los datos reales
p_value_gaps <- mean(simulated_gaps >= observed_gaps)
p_value_gaps

# Gap size analysis

gap_size <- list_rbind(simulated_gaps_size, names_to = "simulation")|>
  group_by(simulation, size)|>
  reframe(f = n())

gap_size_p <- gap_size|>
  group_by(size)|>
  summarise(average = round(mean(f),0), times = length(size), freq = (times+1)/(n_sim+1))

gap_size_p|>
  ggplot(aes(x = size, y = freq))+
  scale_y_continuous(breaks = seq(0,1,0.05))+
  scale_x_continuous(breaks = seq(50,1200,50))+
  geom_hline(yintercept = 0.05, colour = "red", linewidth = 1)+
  ylab("Probability under H0 of uniform distribution")+
  xlab("Gap size (years)")+
  ggtitle("Barbados")+
  geom_line()+
  theme_bw()



