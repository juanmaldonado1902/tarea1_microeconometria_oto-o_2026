# Tarea 1 - Microeconometria Aplicada (ITAM)
# Creado por : Juan Pablo Maldonado
# Revisado y editado por: Arturo Aguilar

# Preparacion
library(tidyverse)
library(sandwich)      
library(lmtest)        
library(stargazer)
library(car)
theme_set(theme_minimal(base_size = 11))

# Definir nombre de variables que se transformaran a dummies
dummies <- c("driveway", "recreation", "fullbase", "gasheat", "aircon", "prefer")

# Abrir la base d edatos y crear las dummies correspondientes
hp <- read_csv("hp.csv", show_col_types = FALSE) %>%
  mutate(across(all_of(dummies), ~ as.integer(.x == "yes")),
         log_price = log(price), log_lotsize = log(lotsize))

# ====/// 3: estadistica descriptiva \\\=====
vars_desc <- c("price", "lotsize", "bedrooms", "bathrooms", "stories", "garage", dummies)

desc <- hp %>% select(all_of(vars_desc)) %>%
  rename_with(~ paste0(.x, "*"), all_of(dummies)) %>% as.data.frame()

stargazer(desc, digits = 3, header = FALSE, table.placement = "H",
          summary.stat = c("n", "mean", "sd", "min", "max"),
          title = "Estadística descriptiva de las variables (546 viviendas, Windsor, Canadá).",
          label = "tab:descriptivas", notes.align = "l",
          notes = "$^{*}$ Variable dicotómica codificada 1 = sí, 0 = no.",
          out = "tablas/tabla1_p3.tex")

# ====/// 4: distribucion de price y log(price) \\\=====

# Histograma de price
ggplot(hp, aes(x = price)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "white") +
  scale_x_continuous(labels = scales::comma) +
  geom_vline(xintercept = mean(hp$price), linetype = "dashed") +
  labs(title = "Distribución del precio de venta", x = "Precio", y = "Frecuencia")
ggsave("figuras/fig1_p4_price.pdf", width = 6, height = 4)

# Histograma de log(price)
ggplot(hp, aes(x = log_price)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "white") +
  geom_vline(xintercept = mean(hp$log_price), linetype = "dashed") +
  labs(title = "Distribución del logaritmo del precio", x = "log(precio)", y = "Frecuencia")
ggsave("figuras/fig2_p4_logprice.pdf", width = 6, height = 4)

# Outliers usando el rango intercuantílico
hp_sin_atip <- hp %>%
  filter(price >= quantile(price, .25) - 1.5 * IQR(price),
         price <= quantile(price, .75) + 1.5 * IQR(price))

nrow(hp) - nrow(hp_sin_atip)   # numero de outliers eliminados

# Histograma de price sin outliers
ggplot(hp_sin_atip, aes(x = price)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "white") +
  geom_vline(xintercept = mean(hp_sin_atip$price), linetype = "dashed") +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Distribución del precio de venta sin outliers", x = "Precio", y = "Frecuencia")
ggsave("figuras/fig3_p4_price_outliers.pdf", width = 6, height = 4)


# ====/// 5: price contra lotsize \\\=====

ggplot(hp, aes(lotsize, price)) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point(alpha = .35) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Tamaño del lote (pies cuadrados)", y = "Precio (dólares canadienses)")

ggsave("figuras/fig4_p5_scatter.pdf", width = 6, height = 4)


# ====/// 6: Regresion no parametrica (puntos extra) \\\=====

hp$nw_1000 <- NA
hp$nw_2000 <- NA
hp$ll_2000 <- NA

for (i in 1:nrow(hp)) {
  dist   <- hp$lotsize - hp$lotsize[i]          # distancia de cada vivienda a la vivienda i
  w_1000 <- pmax(1 - abs(dist / 1000), 0)       # pesos kernel triangular, h = 1,000
  w_2000 <- pmax(1 - abs(dist / 2000), 0)       # pesos kernel triangular, h = 2,000
  
  # Nadaraya-Watson: promedio ponderado de price
  hp$nw_1000[i] <- sum(w_1000 * hp$price) / sum(w_1000)
  hp$nw_2000[i] <- sum(w_2000 * hp$price) / sum(w_2000)
  
  # Lineal local: MCO ponderado de price sobre la distancia; la constante es la prediccion
  hp$ll_2000[i] <- coef(lm(price ~ dist, data = hp, weights = w_2000))[1]
}

# Scatterplot con las tres estimaciones
ggplot(hp, aes(lotsize, price)) +
  geom_point(alpha = .35) +
  geom_line(aes(y = nw_1000, colour = "NW, h = 1,000")) +
  geom_line(aes(y = nw_2000, colour = "NW, h = 2,000")) +
  geom_line(aes(y = ll_2000, colour = "LLR, h = 2,000")) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Tamaño del lote (pies cuadrados)", y = "Precio (dólares canadienses)",
       colour = "Estimación") +
  theme(legend.position = "bottom")
ggsave("figuras/fig5_p6_kernel.pdf", width = 6, height = 4)

# ====/// 7: Tabla de MCO  \\\=====

m1 <- lm(price ~ lotsize, data = hp)
m2 <- lm(log_price ~ lotsize, data = hp)
m3 <- lm(price ~ log_lotsize + bedrooms + bathrooms, data = hp)
m4 <- lm(log_price ~ log_lotsize + bedrooms + bathrooms + stories, data = hp)
m5 <- lm(log_price ~ log_lotsize + bedrooms + bathrooms + stories + 
           driveway + aircon + garage + prefer, data = hp)

se1 <- sqrt(diag(vcovHC(m1,type="HC1")))
se2 <- sqrt(diag(vcovHC(m2,type="HC1")))
se3 <- sqrt(diag(vcovHC(m3,type="HC1")))
se4 <- sqrt(diag(vcovHC(m4,type="HC1")))
se5 <- sqrt(diag(vcovHC(m5,type="HC1")))

stargazer(m1, m2, m3, m4, m5, se = list(se1, se2, se3, se4, se5),
          header = FALSE, table.placement = "H", digits = 3,
          title = "Estimaciones de MCO", label = "tab:tabla2",
          dep.var.caption = "Variable dependiente",
          dep.var.labels = c("price", "log(price)", "price", "log(price)"),
          covariate.labels = c("lotsize", "log(lotsize)", "bedrooms", "bathrooms", "stories",
                               "driveway", "aircon", "garage", "prefer", "Constante"),
          keep.stat = "rsq",
          add.lines = list(c("n", nobs(m1), nobs(m2), nobs(m3), nobs(m4), nobs(m5))),
          notes.label = "Notas:", notes.append = FALSE, notes.align = "l",
          notes = c("Errores estándar robustos (HC1) entre paréntesis.",
                    "$^{*}$ significativo al $10%$; $^{**}$ al $5%$; $^{***}$ al $1%$"),
          out = "tablas/tabla2_p7.tex")

# ====/// 10: Interaccion con zona preferida  \\\=====

m10   <- lm(log_price ~ log_lotsize + bedrooms + bathrooms + stories + driveway +
              aircon + garage + prefer + log_lotsize:prefer, data = hp)
se10 <- sqrt(diag(vcovHC(m10,type="HC1")))
stargazer(m10, se = list(se10),
          header = FALSE, table.placement = "H", digits = 3,
          title = "Heterogeneidad para prefer", label = "tab:tabla10",
          dep.var.caption = "Variable dependiente",
          dep.var.labels = c("log(price)"),
          covariate.labels = c("log(lotsize)", "bedrooms", "bathrooms", "stories",
                               "driveway", "aircon", "garage", "prefer", 
                               "log(lotsize) x Prefer","Constante"),
          keep.stat = "rsq",
          add.lines = list(c("n", nobs(m10))),
          notes.label = "Notas:", notes.append = FALSE, notes.align = "l",
          notes = c("Errores estándar robustos (HC1) entre paréntesis.",
                    "$^{*}$ significativo al $10%$; $^{**}$ al $5%$; $^{***}$ al $1%$"),
          out = "tablas/tabla_p10.tex")

# Elasticidades
V10   <- vcovHC(m10, "HC1")

(theta_NP <- coef(m10)['log_lotsize'])
(var_NP <- vcovHC(m10, "HC1")['log_lotsize','log_lotsize'])
(IC_NP <- c(theta_NP-abs(qnorm(0.05))*sqrt(var_NP),
           theta_NP+abs(qnorm(0.05))*sqrt(var_NP)))

(theta_P <- coef(m10)['log_lotsize']+coef(m10)['log_lotsize:prefer'])
(var_P <- vcovHC(m10, "HC1")['log_lotsize','log_lotsize']+
    vcovHC(m10, "HC1")['log_lotsize:prefer','log_lotsize:prefer']+
    2*vcovHC(m10, "HC1")['log_lotsize','log_lotsize:prefer'])
(IC_P <- c(theta_P-abs(qnorm(0.05))*sqrt(var_P),
           theta_P+abs(qnorm(0.05))*sqrt(var_P)))


# ====/// 11: prediccion, bootstrap y metodo delta  \\\=====

# a. Prediccion
beta5 <- coefficients(m5)
x_pred <- c(1,log(6000),3,2,2,1,1,1,1)
(price_pred <-exp(t(beta5)%*%x_pred))

#b. bootstrap
set.seed(14303)
price_boot <- c()
n_boot <- length(hp$price)

for (n in 1:1000) {
  muestra <- hp[sample(nrow(hp), size = n_boot, replace = TRUE), ]
  m5_boot <- lm(log_price ~ log_lotsize + bedrooms + bathrooms + stories + 
                        driveway + aircon + garage + prefer, data = muestra)
  beta5_boot <- coefficients(m5_boot)
  price_boot <- c(price_boot,exp(t(beta5_boot)%*%x_pred))
}

(IC_boot <- quantile(price_boot, c(.035, .965)))
price_boot <- as.data.frame(price_boot)

ggplot(price_boot, aes(price_boot)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey70", colour = "white") +
  geom_density(colour = "steelblue", linewidth = .7) +
  geom_vline(xintercept = mean(price_boot$price_boot), colour = "firebrick", linewidth = .7) +
  geom_vline(xintercept = IC_boot, linetype = "dashed", colour = "firebrick") +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Distribución bootstrap de la predicción (1,000 réplicas)",
       subtitle = "Línea sólida: predicción exp(x'b); punteadas: percentiles 3.5 y 96.5",
       x = "Precio predicho (dólares canadienses)", y = "Densidad")

ggsave("figuras/fig8_bootstrap.pdf", width = 6, height = 4)

# c. metodo delta

dg <- price_pred*x_pred
s2_md <- t(dg)%*%vcovHC(m5,type="HC1")%*%dg
(IC_md <- c(price_pred-abs(qnorm(0.035))*sqrt(s2_md),
            price_pred+abs(qnorm(0.035))*sqrt(s2_md)))
  
# ====/// 12: lotsize y log(lotsize) en la misma regresion  \\\=====

m12 <- lm(log_price ~ log_lotsize + lotsize + bedrooms + bathrooms + stories + driveway +
            aircon + garage + prefer, data = hp)
se12 <- sqrt(diag(vcovHC(m12,type="HC1")))

stargazer(m12, se = list(se12),
          header = FALSE, table.placement = "H", digits = 3,
          title = "lotsize y log(lotsize)", label = "tab:tabla12",
          dep.var.caption = "Variable dependiente",
          dep.var.labels = c("log(price)"),
          covariate.labels = c("log(lotsize)","lotsize", "bedrooms", "bathrooms", 
                               "stories","driveway", "aircon", "garage", "prefer", 
                               "Constante"),
          keep.stat = "rsq",
          add.lines = list(c("n", nobs(m12))),
          notes.label = "Notas:", notes.append = FALSE, notes.align = "l",
          notes = c("Errores estándar robustos (HC1) entre paréntesis.",
                    "$^{*}$ significativo al $10%$; $^{**}$ al $5%$; $^{***}$ al $1%$"),
          out = "tablas/tabla_p12.tex")

(avg_ch <- coefficients(m12)['log_lotsize']+coefficients(m12)['lotsize']*mean(hp$lotsize))
