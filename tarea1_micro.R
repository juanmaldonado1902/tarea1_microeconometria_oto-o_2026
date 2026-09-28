# Tarea 1 - Microeconometria Aplicada (ITAM)
# Creador : Juan Pablo Maldonado
# Revisor:

library(tidyverse)
library(sandwich)      
library(lmtest)        
library(modelsummary)  
library(knitr)
library(here)

# Rutas relativas a la raiz del proyecto: el script corre igual con Rscript,
# con source() o linea por linea desde RStudio.
fig <- function(x) here("figuras", x)
tab <- function(x) here("tablas", x)
walk(here(c("figuras", "tablas")), dir.create, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 11))

dummies <- c("driveway", "recreation", "fullbase", "gasheat", "aircon", "prefer")

hp <- read_csv(here("hp.csv"), show_col_types = FALSE) %>%
  mutate(across(all_of(dummies), ~ as.integer(.x == "yes")),
         log_price = log(price), log_lotsize = log(lotsize))

tt <- function(x) map_chr(x, ~ str_glue("\\texttt{{{.x}}}"))

envolver <- function(cuerpo, caption, label, nota = NULL) {
  c("\\begin{table}[H]", "\\centering",
    str_glue("\\caption{<caption>}", .open = "<", .close = ">"),
    str_glue("\\label{tab:<label>}", .open = "<", .close = ">"), cuerpo,
    if (!is.null(nota)) c("\\vspace{0.4em}", "\\begin{minipage}{0.96\\textwidth}",
                          paste("\\footnotesize", nota), "\\end{minipage}"),
    "\\end{table}")
}

# ---- Pregunta 3: estadistica descriptiva ------------------------------------
vars_desc <- c("price", "lotsize", "bedrooms", "bathrooms", "stories", "garage", dummies)

desc <- hp %>%
  select(all_of(vars_desc)) %>%
  pivot_longer(everything(), names_to = "var") %>%
  summarise(N = n(), Media = mean(value), Varianza = var(value), `Desv. est.` = sd(value),
            `Mínimo` = min(value), `Máximo` = max(value), .by = var) %>%
  arrange(match(var, vars_desc)) %>%
  mutate(Variable = str_glue("\\texttt{{{var}}}{ifelse(var %in% dummies, '$^{*}$', '')}"),
         .keep = "unused", .before = 1)

print(desc, n = Inf)

writeLines(envolver(
  kable(desc, "latex", booktabs = TRUE, escape = FALSE, align = "lrrrrrr", linesep = "",
        digits = 3, format.args = list(big.mark = ",", scientific = FALSE)),
  "Estadística descriptiva de las variables (546 viviendas, Windsor, Canadá).", "descriptivas",
  "$^{*}$ Variable dicotómica codificada 1 = sí, 0 = no."), tab("tabla1_descriptivas.tex"))

# ---- Pregunta 4: distribucion de price y log(price) -------------------------
asimetria <- function(x) mean((x - mean(x))^3) / mean((x - mean(x))^2)^1.5
atipicos  <- function(x) {
  lim <- quantile(x, c(.25, .75)) + c(-1.5, 1.5) * IQR(x)
  x[x < lim[1] | x > lim[2]]
}

map(list(price = hp$price, `log(price)` = hp$log_price),
    ~ tibble(N = length(.x), Media = mean(.x), Mediana = median(.x), DE = sd(.x),
             Asimetria = asimetria(.x), Min = min(.x), Max = max(.x),
             Tukey_sup = quantile(.x, .75) + 1.5 * IQR(.x), Atipicos = length(atipicos(.x)),
             Atip_min = min(atipicos(.x)), Atip_3sd = sum(.x > mean(.x) + 3 * sd(.x)),
             `(Max-Media)/DE` = (max(.x) - mean(.x)) / sd(.x))) %>%
  bind_rows(.id = "Variable") %>% as.data.frame() %>% print(digits = 8)

hp %>% slice_max(price) %>% select(rownames, price, lotsize, bedrooms, bathrooms, prefer) %>% print()

histograma <- function(x, etiqueta, archivo) {
  ggsave(archivo, width = 6, height = 4,
         ggplot(hp, aes({{ x }})) +
           geom_histogram(bins = 30, fill = "grey70", colour = "white") +
           geom_vline(aes(xintercept = mean({{ x }})), linetype = "dashed", colour = "firebrick") +
           scale_x_continuous(labels = scales::comma) +
           labs(title = str_glue("Distribución de {etiqueta}"),
                subtitle = "La línea punteada marca la media", x = etiqueta, y = "Frecuencia"))
}

histograma(price,     "el precio de venta",      fig("fig1_hist_price.pdf"))
histograma(log_price, "el logaritmo del precio", fig("fig2_hist_log_price.pdf"))

# ---- Pregunta 5: price contra lotsize ---------------------------------------
# Los R2 solo son comparables entre modelos con la misma variable dependiente.
tribble(~forma,        ~f,
        "nivel-nivel", price     ~ lotsize,
        "log-nivel",   log_price ~ lotsize,
        "nivel-log",   price     ~ log_lotsize,
        "log-log",     log_price ~ log_lotsize) %>%
  mutate(ajuste = map(f, lm, data = hp),
         pendiente = map_dbl(ajuste, ~ coef(.x)[2]),
         R2 = map_dbl(ajuste, ~ summary(.x)$r.squared), .keep = "unused") %>%
  print()

print(c(correlacion_price_lotsize = cor(hp$price, hp$lotsize)))

base_scatter <- ggplot(hp, aes(lotsize, price)) +
  geom_point(alpha = .35) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Tamaño del lote (pies cuadrados)", y = "Precio (dólares canadienses)")

ggsave(fig("fig4_scatter_lm.pdf"), width = 6, height = 4,
       base_scatter + geom_smooth(method = "lm", formula = y ~ x, colour = "firebrick") +
         labs(title = "Precio y tamaño del lote", subtitle = "Con la recta de regresión lineal simple"))

# ---- Pregunta 6: regresion no parametrica (puntos extra) --------------------
# Nadaraya-Watson: promedio local ponderado por un kernel gaussiano.
bws <- c(1000, 2500, 5000)
nw  <- function(h, x) approx(ksmooth(hp$lotsize, hp$price, "normal", bandwidth = h,
                                     x.points = hp$lotsize), xout = x, ties = mean)$y

m_lin   <- lm(price ~ lotsize, data = hp)
m_ll    <- lm(log_price ~ log_lotsize, data = hp)
duan_ll <- mean(exp(resid(m_ll)))            # retransformacion de Duan a niveles
f_ll    <- function(x) exp(predict(m_ll, tibble(log_lotsize = log(x)))) * duan_ll

r2 <- function(f) 1 - sum((hp$price - f)^2) / sum((hp$price - mean(hp$price))^2)
c(list(`MCO lineal` = fitted(m_lin), `Log-log retransformado` = f_ll(hp$lotsize)),
  set_names(map(bws, nw, x = hp$lotsize), str_glue("Nadaraya-Watson h={bws}"))) %>%
  imap(~ tibble(Ajuste = .y, R2 = r2(.x))) %>% bind_rows() %>% as.data.frame() %>% print(digits = 8)

rejilla <- tibble(lotsize = seq(min(hp$lotsize), max(hp$lotsize), length.out = 300))

ggsave(fig("fig7_bandwidths.pdf"), width = 6.5, height = 4.2,
       base_scatter +
         geom_line(aes(colour = h), linewidth = .8, data = map(bws, ~ tibble(
           lotsize = rejilla$lotsize, h = str_glue("h = {format(.x, big.mark = ',')}"),
           price = nw(.x, rejilla$lotsize))) %>% bind_rows() %>% drop_na()) +
         scale_colour_manual("Ancho de banda", values = c("firebrick", "darkgreen", "steelblue")) +
         labs(title = "Nadaraya-Watson con tres anchos de banda",
              subtitle = "A menor h, menor sesgo y mayor varianza") +
         theme(legend.position = "bottom"))

# ---- Pregunta 7: las cinco especificaciones de MCO (HC1) --------------------
modelos <- list(
  price     ~ lotsize,
  log_price ~ lotsize,
  price     ~ log_lotsize + bedrooms + bathrooms,
  log_price ~ log_lotsize + bedrooms + bathrooms + stories,
  log_price ~ log_lotsize + bedrooms + bathrooms + stories + driveway + aircon + garage + prefer) %>%
  map(lm, data = hp)

rob <- map(modelos, ~ coeftest(.x, vcov. = vcovHC(.x, "HC1")))
walk(rob, print)

# Decimales por magnitud: 3 si pasa de mil, 7 si es muy chico (lotsize, que con
# el redondeo usual apareceria como 0.0001), 4 en el resto.
fmt_mag <- function(x) map_chr(x, function(v) if (is.na(v)) NA_character_ else
  formatC(v, format = "f", big.mark = ",",
          digits = if (abs(v) >= 1000) 3 else if (abs(v) < 1e-4) 7 else 4))

gof <- tribble(~raw, ~clean, ~fmt, "nobs", "Observaciones", 0, "r.squared", "$R^2$", 4)

etiquetas <- c(lotsize = "lotsize", log_lotsize = "log(lotsize)", bedrooms = "bedrooms",
               bathrooms = "bathrooms", stories = "stories", driveway = "driveway",
               aircon = "aircon", garage = "garage", prefer = "prefer")

# modelsummary hace la estimacion y el formato; aqui solo se arma el LaTeX booktabs.
tabla_mco <- function(mods, etiq, encabezado, caption, label, archivo, nota = NULL) {
  d <- msummary(mods, output = "data.frame", vcov = "HC1", fmt = fmt_mag, gof_map = gof,
                coef_map = c(tt(etiq), `(Intercept)` = "Constante"),
                stars = c("*" = .1, "**" = .05, "***" = .01)) %>%
    # raya donde la variable no entra; el resto, a notacion LaTeX
    mutate(across(-c(part, term, statistic), ~ if_else(.x == "",
             if_else(statistic == "estimate", "---", ""),
             .x %>% str_replace("^-", "$-$") %>% str_replace("(\\*+)$", "$^{\\1}$"))),
           term = if_else(statistic == "std.error", "", term))

  filas <- d %>% select(-part, -statistic) %>% pmap_chr(~ paste(c(...), collapse = " & ")) %>%
    paste0(if_else(d$statistic == "std.error" & d$part == "estimates", " \\\\[0.4em]", " \\\\"))

  writeLines(envolver(c(
    "\\small", str_glue("\\begin{tabular}{l <strrep('c', length(mods))>}", .open = "<", .close = ">"),
    "\\toprule", encabezado, "\\midrule", filas[d$part == "estimates"],
    "\\midrule", filas[d$part == "gof"], "\\bottomrule", "\\end{tabular}"),
    caption, label, nota), archivo)
}

tabla_mco(modelos, etiquetas,
  c("& \\multicolumn{5}{c}{\\textit{Variable dependiente}} \\\\", "\\cmidrule(lr){2-6}",
    paste("&", paste(tt(c("price", "log(price)", "price", "log(price)", "log(price)")),
                     collapse = " & "), "\\\\"),
    paste("&", paste(sprintf("(%d)", 1:5), collapse = " & "), "\\\\")),
  "Estimaciones de MCO", "tabla2", tab("tabla2.tex"),
  paste("\\textit{Notas:} errores estándar robustos a heterocedasticidad (HC1) entre paréntesis.",
        "Los asteriscos indican significancia al $^{*}10\\%$, $^{**}5\\%$ y $^{***}1\\%$.",
        "Una raya indica que la variable no se incluye en esa especificación."))

# ---- Pregunta 10: la elasticidad del lote en zona preferida -----------------
m10   <- lm(log_price ~ log_lotsize + bedrooms + bathrooms + stories + driveway +
              aircon + garage + prefer + log_lotsize:prefer, data = hp)
V10   <- vcovHC(m10, "HC1")
rob10 <- coeftest(m10, vcov. = V10)
gl10  <- df.residual(m10)
print(rob10)

# Var(b1 + b3) = V11 + V33 + 2*V13
elasticidad <- function(comb) {
  est <- sum(coef(m10)[comb]); se <- sqrt(sum(V10[comb, comb]))
  tibble(Elasticidad = est, EE = se, t = est / se,
         li90 = est - qt(.95, gl10) * se, ls90 = est + qt(.95, gl10) * se)
}

bind_rows(`Fuera de zona preferida` = elasticidad("log_lotsize"),
          `En zona preferida` = elasticidad(c("log_lotsize", "log_lotsize:prefer")),
          .id = "Grupo") %>% as.data.frame() %>% print(digits = 8)

print(c(t = rob10["log_lotsize:prefer", 3], p = rob10["log_lotsize:prefer", 4],
        critico_5pct = qt(.975, gl10)))

tabla_mco(list(m10), c(etiquetas[-1], `log_lotsize:prefer` = "log(lotsize)$\\times$prefer"),
  paste("&", tt("log(price)"), "\\\\"),
  "Especificación (5) con interacción entre el tamaño del lote y la zona preferida.",
  "interaccion", tab("tabla4_interaccion.tex"))

# ---- Pregunta 11: prediccion con la especificacion (5) ----------------------
m5 <- modelos[[5]]
x0 <- tibble(lotsize = 6000, bedrooms = 3, bathrooms = 2, stories = 2,
             driveway = 1, aircon = 1, garage = 1, prefer = 1) %>%
  mutate(log_lotsize = log(lotsize))

xb     <- predict(m5, x0)
sigma2 <- summary(m5)$sigma^2
duan   <- mean(exp(resid(m5)))   # smearing: no supone normalidad

tibble(Concepto = c("x'b (log)", "sigma^2", "factor de Duan", "exp(sigma^2/2)",
                    "Ingenua exp(x'b)", "Normalidad", "Duan (reportada)"),
       Valor = c(xb, sigma2, duan, exp(sigma2 / 2),
                 exp(xb), exp(xb + sigma2 / 2), exp(xb) * duan)) %>%
  mutate(Valor = formatC(Valor, format = "f", digits = 4, big.mark = ",")) %>%
  as.data.frame() %>% print()

# Bootstrap no parametrico: se remuestrean viviendas completas y se reestima (5).
# La semilla va aqui y no al inicio, para que el resultado no dependa de lo que
# se haya corrido antes en el script.
set.seed(14303)
boot <- map(1:1000, function(i) {
  ajuste <- lm(formula(m5), data = slice_sample(hp, prop = 1, replace = TRUE))
  p <- unname(exp(predict(ajuste, x0)))
  c(duan = p * mean(exp(resid(ajuste))), ingenua = p)
}) %>% bind_rows()

ic_boot  <- quantile(boot$duan, c(.035, .965))
ic_boot0 <- quantile(boot$ingenua, c(.035, .965))

# IC analitico: el EE robusto va a mano porque predict() usa la matriz clasica
X0     <- model.matrix(delete.response(terms(m5)), x0)
V5     <- vcovHC(m5, "HC1")
se_rob <- sqrt(drop(X0 %*% V5 %*% t(X0)))
ic_log <- xb + c(-1, 1) * qt(.965, df.residual(m5)) * se_rob

tibble(Metodo = c("Bootstrap (Duan)", "Bootstrap (exp(x'b))", "Analitico (exp(x'b))"),
       li = c(ic_boot[1], ic_boot0[1], exp(ic_log[1])),
       ls = c(ic_boot[2], ic_boot0[2], exp(ic_log[2]))) %>%
  mutate(ancho = ls - li, centro = (ls + li) / 2) %>% as.data.frame() %>% print(digits = 8)

print(c(boot_media = mean(boot$duan), boot_DE = sd(boot$duan)))
print(c(ee_robusto = se_rob, ee_clasico = predict(m5, x0, se.fit = TRUE)$se.fit))

ggsave(fig("fig8_bootstrap.pdf"), width = 6.5, height = 4,
       ggplot(boot, aes(duan)) +
         geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey70", colour = "white") +
         geom_density(colour = "steelblue", linewidth = .7) +
         geom_vline(xintercept = exp(xb) * duan, colour = "firebrick", linewidth = .7) +
         geom_vline(xintercept = ic_boot, linetype = "dashed", colour = "firebrick") +
         scale_x_continuous(labels = scales::comma) +
         labs(title = "Distribución bootstrap de la predicción (1,000 réplicas)",
              subtitle = "Línea sólida: predicción de Duan; punteadas: percentiles 3.5 y 96.5",
              x = "Precio predicho (dólares canadienses)", y = "Densidad"))

# ---- Pregunta 12: lotsize y log(lotsize) en la misma regresion --------------
m12 <- lm(log_price ~ log_lotsize + lotsize + bedrooms + bathrooms + stories + driveway +
            aircon + garage + prefer, data = hp)
V12 <- vcovHC(m12, "HC1")
print(coeftest(m12, vcov. = V12))

print(c(coeficientes_NA = sum(is.na(coef(m12))), rango_X = qr(model.matrix(m12))$rank,
        columnas_X = ncol(model.matrix(m12))))
print(c(cor_lotsize_log = cor(hp$lotsize, hp$log_lotsize),
        ee_log_lotsize_en_5 = sqrt(V5["log_lotsize", "log_lotsize"]),
        ee_log_lotsize_en_12 = sqrt(V12["log_lotsize", "log_lotsize"])))

bind_rows(`(12)` = enframe(car::vif(m12)), `(5)` = enframe(car::vif(m5)), .id = "Modelo") %>%
  pivot_wider(names_from = Modelo) %>% as.data.frame() %>% print(digits = 6)

b12 <- coef(m12)
tibble(Punto = c("Mediana", "Media", "Máximo"),
       lotsize = c(median(hp$lotsize), mean(hp$lotsize), max(hp$lotsize))) %>%
  mutate(dlogp_dlot = b12["lotsize"] + b12["log_lotsize"] / lotsize,
         pct_1000_exacto = 100 * (exp(1000 * dlogp_dlot) - 1),
         dolares_1000 = (exp(1000 * dlogp_dlot) - 1) * median(hp$price),
         elasticidad = b12["log_lotsize"] + b12["lotsize"] * lotsize) %>%
  as.data.frame() %>% print(digits = 8)

tabla_mco(list(m12), tt(etiquetas[c(2, 1, 3:9)]) %>% set_names(names(etiquetas)[c(2, 1, 3:9)]),
  paste("&", tt("log(price)"), "\\\\"),
  "Especificación (5) con \\texttt{lotsize} y $\\log(\\texttt{lotsize})$ simultáneamente.",
  "lotsize_doble", tab("tabla5_lotsize_doble.tex"))
