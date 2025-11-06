
library(tidyverse)
library(here)
library(infer)
library(randomForest)
library(tidymodels)
library(modelr)
library(yardstick)

source(here("scripts", "preprocessing.R"))
load(here("data", "biomarker-clean.RData"))
glimpse(biomarker_clean)

test_fn <- function(.df){
  infer::t_test(.df,
                formula = level ~ group,
                order = c("ASD", "TD"),
                alternative = "two-sided",
                var.equal = FALSE
  )
}

ttests_out <- biomarker_clean %>%
  select(-ados) %>%
  pivot_longer(-group, names_to = "protein", values_to = "level") %>%
  nest(data = c(level, group)) %>%
  mutate(ttest = map(data, test_fn)) %>%
  unnest(ttest) %>%
  arrange(p_value) %>%
  mutate(
    m = n(),
    hm = log(m) + 1/(2*m) - digamma(1),
    rank = row_number(),
    p.adj = m * hm * p_value / rank
  )

proteins_s1 <- ttests_out %>% slice_min(p.adj, n = 10) %>% pull(protein)

predictors <- biomarker_clean %>% select(-c(group, ados))
response   <- biomarker_clean %>% pull(group) %>% factor()

set.seed(101422)
rf_out <- randomForest(
  x = predictors,
  y = response,
  ntree = 1000,
  importance = TRUE
)

print(rf_out$confusion)

proteins_s2 <- rf_out$importance %>%
  as_tibble() %>%
  mutate(protein = rownames(rf_out$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)

proteins_sstar <- intersect(proteins_s1, proteins_s2)

biomarker_sstar <- biomarker_clean %>%
  select(group, any_of(proteins_sstar)) %>%
  mutate(class = if_else(group == "ASD", "ASD", "TD")) %>%
  select(-group)

set.seed(101422)
biomarker_split <- initial_split(biomarker_sstar, prop = 0.8)
train_data <- training(biomarker_split)
test_data  <- testing(biomarker_split)

train_data <- train_data %>% mutate(class = factor(class, levels = c("TD", "ASD")))
test_data  <- test_data  %>% mutate(class = factor(class, levels = c("TD", "ASD")))

fit <- glm(class ~ ., data = train_data, family = binomial(link = "logit"))

test_data <- test_data %>%
  mutate(
    pred_prob  = as.numeric(predict(fit, newdata = test_data, type = "response")),
    pred_class = factor(if_else(pred_prob > 0.5, "ASD", "TD"), levels = c("TD", "ASD"))
  )

cls_metrics <- metric_set(sensitivity, specificity, accuracy)
cls_out <- cls_metrics(
  test_data,
  truth    = class,
  estimate = pred_class,
  event_level = "second"     
)

auc_out <- roc_auc(
  test_data,
  truth = class,
  pred_prob,
  event_level = "second"
)

metrics_out <- bind_rows(
  cls_out,
  tibble(.metric = "roc_auc", .estimator = "binary", .estimate = as.numeric(auc_out$.estimate))
)

print(metrics_out)

cat("\nProteins in intersection (used for model):\n")
print(proteins_sstar)
