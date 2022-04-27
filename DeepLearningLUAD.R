library(keras)
library(tfdatasets)
library(readr)
library(tidyverse)
library(rsample)
library(tensorflow)

datasetLUADLit <- read_csv("../data/LUADGenes.csv")
datasetLUAD <- read_csv("../data/GENESLUADSELEC.csv")

datasetLUADLit <- subset(datasetLUADLit, select= c(TP53,NF1,ARID1A,RB1,CDKN2A,PIK3CA,KRAS,KEAP1,EGFR,STK11,SMARCA4,RBM10,BRAF,ERBB2,SETD2,MGA,MET,ATM,U2AF1L4,RIT1,SMAD4,CTNNB1,APC,RAF1,NRAS,MAP2K1,vital_status) )
datasetLUAD34 <- subset(datasetLUAD, select= c(SI,FAT4,CDH10,PTPRD,RP1L1,ADGRG4,DNAH9,PAPPA2,TNR,KEAP1,DAMTS12,RYR3,APOB,MUC17,PCDH15,PCLO,ANK2,CSMD1,ZNF536,COL11A1,NAV3,FAT3,SPTA1,FLG,XIRP2,KRAS,ZFHX4,USH2A,LRP1B,RYR2,CSMD3,MUC16,TP53,TTN,vital_status) )
datasetLUAD14 <-  subset(datasetLUAD, select= c(TP53,MUC16,LRP1B,KRAS,MUC17,KEAP1,ADGRG4,PTPRD,CDH10,PCLO,FLG,PAPPA2,PCDH15,XIRP2,vital_status) )

# first we split between training and testing sets
write_csv(datasetLUADLit, file = "luad_train.csv")
write_csv(validacaoLUSC15, file = "luad_test.csv")
View(dataset1)

# datasetLUADLit
# first we split between training and testing sets
split <- initial_split(datasetLUADLit, prop = 4/5)
train <- training(split)
write_csv(train, file = "luad_train.csv")
test <- testing(split)
write_csv(test, file = "luad_test.csv")


TRAIN_DATA_URL<- "C:/Users/debor/Documents/UFRN/Dissertacao/LUAD/luad_train.csv"
TEST_DATA_URL<- "C:/Users/debor/Documents/UFRN/Dissertacao/LUAD/luad_test.csv"

train_file_path <- get_file("luad_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("luad_test.csv",TEST_DATA_URL )

train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1
)

train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

spec <- feature_spec(train_dataset, vital_status ~ .)

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)
layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()

model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 30, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()
predict(model, batch)

plot(history)

##prediction
test_predictions <- model %>% predict(test, select=c(vital_status))
test_predictions[ , 1]

real <- subset(test, select= c(vital_status)) 
real[,1]
View(pred_dataset )
pred_dataset = round(test_predictions)

#Display table with all predicted and actual values
final_statusLUADLit=cbind (real, pred_dataset )
colnames(final_statusLUADLit) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusLUADLit$Predição, final_statusLUADLit$Real))
print(cfm)

# datasetLUAD34
# first we split between training and testing sets
split <- initial_split(datasetLUAD34, prop = 4/5)
train <- training(split)
write_csv(train, file = "luad_train.csv")
test <- testing(split)
write_csv(test, file = "luad_test.csv")


TRAIN_DATA_URL<- "C:/Users/debor/Documents/UFRN/Dissertacao/LUAD/luad_train.csv"
TEST_DATA_URL<- "C:/Users/debor/Documents/UFRN/Dissertacao/LUAD/luad_test.csv"

train_file_path <- get_file("luad_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("luad_test.csv",TEST_DATA_URL )

train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1
)

train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

spec <- feature_spec(train_dataset, vital_status ~ .)

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)
layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()

model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 100, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()
predict(model, batch)

plot(history)

##prediction
test_predictions <- model %>% predict(test, select=c(vital_status))
test_predictions[ , 1]

real <- subset(test, select= c(vital_status)) 
real[,1]
View(pred_dataset )
pred_dataset = round(test_predictions)

final_statusLUAD34=cbind (real, pred_dataset )
colnames(final_statusLUAD34) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusLUAD34$Predição, final_statusLUAD34$Real))
print(cfm)

#datasetLUAD14
# first we split between training and testing sets
split <- initial_split(datasetLUAD14, prop = 4/5)
train <- training(split)
write_csv(train, file = "luad_train.csv")
test <- testing(split)
write_csv(test, file = "luad_test.csv")


TRAIN_DATA_URL<- "C:/Users/debor/Documents/UFRN/Dissertacao/LUAD/luad_train.csv"
TEST_DATA_URL<- "C:/Users/debor/Documents/UFRN/Dissertacao/LUAD/luad_test.csv"

train_file_path <- get_file("luad_train.csv",TRAIN_DATA_URL )
test_file_path <- get_file("luad_test.csv",TEST_DATA_URL )

train_dataset <- make_csv_dataset(
  train_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1,
)

test_dataset <- train_dataset <- make_csv_dataset(
  test_file_path, 
  field_delim = ",",
  batch_size = 5, 
  num_epochs = 1
)

train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()

spec <- feature_spec(train_dataset, vital_status ~ .)

spec <- feature_spec(train_dataset, vital_status ~ .) %>% 
  step_numeric_column(all_numeric(), normalizer_fn = scaler_standard()) %>% 
  step_categorical_column_with_vocabulary_list(all_nominal()) %>% 
  step_indicator_column(all_nominal())

spec <- fit(spec)
layer <- layer_dense_features(feature_columns = dense_features(spec))
train_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  layer()

model <- keras_model_sequential() %>% 
  layer_dense_features(feature_columns = dense_features(spec)) %>%
  layer_dense(units = 100, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 50, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 20, activation = "relu") %>% 
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = 1, activation = "sigmoid")

model %>% compile(
  loss = "binary_crossentropy",
  optimizer = "adam",
  metrics = "accuracy"
)

history <- model %>% 
  fit(
    train_dataset %>% dataset_use_spec(spec) %>% dataset_shuffle(500),
    epochs = 100,
    validation_data = test_dataset %>% dataset_use_spec(spec),
    verbose = 2,
  )
summary(model)

model %>% evaluate(test_dataset %>% dataset_use_spec(spec), verbose = 0)

batch <- test_dataset %>% 
  reticulate::as_iterator() %>% 
  reticulate::iter_next() %>% 
  reticulate::py_to_r()
predict(model, batch)

plot(history)

##prediction
test_predictions <- model %>% predict(test, select=c(vital_status))
test_predictions[ , 1]

real <- subset(test, select= c(vital_status)) 
real[,1]
View(pred_dataset )
pred_dataset = round(test_predictions)

final_statusLUAD14=cbind (real, pred_dataset )
colnames(final_statusLUAD14) = c("Real","Predição")

library(caret)
cfm=caret::confusionMatrix(table(final_statusLUAD14$Predição, final_statusLUAD14$Real))
print(cfm)

# ROC curve Figure 4A
library(pROC)
rocobjLUAD1 <- plot.roc(final_statusLUADLit$Real, final_statusLUADLit$Predição,
                     percent=TRUE, col="#000000",
                    ci=TRUE)
legend("bottomright", legend=c("LUAD lit","LUAD 35","LUAD 14"), 
       col=c("#000000","#008600","#0000FF"), lwd=2)

plot.roc (final_statusLUAD35$Real, final_statusLUAD35$Predição, percent=TRUE, col="#008600",ci=TRUE)
plot.roc(final_statusLUAD14$Real, final_statusLUAD14$Predição, percent=TRUE, col="#0000FF" ,ci=TRUE)
