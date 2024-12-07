library(dplyr)
library(plyr)

# STEP ZERO - INITILIZIATION
# caricare tutti i dati

#CLINICAL DATA
clinical <- read.table("./Dati_Cilia/brca_clinical.txt", header=TRUE, sep ="\t")
#DATI DELLA MATRICE DI ESPRESSIONE DEI GENI
exprs <- as.matrix(read.table("./Dati_Cilia/BRCA_expressions.txt", header=TRUE, sep = "\t"))
#substep 1.0 estraggo dalla matrice solo gli indici delle righe che corrispondono ai geni
exprs_rows <- exprs[,0] # matrice senza colonne solo righe
# nomi delle colonne per estrarre i nomi dei pazienti (servira' dopo)
c_matrix <- colnames(exprs)
# NB. sono nomi quindi saranno in formato chr
g_matrix <- as.integer(rownames(exprs))
# per estrarre e quindi visualizzare
# il nome di ciascuna riga e quindi il gene

#STEP ONE - creare la pathway per ciascun paziente
nodes <- read.table("./Dati_Cilia/metapathway_nodes.txt", header=TRUE, sep ="\t")
#substep 1.1 estraggo dai nodi solo gli id
g_nodes <- nodes[,1]
edges <- read.table("./Dati_Cilia/metapathway_edges.txt", header=TRUE, sep ="\t")

# substep 1.2 dopo aver estratto id-nodes e id-matrix-exprs faccio l'intersezione
# tra i due vettori in modo tale da individuare solo gli elementi in comune
# cio' che e' presente nei nodi della metapathway ma non e' presente nella matrice
# di espressione non verra' preso in considerazione
g_match <- intersect(g_matrix, g_nodes)

# substep 1.3 la matrice di espressione filtrata a cui fare riferimento e' la seguente:
filtered_exprs <- exprs[g_matrix%in%g_match,]
# vengono rimossi i nodi che non sono presenti nella matrice di espressione
filtered_nodes <- nodes[g_nodes%in%g_match,]

# substep 1.4 nella matrice di espressione devono essere rimossi per ciascun paziente
# i geni che sono non significativamente espressi compresi tra i valori -1 e 1 SOSTITUISCO CON NA

# substep 1.5 NA per i geni non significativamente sovra o sotto espressi
filtered_exprs[(filtered_exprs > -1) & (filtered_exprs < 1)] <- NA

# substep 1.6 considero gli archi i cui nodi di partenza/destinazione sono 
# effettivamente presenti all'interno di filtered_nodes 
edges_to_keep1 <- c(which(edges$Start %in% filtered_nodes$Id))
edges_to_keep2 <- c(which(edges$End %in% filtered_nodes$Id))

edges_to_keep <- intersect(edges_to_keep1, edges_to_keep2)

edges <- edges[edges_to_keep, ]

# idea: costruire le reti per OGNI paziente
# fare un ciclo che visualizza tutte le colonne con le varie espressioni
# dei geni che vengono filtrate come sovra-espressi/sotto-espressi
# quindi prendere il nome del paziente (colonna j-esima fissata ad ogni
# iterazione) e i valori di i per i che varia dal primo all'ultimo valore
# della riga e tenere conto anche del nome della riga che corrisponde al gene
# in modo tale che per ogni paziente ho come informazione l'id del gene
# e faccio l'intersezione con filtered_nodes e dall'insieme edge 
# estraggo per ogni paziente solo gli archi dei nodi effettivamente presenti
# per quel paziente. questo dovrebbe essere salvato in output 
# un file per ciascun paziente, praticamente una rete in formato .graph per ogni
# paziente

m <- ncol(filtered_exprs) # numero colonne
n <- nrow(filtered_exprs) # numero righe

# df temporaneo per visualizzare espressione del gene e il corrispettivo id
tmp_patient_nodes <- data.frame(matrix(ncol=2, nrow=0, dimnames = list(NULL, c("Espressione", "Id"))))
# df temporaneo contenente gli archi e quindi le connessioni tra i geni presenti per quel paziente
tmp_edges <- edges
setwd("/media/denise-pc/new")


b <- length(c_matrix)
for(q in 1:b)
  c_matrix[q] <- gsub("\\.", "-", c_matrix[q])

for(j in 1:1){ # ciclo sulle colonne, ciascun paziente viene preso in esame
  tmp_patient_nodes <- data.frame(matrix(ncol=2, nrow=0, dimnames = list(NULL, c("Espressione", "Id"))))
  tmp_edges <- edges
  for(i in 1:n){ # ciclo sulle righe, esaminiamo il valore di ogni singolo gene presente nella matrice di exprs
    if(is.na(filtered_exprs[i,j])) next # se incontriamo un NA che indica un gene non significativamente sovra o sotto espresso continuiamo nella prossima iterazione 
    
    # se il gene è sottoespresso aggiungiamo nel df una riga con il valore sottoespresso e l'id del gene
    if(filtered_exprs[i,j] <= -1){
      new_row <- c("sottoespresso", g_match[i])
      tmp_patient_nodes[nrow(tmp_patient_nodes) + 1, ] <- new_row
    }
    # viceversa lo stesso se è sovraespresso
    if(filtered_exprs[i,j] >= 1){
      new_row <- c("sovraespresso", g_match[i])
      tmp_patient_nodes[nrow(tmp_patient_nodes) + 1, ] <- new_row
    }
  }
  tot_tmp_nodes <- as.character(nrow(tmp_patient_nodes)) # cast a chr per utilizzare la fne writeLines
  
  # condizioni per tenere traccia degli archi che effettivamente collegano i geni presenti per quel determinato
  # paziente in esame 
  edges_to_keep1 <- c(which((tmp_edges$Start %in% tmp_patient_nodes$Id) & (tmp_edges$End %in% tmp_patient_nodes$Id)))
  tmp1_edges <- tmp_edges[edges_to_keep1, ]
  
  # output che segue il formato network da dare in input a MultiMotif
  name <- paste(c_matrix[j], ".txt", sep = "")
  fileP <- file(name)
  writeLines(c("directed", tot_tmp_nodes), fileP)
  write.table(tmp_patient_nodes[,1], name, append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE)
  # substep 1.7 sappiamo che la numerazione all'interno dell'algoritmo MultiMotif va da 0 a N-1
  # dove 0 indica il primo nodo ed N-1 l'ultimo nodo presente, pertanto bisogna cambiare l'identificazione
  # dei nodi, che nel caso del df tmp1_edges utilizza come indici gli Id dei geni per riferirsi ai nodi.
  # Per farlo individuiamo la posizione di ciascun gene all'interno del df tmp_patient_nodes e sottraiamo 1
  # Esempio: il gene identificato con il numero 2 è il primo gene della lista, dunque il primo nodo;
  # L'algoritmo MultiMotif lo vedrà dunque come nodo 0 di partenza. Dai successivi confronti 
  # vediamo immediatamente che l'Id del gene 2 è nella prima posizione. 1-1 = 0 
  # Per quanto riguarda il corrispondente nodo destinazione 5664 sappiamo invece che si trova
  # alla posizione 103 nel dataframe, sarà dunque il nostro 103-esimo nodo. Per mantenere 
  # la coerenza del formato input avremo 103-1 = 102.
  # L'arco verrà inserito tra il nodo 0 e il nodo 102 e sara' un arco etichettato da INHIBITION.
  # NB: chiaramente la definizione dell'indice varia da paziente a paziente, strettamente legata
  # al numero di nodi/geni che sono presenti/che appaiono per quel paziente
  for(i in 1:nrow(tmp1_edges)){
    tmp1_edges$Start[i] <- (which(tmp1_edges$Start[i] == tmp_patient_nodes$Id)-1)
    tmp1_edges$End[i] <- (which(tmp1_edges$End[i] == tmp_patient_nodes$Id)-1)
  }
  write.table(tmp1_edges, name, append = TRUE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  close(fileP)
}

#STEP 2 applicato. Ricevo in output m file con conteggi dei motivi presenti nelle reti di ciascun paziente


#STEP 3 caricare gli output in un unico dataset che verra' suddiviso in 75% training set e 25% test set
# testare quale dei classificatori ha il livello di accuratezza piu' alto. 



motifset <- data.frame(Motif_Nodes = factor(), 
                       Motif_edges = factor())

breast_cancer_df <- data.frame(matrix(ncol=16249, nrow=0))

for(i in 1:m){
  name <- paste("./results/", c_matrix[i], ".txt", sep = "")
  tmp_motifset <- read.table(name, header=TRUE, sep ="\t")
  tmp_motifset <- tmp_motifset[,-3]
  motifset <- bind_rows(motifset, tmp_motifset)
}
#aggiungere ulteriore colonna di identificativo
motifset <- motifset[!duplicated(motifset),]
motifset$Id <- 0
for(i in 1:nrow(motifset)){
  motifset[i,3] <- i
}
#fare poi l'intersezione con ogni file cosi' da tenere solo i motivi effettivamente presenti
#estrarre il conteggio e aggiungerlo in breast_cancer_df nella colonna corrispondente all'ID

labels <- data.frame(AJCC_PATHOLOGIC = character(),
                     STAGE = character(),
                     stringsAsFactors = FALSE)

for(i in 1:m){
  name <- paste("./results/", c_matrix[i], ".txt", sep = "")
  tmp_motifset <- read.table(name, header=TRUE, sep ="\t")
  tmp_breast_cancer_df <- tmp_motifset
  tmp_breast_cancer_df <- inner_join(tmp_breast_cancer_df, motifset)
  tmp_breast_cancer_df$Num_occ_input_graph <- as.integer(tmp_breast_cancer_df$Num_occ_input_graph)
  for(j in 1:nrow(tmp_breast_cancer_df)){
    breast_cancer_df[i, tmp_breast_cancer_df[j, 4]] <- tmp_breast_cancer_df[j, 3]
  }
  
  labels[i, 1] <- as.character(clinical[which(c_matrix[i] == clinical$SUBTYPE), 2])
  labels[i, 2] <- as.character(clinical[which(c_matrix[i] == clinical$SUBTYPE), 3])
}

#modificare il nome delle colonne e togliere la X
colnam <- colnames(breast_cancer_df)
b <- length(colnam)
for(q in 1:b)
  colnam[q] <- gsub("X", "Motivo", colnam[q])

colnames(breast_cancer_df) <- colnam


#rimozione di tutti i valori NA
breast_cancer_df <- as.matrix(breast_cancer_df)
breast_cancer_df[is.na(breast_cancer_df)] <- 0
breast_cancer_df <- data.frame(breast_cancer_df, stringAsFactor = FALSE)
breast_cancer_df <- breast_cancer_df[,1:(ncol(breast_cancer_df)-1)]

labels <- as.matrix(labels)
labels[is.na(labels)] <- 0 # azzera i valori NA su AJCC_PATHOLOGIC
labels[labels == ""] <- 0 # azzera i valori nulli su STAGE
# raggruppare gli stadi cosi' da lasciare solamente STAGE I, STAGE II, etc.
labels <- data.frame(labels)
stage_cleaned <- labels[2]
for(i in 1:nrow(stage_cleaned)){
  if(grepl("B", stage_cleaned[i, 1]))
    stage_cleaned[i, 1] <- gsub("B", "", stage_cleaned[i, 1])
  else if(grepl("C", stage_cleaned[i, 1]))
    stage_cleaned[i, 1] <- gsub("C", "", stage_cleaned[i, 1])
  else if(substr(stage_cleaned[i, 1], nchar(stage_cleaned[i, 1]), nchar(stage_cleaned[i, 1])) == "A")
    stage_cleaned[i, 1] <- substr(stage_cleaned[i, 1], 1, (nchar(stage_cleaned[i, 1])-1))
}
labels <- labels[1]
labels <- cbind(labels, stage_cleaned)

# conversione a factor
# i livelli di STAGE vengono ridotti da 13 a 6, compreso di livello nullo
labels$STAGE <- factor(labels$STAGE)
labels$AJCC_PATHOLOGIC <- factor(labels$AJCC_PATHOLOGIC)

breast_cancer_df <- breast_cancer_df[,1:16249]

breast_cancer_df <- cbind(breast_cancer_df, labels)

#breast_cancer_df_test <- breast_cancer_df_test[1:100, ]

perc.splitting <- 0.75
nobs.training <- round(perc.splitting * nrow(breast_cancer_df))
sampled.pos <- sample(1:nrow(breast_cancer_df), nobs.training)

first_breast_training <- breast_cancer_df[sampled.pos, -16251]
second_breast_training <- breast_cancer_df[sampled.pos, -16250]

first_breast_testing <- breast_cancer_df[-sampled.pos, -16251]
second_breast_testing <- breast_cancer_df[-sampled.pos, -16250]

first_true_classes <- first_breast_testing[, 16250]
second_true_classes <- second_breast_testing[, 16250]

first_breast_testing <- first_breast_testing[, -16250]
second_breast_testing <- second_breast_testing[, -16250]

#ALLENARE SU AJCC_PATHOLOGIC E STAGE SEPARATAMENTE E FORNIRE I RISULTATI

#VALIDAZIONE DEI CLASSIFICATORI pacchetto caret
library(caret)

#-----AJCC_PATHOLOGIC-----#
train.control <- trainControl(method="cv", number=4)
first.breast.caret.cart <- train(AJCC_PATHOLOGIC ~ ., data = first_breast_training, trControl=train.control, 
                                 method="rpart", na.action = na.omit)
first.breast.caret.rf <- train(AJCC_PATHOLOGIC ~ ., data = first_breast_training, trControl=train.control, 
                               method="rf", na.action = na.omit)
first.breast.caret.nb <- train(AJCC_PATHOLOGIC ~ ., data = first_breast_training, trControl=train.control, 
                               method="nb", na.action = na.omit)
first.breast.caret.svm <- train(AJCC_PATHOLOGIC ~ ., data = first_breast_training, trControl=train.control, 
                                method="svmRadial", na.action = na.omit)
#--------------------------------#
predict.caret.cart <- predict(first.breast.caret.cart, first_breast_testing)
ajcc.caret.cart.results <- data.frame(real = first_true_classes, predicted = predict.caret.cart)
first.caret.cart.accuracy <- sum(ajcc.caret.cart.results$real == ajcc.caret.cart.results$predicted)/nrow(ajcc.caret.cart.results)

predict.caret.rf <- predict(first.breast.caret.rf, first_breast_testing)
ajcc.caret.rf.results <- data.frame(real = first_true_classes, predicted = predict.caret.rf)
first.caret.rf.accuracy <- sum(ajcc.caret.rf.results$real == ajcc.caret.rf.results$predicted)/nrow(ajcc.caret.rf.results)

predict.caret.nb <- predict(first.breast.caret.nb, first_breast_testing)
ajcc.caret.nb.results <- data.frame(real = first_true_classes, predicted = predict.caret.nb)
first.caret.nb.accuracy <- sum(ajcc.caret.nb.results$real == ajcc.caret.nb.results$predicted)/nrow(ajcc.caret.nb.results)

predict.caret.svm <- predict(first.breast.caret.svm, first_breast_testing)
ajcc.caret.svm.results <- data.frame(real = first_true_classes, predicted = predict.caret.svm)
first.caret.svm.accuracy <- sum(ajcc.caret.svm.results$real==ajcc.caret.svm.results$predicted)/nrow(ajcc.caret.svm.results)

capture.output(confusionMatrix(ajcc.caret.cart.results$predicted, ajcc.caret.cart.results$real), file = "./cart-results-ajcc.txt")
capture.output(confusionMatrix(ajcc.caret.rf.results$predicted, ajcc.caret.rf.results$real), file = "./rf-results-ajcc.txt")
capture.output(confusionMatrix(ajcc.caret.nb.results$predicted, ajcc.caret.nb.results$real), file = "./nb-results-ajcc.txt")
capture.output(confusionMatrix(ajcc.caret.svm.results$predicted, ajcc.caret.svm.results$real), file = "./svm-results-ajcc.txt")

#-----STAGE----#
train.control <- trainControl(method="cv", number = 4)
#train.control.nb <- trainControl(method="cv", sampling = "up", number = 4)
second.breast.caret.cart <- train(STAGE ~ ., data = second_breast_training, trControl=train.control, 
                                  method="rpart", na.action = na.omit)
second.breast.caret.rf <- train(STAGE ~ ., data = second_breast_training, trControl=train.control, 
                                method="rf", na.action = na.omit)
second.breast.caret.nb <- train(STAGE ~ ., data = second_breast_training, trControl=train.control, 
                                method="nb", na.action = na.omit)
second.breast.caret.svm <- train(STAGE ~ ., data = second_breast_training, trControl=train.control, 
                                 method="svmRadial", na.action = na.omit)
#--------------------------------#
predict.caret.cart <- predict(second.breast.caret.cart, second_breast_testing)
stage.caret.cart.results <- data.frame(real = second_true_classes, predicted = predict.caret.cart)
first.caret.cart.accuracy <- sum(stage.caret.cart.results$real == stage.caret.cart.results$predicted)/nrow(stage.caret.cart.results)

predict.caret.rf <- predict(second.breast.caret.rf, second_breast_testing)
stage.caret.rf.results <- data.frame(real = second_true_classes, predicted = predict.caret.rf)
second.caret.rf.accuracy <- sum(stage.caret.rf.results$real == stage.caret.rf.results$predicted)/nrow(stage.caret.rf.results)

predict.caret.nb <- predict(second.breast.caret.nb, second_breast_testing)
stage.caret.nb.results <- data.frame(real = second_true_classes, predicted = predict.caret.nb)
second.caret.nb.accuracy <- sum(stage.caret.nb.results$real == stage.caret.nb.results$predicted)/nrow(stage.caret.nb.results)

predict.caret.svm <- predict(second.breast.caret.svm, second_breast_testing)
stage.caret.svm.results <- data.frame(real = second_true_classes, predicted = predict.caret.svm)
second.caret.svm.accuracy <- sum(stage.caret.svm.results$real == stage.caret.svm.results$predicted)/nrow(stage.caret.svm.results)

capture.output(confusionMatrix(stage.caret.cart.results$predicted, stage.caret.cart.results$real), file = "./cart-results-stage.txt")
capture.output(confusionMatrix(stage.caret.rf.results$predicted, stage.caret.rf.results$real), file = "./rf-results-stage.txt")
capture.output(confusionMatrix(stage.caret.nb.results$predicted, stage.caret.nb.results$real), file = "./nb-results-stage.txt")
capture.output(confusionMatrix(stage.caret.svm.results$predicted, stage.caret.svm.results$real), file = "./svm-results-stage.txt")