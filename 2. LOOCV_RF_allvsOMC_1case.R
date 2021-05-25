ntree = 1000
mth = "RF"
output.type = "BestResCat"

library(randomForest)

library(doParallel)
if(!exists("cl")) {
  cl <- makeCluster(20)
  registerDoParallel(cl)
}

library(caret)

if(all(dir("../" ) != paste(mth,".LOOCV2_predicted", sep = ""))) dir.create(paste("../",mth,".LOOCV2_predicted", sep = ""))
if(all(dir("../" ) != "results_selfeats")) dir.create("../results_selfeats/")

## LOADING DATA
treatments = read.csv("../inp/pdxe_treatments50.csv", sep = ",")
ntreatments = nrow(treatments)
cancer.type = read.csv("../inp/pdxe_cancertype.csv", sep= ",")

## INPUT BY USERS
cancertype = "BRCA"
feat.type = "SNV"
treatmentName = "paclitaxel"


## MODEL BUILDING
sel.meth = ifelse(feat.type %in% c("SNV", "CNA"),"fisher","ttest")

treatmentid = treatments$Treatment_ID[which(treatments$Treatment_Details == treatmentName)]
cat(paste("Treatment ",treatmentid," - ",treatmentName, "\n",sep = "")) 

inpfile = paste("../processedData.",feat.type,"/", cancertype, "_Treatment",treatmentid,"_",output.type,"_",feat.type,".csv",sep = "")    
data = read.csv(file = inpfile,sep = ",", header = TRUE)
dataBin = ifelse(data$BestResCategory %in% c("PR","CR","SD"),"Sen","Res")
data.output = data[,1]
count = length(unique(data$Model))
maxfeats = round(count/2)  ## Set Maxfeats = N/2

if(all(dir(paste("../",mth,".LOOCV2_predicted", sep = "")) != cancertype)) dir.create(paste("../",mth,".LOOCV2_predicted/",cancertype, sep = ""))

if(all(dir(paste("../",mth,".LOOCV2_predicted/",cancertype, sep = "")) != paste(treatmentName,feat.type,sep = "."))) {
dir.create(paste("../",mth,".LOOCV2_predicted/",cancertype,"/",paste(treatmentName,feat.type,sep = "."), sep = ""))}
if(all(dir(paste("../",mth,".LOOCV2_predicted/",cancertype, sep = "")) != paste(treatmentName,feat.type,"allfts",sep = "."))) {
dir.create(paste("../",mth,".LOOCV2_predicted/",cancertype,"/",paste(treatmentName,feat.type,"allfts",sep = "."), sep = ""))}

SPEC.tst = c()
RECALL.tst = c()
PREC.tst = c()
F1.tst = c()
nSens.tst = c()
best.feat.med = c()
MCC.tst = c()
MCC.tst.all = c()
# seed1 = 5678
seeds = c(5678,1111,2222,3333,4444,5555,6666,7777,8888,9999)
# n = 1
M2 = foreach (n = 1:length(seeds),.packages=c('randomForest','caret','doParallel'),.combine = 'rbind', .inorder = TRUE )%do%{
  best.feats.list = c()
  allpred.all = c()
  allpred = c()
  seed1 = seeds[n]
  folds.outer <- caret::createFolds(dataBin,k = nrow(data))
  print.pred = data.frame(matrix(NA, ncol = 4, nrow = 0))
  colnames(print.pred) = c("Fold", "PDXindex","Res","Sen")
  print.pred.all = data.frame(matrix(NA, ncol = 4, nrow = 0))
  colnames(print.pred.all) = c("Fold", "PDXindex","Res","Sen")
  
  # j = 1
  for(j in 1:nrow(data)){
    train.data = data[-folds.outer[[j]],] #Set the training set
    test.data = data[folds.outer[[j]],] #Set the validation set
    
    train.output = train.data[,1]  
    test.output = test.data[,1]  
    
    trainBin = as.factor(ifelse(train.output %in% c("PR","CR","SD"),"Sen","Res"))
    train.data1 = train.data[,-c(1,ncol(train.data))]
    test.data1 = test.data[,-c(1,ncol(test.data))]
    
    mtry.all = round(sqrt(ncol(train.data1)))
    sens = length(which(train.output %in% c("CR","PR","SD")))
    res = length(which(train.output == "PD"))
    prop.sens = sens/nrow(train.data1)
    prop.res  = res/nrow(train.data1)
    
    set.seed(seed1)
    rf.all = randomForest(train.data1,trainBin,mtry = mtry.all,importance = TRUE,ntree = ntree, classwt = c(Res = prop.sens,Sen = prop.res)) ## class weighting
    # rf2 = randomForest(train.data,trainBin,mtry = mtry1,importance = TRUE,ntree = ntree, sampsize = c(min(c(sens,res)),min(c(sens,res)))) ## class weighting
    
    # On test data
    pred.tst.prob.all = predict(rf.all,test.data1, type = "prob")
    pred.tst.all = c()
    thres = 0.5
    for (r in 1:nrow(pred.tst.prob.all)){
      if(pred.tst.prob.all[r,1] > thres) {
        pred.tst.all[r] = "Res"
      }else{
        if(pred.tst.prob.all[r,1] == thres) {
          set.seed(5678)
          pred.tst.all[r] = sample(c("Res","Sen"),1)
        } else {
          pred.tst.all[r] = "Sen"}
      }
    }
    names(pred.tst.all) = rownames(test.data)
    rownames(pred.tst.prob.all)= rownames(test.data)
    allpred.all = c (allpred.all,pred.tst.all)
    
    pred.tst.prob.all = cbind(rep(j, nrow(pred.tst.prob.all)),rownames(pred.tst.prob.all),pred.tst.prob.all)
    colnames(pred.tst.prob.all) = c("Fold", "PDXindex","Res","Sen")
    print.pred.all = rbind(print.pred.all,pred.tst.prob.all)
    
    folds.inner <- caret::createFolds(trainBin,k = nrow(train.data))
    
    MCC.tst.inner = c()
    M3 = foreach(l = 1:nrow(train.data),.packages=c('randomForest','caret'),.combine = 'rbind', .inorder = TRUE) %dopar%{
      train.data.in = train.data[-folds.inner[[l]],] #Set the training set
      test.data.in = train.data[folds.inner[[l]],] #Set the validation set
      train.output.in = train.data.in[,1]  
      test.output.in = test.data.in[,1]  
      
      trainBin.in = as.factor(ifelse(train.output.in %in% c("PR","CR","SD"),"Sen","Res"))
      
      train.data.in = train.data.in[,-c(1,ncol(train.data.in))]
      test.data.in = test.data.in[,-c(1,ncol(test.data.in))]
      
      oob.list = c()
      nfeats.list = c()
      mtry1 = round(sqrt(ncol(train.data.in)))
      sens = length(which(train.output.in %in% c("CR","PR","SD")))
      res = length(which(train.output.in == "PD"))
      prop.sens.in = sens/nrow(train.data.in)
      prop.res.in  = res/nrow(train.data.in)
      set.seed(seed1)
      rf1 = randomForest(train.data.in,trainBin.in,mtry = mtry1,importance = TRUE,ntree = ntree, classwt = c(Res = prop.sens.in,Sen = prop.res.in)) ## class weighting
      # rf2 = randomForest(train.data,trainBin,mtry = mtry1,importance = TRUE,ntree = ntree, sampsize = c(min(c(sens,res)),min(c(sens,res)))) ## class weighting
      
      pred.tst.prob.all.in = predict(rf1,test.data.in, type = "prob")
      pred.tst.all.in = c()
      thres = 0.5
      for (r in 1:nrow(pred.tst.prob.all.in)){
        if(pred.tst.prob.all.in[r,1] > thres) {
          pred.tst.all.in[r] = "Res"
        }else{
          if(pred.tst.prob.all.in[r,1] == thres) {
            set.seed(5678)
            pred.tst.all.in[r] = sample(c("Res","Sen"),1)
          } else {
            pred.tst.all.in[r] = "Sen"}
        }
      }
      
      u=apply(train.data.in, 2, function(col) length(unique(col)))
      selCols=which(u != 1)
      
      train.data.in = train.data.in[,selCols]
      
      if(sel.meth == "fisher"){
        pvals.ft = apply(train.data.in,2,function(x){
          contigency.tab = matrix(c(length(x[which(trainBin.in == "Sen" & x == 1)]),
                                    length(x[which(trainBin.in == "Sen" & x == 0)]),
                                    length(x[which(trainBin.in == "Res" & x == 1)]),
                                    length(x[which(trainBin.in == "Res" & x == 0)])),
                                  nrow = 2,
                                  dimnames = list(Gene = c("Mutated", "WT"),
                                                  Response = c("Sen", "Res")))
          fisher.test(contigency.tab, alternative = "two.sided")$p.value
        })
      }else{
        pvals.tt = apply(train.data.in,2,
                         function(x){
                           # check if the arrays have a single value in all cases and 
                           # if these values are different between sensitive and resistant cases
                           if (length(unique(x[which(trainBin.in == "Sen")])) == 1 & 
                               length(unique(x[which(trainBin.in == "Res")])) == 1){
                             if (unique(x[which(trainBin.in == "Sen")]) != unique(x[which(trainBin.in == "Res")])){
                               -1 } else { 1 }
                           }else {
                             a = t.test(x[which(trainBin.in == "Sen")],x[which(trainBin.in == "Res")],paired = F)
                             a$p.value}
                         })
      }
      
      pred.tst.in = c()
      for (nfeats in 2:maxfeats){
        if(sel.meth == "fisher") topfeats = names(pvals.ft[order(pvals.ft, decreasing = F)][1:nfeats]) else topfeats = names(pvals.tt[order(pvals.tt, decreasing = F)][1:nfeats])
        
        train.data.in2 = train.data.in[,which(colnames(train.data.in) %in% topfeats)]
        test.data.in2  = test.data.in[,which(colnames(test.data.in) %in% topfeats)]
        
        mtry2 = round(sqrt(ncol(train.data.in2)))
        set.seed(seed1)
        rf2 = randomForest(train.data.in2,trainBin.in,mtry = mtry2,importance = TRUE,ntree = ntree, classwt = c(Res = prop.sens.in,Sen = prop.res.in)) ## class weighting
        # rf2 = randomForest(train.data,trainBin,mtry = mtry1,importance = TRUE,ntree = ntree, sampsize = c(min(c(sens,res)),min(c(sens,res)))) ## class weighting
        
        pred.tst.prob.sel = predict(rf2,test.data.in2,type = "prob")
        
        pred.tst.sel = c()
        thres = 0.5
        for (r in 1:nrow(pred.tst.prob.sel)){
          if(pred.tst.prob.sel[r,1] > thres) {
            pred.tst.sel[r] = "Res"
          }else{
            if(pred.tst.prob.sel[r,1] == thres) {
              set.seed(5678)
              pred.tst.sel[r] = sample(c("Res","Sen"),1)
            } else {
              pred.tst.sel[r] = "Sen"}
          }
        }
        pred.tst.in = c(pred.tst.in,pred.tst.sel)
      }
      c(pred.tst.all.in,pred.tst.in,rownames(test.data.in))
    }
    calculate.mcc = function(predict, observe){
      TP = length(which(predict == "Sen" & observe %in% c("PR","CR","SD")))    # cells predicted sensitive which are actually sensitive (cell sensitive is +ve)
      TN = length(which(predict == "Res" & observe == "PD"))  # cells predicted resistant which are actually resistant (cell resistant is -ve)
      FP = length(which(predict == "Sen" & observe == "PD")) # cells predicted sensitive which are actually resistant (wrong prediction)
      FN = length(which(predict == "Res" & observe %in% c("PR","CR","SD")))       # cells predicted resistant which are actually sensitive (wrong prediction)
      dum1 = ((TP*TN)-(FP*FN))
      dum2 = sqrt(1.0*(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) # 1.0 makes the multiplication float, thus avoiding integer overflow
      if(dum2 == 0) dum2 = 1
      mcc = dum1/dum2
    } 
    train.output = train.output[match(rownames(train.data),as.numeric(M3[,ncol(M3)]))]
    med.MCC = apply(M3[,-ncol(M3)], 2, function(x) calculate.mcc(x,train.output))
    
    nfeats = if (med.MCC[1] <= max(med.MCC[-1])) which.max(med.MCC[-1])+1 else ncol(data)-2 
    best.feats.list = c(best.feats.list,nfeats)
    train.data = train.data[,-c(1,ncol(train.data))]
    test.data = test.data[,-c(1,ncol(test.data))]
    
    u=apply(train.data, 2, function(col) length(unique(col)))
    selCols = which(u != 1)
    
    train.data = train.data[,selCols]
    
    if(sel.meth == "fisher"){
      pvals.ft = apply(train.data,2,function(x){
        contigency.tab = matrix(c(length(x[which(trainBin == "Sen" & x == 1)]),
                                  length(x[which(trainBin == "Sen" & x == 0)]),
                                  length(x[which(trainBin == "Res" & x == 1)]),
                                  length(x[which(trainBin == "Res" & x == 0)])),
                                nrow = 2,
                                dimnames = list(Gene = c("Mutated", "WT"),
                                                Response = c("Sen", "Res")))
        fisher.test(contigency.tab, alternative = "two.sided")$p.value
      })
    }else{
      pvals.tt = apply(train.data,2,
                       function(x){
                         # check if the arrays have a single value in all cases and 
                         # if these values are different between sensitive and resistant cases
                         if (length(unique(x[which(trainBin == "Sen")])) == 1 & 
                             length(unique(x[which(trainBin == "Res")])) == 1){
                           if (unique(x[which(trainBin == "Sen")]) != unique(x[which(trainBin == "Res")])){
                             -1 } else { 1 }
                         }else {
                           a = t.test(x[which(trainBin == "Sen")],x[which(trainBin == "Res")],paired = F)
                           a$p.value}
                       })
    }
    if (nfeats == ncol(data) -2){
      pred.tst.prob = predict(rf.all,test.data, type = "prob")
      pred.tst = pred.tst.all
    } else {
      if(sel.meth == "fisher") topfeats = names(pvals.ft[order(pvals.ft, decreasing = F)][1:nfeats]) else topfeats = names(pvals.tt[order(pvals.tt, decreasing = F)][1:nfeats])
      
      train.data2 = train.data[,which(colnames(train.data) %in% topfeats)]
      test.data2  = test.data[,which(colnames(test.data) %in% topfeats)]
      
      mtry3 = round(sqrt(ncol(train.data2)))
      sens = length(which(train.output %in% c("CR","PR","SD")))
      res = length(which(train.output == "PD"))
      prop.sens = sens/nrow(train.data2)
      prop.res  = res/nrow(train.data2)
      set.seed(seed1)
      rf3 = randomForest(train.data2,trainBin,mtry = mtry3,importance = TRUE,ntree = ntree, classwt = c(Res = prop.sens,Sen = prop.res)) ## class weighting
      
      set.seed(5678)
      pred.tst.prob = predict(rf3,test.data2, type = "prob")
      pred.tst = c()
      thres = 0.5
      for (r in 1:nrow(pred.tst.prob)){
        if(pred.tst.prob[r,1] > thres) {
          pred.tst[r] = "Res"
        }else{
          if(pred.tst.prob[r,1] == thres) {
            set.seed(5678)
            pred.tst[r] = sample(c("Res","Sen"),1)
          } else {
            pred.tst[r] = "Sen"}
        }
      }
      names(pred.tst) = rownames(test.data) 
      rownames(pred.tst.prob)= rownames(test.data)
    }
    pred.tst.prob = cbind(rep(j, nrow(pred.tst.prob)),rownames(pred.tst.prob),pred.tst.prob)
    colnames(pred.tst.prob) = c("Fold", "PDXindex","Res","Sen")
    print.pred = rbind(print.pred,pred.tst.prob)
    allpred = c(allpred,pred.tst)
    
  }
  write.csv(print.pred, row.names = F,
            file = paste("../",mth,".LOOCV2_predicted/",cancertype,"/",paste(treatmentName,feat.type,sep = "."),"/Predicted_Replicate",n,".csv", sep = ""))  
  write.csv(print.pred.all, row.names = F,
            file = paste("../",mth,".LOOCV2_predicted/",cancertype,"/",paste(treatmentName,feat.type,"allfts",sep = "."),"/Predicted_Replicate",n,".csv", sep = ""))  
  
  best.feat.med = c(best.feat.med,median(best.feats.list))
  
  allpred = allpred[match(rownames(data),names(allpred))]
  allpred.all = allpred.all[match(rownames(data),names(allpred.all))]
  
  TP = length(which(allpred == "Sen" & data.output %in% c("PR","CR","SD")))    # cells predicted sensitive which are actually sensitive (cell sensitive is +ve)
  TN = length(which(allpred == "Res" & data.output == "PD"))  # cells predicted resistant which are actually resistant (cell resistant is -ve)
  FP = length(which(allpred == "Sen" & data.output == "PD")) # cells predicted sensitive which are actually resistant (wrong prediction)
  FN = length(which(allpred == "Res" & data.output %in% c("PR","CR","SD")))       # cells predicted resistant which are actually sensitive (wrong prediction)
  dum1 = ((TP*TN)-(FP*FN))
  dum2 = sqrt(1.0*(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) # 1.0 makes the multiplication float, thus avoiding integer overflow
  if(dum2 == 0) mcc = NaN else {
    mcc = dum1/dum2  
  }
  MCC.tst = c(MCC.tst,mcc)
  precision = TP/(TP+FP)
  recall = TP/(TP+FN)
  specif = TN/(TN+FP)
  f1 = 2*precision*recall/(precision+recall)
  nsens = sum(pred.tst == "Sen")
  
  TPa = length(which(allpred.all == "Sen" & data.output %in% c("PR","CR","SD")))    # cells predicted sensitive which are actually sensitive (cell sensitive is +ve)
  TNa = length(which(allpred.all == "Res" & data.output == "PD"))  # cells predicted resistant which are actually resistant (cell resistant is -ve)
  FPa = length(which(allpred.all == "Sen" & data.output == "PD")) # cells predicted sensitive which are actually resistant (wrong prediction)
  FNa = length(which(allpred.all == "Res" & data.output %in% c("PR","CR","SD")))       # cells predicted resistant which are actually sensitive (wrong prediction)
  dum1a = ((TPa*TNa)-(FPa*FNa))
  dum2a = sqrt(1.0*(TPa+FPa)*(TPa+FNa)*(TNa+FPa)*(TNa+FNa)) # 1.0 makes the multiplication float, thus avoiding integer overflow
  if(dum2a == 0) mcc.all = NaN else {
    mcc.all = dum1a/dum2a  
  }
  MCC.tst.all = c(MCC.tst.all,mcc.all)  
  precision.all = TPa/(TPa+FPa)
  recall.all = TPa/(TPa+FNa)
  f1.all = 2*precision.all*recall.all/(precision.all+recall.all)
  c(median(best.feats.list),mcc,precision,recall,f1,
    mcc.all, precision.all,recall.all,f1.all)
}
M2 = data.frame(rbind(M2, c(" ",apply(M2[,-1],2,median))))
M2 = data.frame(cbind(c(paste("Rep", seq(10)),"Median"), M2))
colnames(M2) = c("Replicates","Median.nfeats","MCC.OMC", "PREC.OMC", "RECALL.OMC","F1.OMC",
                 "MCC.all", "PREC.all", "RECALL.all","F1.all")
write.csv(M2, file = paste("../results_selfeats/single.", mth, "_",feat.type,"_",cancertype,"_", treatmentName,".nestedLOOCV.10x.csv", sep=""), row.names = F)
