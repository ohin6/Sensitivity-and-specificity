#import packages
require(ggplot2)
require(pROC)

##read in data

RawData = read.csv(file.choose(), header = FALSE)

# read in dataframe without the first four lines-just the data
data_raw<-read.csv(file.choose(), skip =3)

#Set column names in data_raw from copying row 2 in RawData
colnames(data_raw)=RawData[2,]

#replace 0 with NA
#data_raw[data_raw == 0] <- NA

#remove empty cells
data_raw<-data_raw[1:800,1:5]

#Create dataframe

df = data_raw


#calculating canal paresis 
# CP = ((WR+CR)-(WL+CL))/(WR+WL+CR+CL)*100
df$canal_paresis = ((df$WR+df$CR)-(df$WL+df$CL))/(df$WR+df$WL+df$CR+df$CL)*100

#Calculating directional preponderance
#DP = ((WR+CL) - (WL+CR)) / (WR+WL+CR+CL) *100
df$dir_pro = ((df$WR+df$CL)-(df$WL+df$CR))/(df$WR+df$WL+df$CR+df$CL)*100


# Create monothermal screen column
df$monothermal = NA


#Create function specific depending on order ie warm or Cold
Func_W = which(df$Order=="W")
Func_C = which(df$Order=="C")


# calculate monothermal screen for order Warm 
df[Func_W,"monothermal"] = (df[Func_W, "WR"]-df[Func_W,"WL"])/(df[Func_W, "WR"]+df[Func_W,"WL"])*100

# calculate monothermal screen for order Cold 
df[Func_C,"monothermal"] = (df[Func_C, "CR"]-df[Func_C,"CL"])/(df[Func_C, "CR"]+df[Func_C,"CL"])*100

#####Create column to determine whether significant (True) not significant (False)####

#Significance for CP (using absolute values to compensate for negative values)
df$CP_signif=ifelse(abs(df$canal_paresis)>20, TRUE, FALSE)

#Significance for DP (using absolute values to compensate for negative values)
df$DP_signif=ifelse(abs(df$canal_paresis)>20, TRUE, FALSE)


#Determine whether "gold" standard method was significant or non significant 

df$gold = ifelse(df$DP_signif==FALSE & df$CP_signif==FALSE,FALSE,TRUE)

#create column for TN, FN, TP, TN
df$rating = NA

#Create column for monosignificance
df$mono_sig = NA


######Creates a function that calculates sensitivity and specificity with given thresholds of monothermal screen threshold (signif)#####
Sen_Spec=function(df,signif){

  #Boolean to return whether absolute value of monothermal screen
  df$mono_sig = abs(df$monothermal)>signif
  
  
  #Input TN, FN, TP, TN into "rating" column
  #Input TP (true positive) in ratings column where monothermal = True and gold = True
  df[which(df$gold&df$mono_sig),"rating"]<-"TP"
  #Input TN (true negative) in ratings column where monothermal = True and gold = True
  df[which(!df$gold&!df$mono_sig),"rating"]<-"TN"
  #Input TP (true positive) in ratings column where monothermal = True and gold = True
  df[which(!df$gold & df$mono_sig),"rating"]<-"FP"
  #Input TN (true negative) in ratings column where monothermal = True and gold = True
  df[which(df$gold & !df$mono_sig ),"rating"]<-"FN"
  
  #calculate sensitivity and specificity where
  ##include na.rm = T otherwise does not work - removes Na's
  # sensitivity = sum of TP/(sum of TP + sum of FN)
  sensitivity = sum(df$rating=="TP",na.rm = T)/(sum(df$rating=="TP",na.rm = T)+sum(df$rating=="FN",na.rm = T))
  # specificity = sum of TN/(sum of TN + sum of FP)
  Specificity = sum(df$rating=="TN",na.rm = T)/(sum(df$rating=="TN",na.rm = T)+sum(df$rating=="FP",na.rm = T))
  
  
  #Return specificity and sensitivity
  return(c(sensitivity, Specificity))
  
}



#set min and max significance 
min_sig=0
max_sig=100

#sequence to iterate over
sig_range=seq(from=min_sig,to=max_sig,length.out=100)
#round to the nearest whole number to make the labels nicer!
sig_range=round(sig_range, 0)

#make new dataset to match the number of values we're using
sens_spec_dataset=data.frame(matrix(nrow=(length(sig_range)),ncol=3))
colnames(sens_spec_dataset)=c("significance","sensitivity","specificity")

#first column in dataset is the significance value
sens_spec_dataset[,1]<-sig_range

#loop through length of sig_range
for (i in 1:length(sig_range)){
  #call function and set output to res
  res=Sen_Spec(df,sig_range[i])
  #store sensitivity in dataset for current significance value
  sens_spec_dataset[i,2]=res[1]
  #store specificity in dataset
  sens_spec_dataset[i,3]=res[2]
}

#false positive rate is 1 - specificity
sens_spec_dataset$false_pos_rate<- 1 - sens_spec_dataset$specificity

#make ROC plot (ggplot)
g<-ggplot(data=sens_spec_dataset,aes(false_pos_rate,sensitivity))+
  geom_text(label=sens_spec_dataset$significance,nudge_x = 0.05, nudge_y = 0.05, check_overlap = T)+
  geom_point()+theme_bw()+ylab("True Positive Rate (specificity)")+xlab("False Positive Rate (1-sensitivity)")
g

#pROC package prefers true and false to be 1s and 0s, multiplying by one does this
df$gold<-df$gold*1
#use roc function from the overall significance and the results of the screen
roc_obj<-roc(df$gold,abs(df$monothermal),na.rm=TRUE)
#get AUC- EXCELLENT
auc(roc_obj)
#plot ROC
plot(roc_obj)
  
Sen_Spec(df,20)

#Delete columns
#df = subset(df, select = -c(CP_significance))
