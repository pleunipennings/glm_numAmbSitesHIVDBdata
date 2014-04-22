Data<-read.csv("Pleuni.data2.txt",sep=" ")

Data$PercAmb = Data$numAmbig/Data$RTlen
Data$numNotAmbig=Data$RTlen - Data$numAmbig
DRMcolums=which(substr(names(Data),1,3)=="pos")
Data$numDRMs <- apply(Data[,DRMcolums], 1, function(x) sum(x))
#make binary whether they have a DRM or not
Data$numDRMbinary <- 0
Data$numDRMbinary[which(Data$numDRMs>0)]<-1

#PtID is a unique identifier for this dataset.

#I think I need to code the different regimens.
ListDrugnames<-sort(unique(unlist( strsplit(as.character(unique(Data$Regimen)), "and"))))

for (Drugname in ListDrugnames){

newcolumn=length(Data[1,])+1
Data[,newcolumn]<-0
names(Data)[newcolumn]<-Drugname
if (Drugname=="3TC") names(Data)[newcolumn]<-"d3TC"
Data[,newcolumn][grep(Drugname,Data$Regimen)]<-1

}

###############################
#organize in classes of drugs
#also count number of drugs
if (TRUE){
    
newcolumn=length(Data[1,])+1
Data[,newcolumn]<-0
names(Data)[newcolumn]<-"ClassNRTI"
#length(which(Data$ClassNRTI==1))

for (d in c("d3TC",
"ABC",
"AZT",
"D4T",
"DDC",
"DDI",
"FTC",
"NRTI",
"TDF")){
drugcolumn=which(names(Data)==d)
Data$ClassNRTI[which(Data[,drugcolumn]==1)]<-1}



newcolumn=length(Data[1,])+1
Data[,newcolumn]<-0
names(Data)[newcolumn]<-"ClassNNRTI"

for (d in c("DLV",
"NNRTI",
"EFV",
"NVP",
"RPV")){
    drugcolumn=which(names(Data)==d)
    Data$ClassNNRTI[which(Data[,drugcolumn]==1)]<-1}


newcolumn=length(Data[1,])+1
Data[,newcolumn]<-0
names(Data)[newcolumn]<-"ClassPI"

for (d in c("APV",
"ATV",
"DRV",
"FPV",
"IDV",
"LPV",
"NFV",
"RTV",
"SQV",
"PI")){
    drugcolumn=which(names(Data)==d)
    Data$ClassPI[which(Data[,drugcolumn]==1)]<-1}


newcolumn=length(Data[1,])+1
Data[,newcolumn]<-0
names(Data)[newcolumn]<-"ClassINT"

for (d in c("RAL")){
    drugcolumn=which(names(Data)==d)
    Data$ClassINT[which(Data[,drugcolumn]==1)]<-1}


newcolumn=length(Data[1,])+1
Data[,newcolumn]<-0
names(Data)[newcolumn]<-"ClassUnkown"

for (d in c(
"RTI",
"aAPA",
"Unknown")){
    drugcolumn=which(names(Data)==d)
    Data$ClassUnkown[which(Data[,drugcolumn]==1)]<-1}


newcolumn=length(Data[1,])+1
Data[,newcolumn]<-0
names(Data)[newcolumn]<-"NumDrugs"
firstdrugcolumn = which(names(Data)=="d3TC")
lastdrugcolumn = which(names(Data)=="Unknown")
DrugColums<-firstdrugcolumn:lastdrugcolumn
Data$NumDrugs=apply(Data[,DrugColums], 1, function(x) sum(x))

hist(Data$NumDrugs,breaks=c(-0.5,0.5,1.5,2.5,3.5,4.5,5.5,10))
}
#############################


#FOCUS only on patients without PI's
DataRTI<-Data[Data$ClassPI==0,]
DataRTI<-DataRTI[DataRTI$ClassUnkown==0,]
DataRTI<-DataRTI[!is.na(DataRTI$numDRMs),]
#DataRTI<-DataRTI[DataRTI$IsolateYear>1994,]

#DataRTI$IsolateYear<- ordered(DataRTI$IsolateYear, levels = sort(unique(DataRTI$IsolateYear)))

y<-cbind(DataRTI$numAmbig,DataRTI$numNotAmbig)

attach(DataRTI)

if (FALSE){
model002<-glm(y~numDRMbinary*IsolateYear*ClassNNRTI*NumDrugs,quasibinomial)
summary(model002)

model003<-update(model002,~.-numDRMbinary:IsolateYear:ClassNNRTI:NumDrugs)
anova(model003,model002,test="F")
#no difference

model004<-update(model003,~.-numDRMbinary:IsolateYear:ClassNNRTI)
anova(model003,model004,test="F")
#no difference

model005<-update(model004,~.-IsolateYear:ClassNNRTI)
anova(model005,model004,test="F")
#no difference

pdf("Try.pdf")

par(mfrow=c(2,1))
years=seq(min(DataRTI$IsolateYear),max(DataRTI$IsolateYear)-1,by=2)
reps=length(years)

for (numdrugs in c(1,2,3)){

for (numMuts in 0:1){
    plot(predict(model005,list(numDRMbinary=rep(numMuts,reps),IsolateYear=years,ClassNNRTI=rep(0,reps),NumDrugs=rep(numdrugs,reps)),type="response"),main=paste("nummuts=",numMuts),ylim=c(0,0.1),xaxt="n",ylab="AmbReads",xlab="years")
    points(predict(model005,list(numDRMbinary=rep(numMuts,reps),IsolateYear=years,ClassNNRTI=rep(1,reps),NumDrugs=rep(numdrugs,reps)),type="response"),col=2)
    
    for (TreatedwNNRTI in c(0,1)){
        for (y in years){
        rawdata=DataRTI$PercAmb[which(DataRTI$numDRMbinary==numMuts&DataRTI$IsolateYear==y&DataRTI$ClassNNRTI==TreatedwNNRTI&DataRTI$NumDrugs==numdrugs)]
        if (length(rawdata)>0){
            print(paste("year=",y,"numMuts",numMuts,"ClassNNRTI",0,"NumDrugs",numdrugs,"MeanNumAmbReads",mean(rawdata)))
            points(which(years==y), mean(rawdata),col=TreatedwNNRTI+1,pch=4)}}
    }
}}

    dev.off()}


model000<-glm(y~numDRMbinary*IsolateYear,quasibinomial)
summary(model000)

pdf("SimplestModel.pdf")

par(mfrow=c(2,1))

years=seq(min(DataRTI$IsolateYear),max(DataRTI$IsolateYear)-1,by=2)
years=min(DataRTI$IsolateYear):max(DataRTI$IsolateYear)
reps=length(years)

    for (numMuts in 0:1){
        plot(predict(model000,list(numDRMbinary=rep(numMuts,reps),IsolateYear=years),type="response"),main=paste("nummuts =",numMuts),ylim=c(0,0.05),xaxt="n",ylab="AmbReads",xlab="years",pch=16)
        axis(1,at=1:reps,labels=years)
        
        for (y in years){
            rawdata=DataRTI$PercAmb[which(DataRTI$numDRMbinary==numMuts&DataRTI$IsolateYear==y)]
            if (length(rawdata)>0){
                print(paste("year=",y,"numMuts",numMuts,"MeanNumAmbReads",mean(rawdata)))
                points(which(years==y), mean(rawdata),col=2,pch=4)
                lines(c(which(years==y),which(years==y)),c(summary(rawdata)[2],summary(rawdata)[5]))
            
            }}
    }

dev.off()


