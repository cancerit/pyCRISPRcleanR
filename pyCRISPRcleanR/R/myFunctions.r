library(pROC)

ROC_Curve<-function(FCsprofile,positives,negatives,display=TRUE,FDRth=NULL){

    FCsprofile<-FCsprofile[intersect(c(positives,negatives),names(FCsprofile))]

    predictions<-FCsprofile
    observations<-is.element(names(FCsprofile),positives)+0
    names(observations)<-names(predictions)

    RES<-roc(observations,predictions,direction = '>')

    if(display){
        plot(RES,col='blue',lwd=3,xlab='TNR',ylab='Recall')
    }

    SENS<-NULL
    threshold<-NULL
    COORS<-coords(RES,'all',ret = c('threshold','ppv','sensitivity','specificity'))
    if(length(FDRth)>0){


        FDR5percTh<-max(COORS['threshold',which(COORS['ppv',]>=(1-FDRth))])

        threshold<-COORS['threshold',min(which(COORS['threshold',]<=FDR5percTh))]

        SENS<-COORS['sensitivity',min(which(COORS['threshold',]<=FDR5percTh))]
        SPEC<-COORS['specificity',min(which(COORS['threshold',]<=FDR5percTh))]
        if(display){
            abline(h=SENS,lty=2)
        }
    }

    if(display){
        if(length(SENS)==0){
            legend('bottomright',paste('AUC = ',format(RES$auc,digits=3)),bty = 'n')
        }else{
            legend('bottomright',c(paste('Recall ',100*FDRth,'%FDR = ',format(SENS,digits=3),sep=''),
                                   paste('AUC = ',format(RES$auc,digits=3))),bty = 'n')
        }

    }


    COORS<-t(COORS[c('specificity','sensitivity','threshold'),])
    RES<-list(AUC=RES$auc,Recall=SENS,sigthreshold=threshold,curve=COORS)
    ### threshold, and recall at fixed FDR to be returned
    return(RES)
}
PrRc_Curve<-function(FCsprofile,positives,negatives,display=TRUE,FDRth=NULL){

    FCsprofile<-FCsprofile[intersect(c(positives,negatives),names(FCsprofile))]

    predictions<- -FCsprofile
    observations<-is.element(names(FCsprofile),positives)+0
    names(observations)<-names(predictions)


    prc<-pr.curve(scores.class0 = predictions,weights.class0 = observations,
                  curve = TRUE,sorted = TRUE)

    PRECISION<-prc$curve[,2]
    RECALL<-prc$curve[,1]

    if(display){
        plot(RECALL,PRECISION,col='blue',lwd=3,xlab='Recall',ylab='Precision',type='l',xlim=c(0,1),ylim=c(0,1))
    }

    SENS<-NULL
    threshold<-NULL
    if(length(FDRth)>0){

        FDR5percTh<- -prc$curve[min(which(prc$curve[,2]>= 1-FDRth)),3]
        SENS<- prc$curve[min(which(prc$curve[,2]>= 1-FDRth)),1]
        threshold<-FDR5percTh
        if(display){
            abline(h=1-FDRth,lty=2)

            abline(v=SENS,lty=1)
        }
    }

    if(display){
        if(length(SENS)==0){
            legend('bottomleft',paste('AUC = ',format(prc$auc.integral,digits=3)),bty = 'n')
        }else{
            legend('bottomleft',c(paste('Recall ',100*FDRth,'%FDR = ',format(SENS,digits=3),sep=''),
                                  paste('AUC = ',format(prc$auc.integral,digits=3))),bty = 'n')
        }

        abline(h=sum(observations)/length(observations))
    }
    #
    curve<-prc$curve
    colnames(curve)<-c('recall','precision','threshold')
    RES<-list(AUC=prc$auc.integral,Recall=SENS,sigthreshold=threshold,curve=curve)
    # ### threshold, and recall at fixed FDR to be returned
    return(RES)
}