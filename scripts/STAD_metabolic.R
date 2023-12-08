
if (T) {
  dir.create("scripts")
  dir.create("files")
  dir.create("PDFs")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
  dir.create("results/pdf",recursive = T)
}
library(reshape2)
library(ggpubr)
library(ggsci)
library(maftools)
library(tidyr)
library(pheatmap)
library(clusterProfiler)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(corrplot)
library(survminer)
library(survival)
options(stringsAsFactors = F)
source('mg_base.R')
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5)
    )
}
sig_boxplot_t<-function(dat,leg,ylab,xlab='',palette=ggsci::pal_lancet()(10)[3:4]){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="t.test",label = "p.signif")+
    ylab(ylab)+xlab(xlab)+labs(color=leg)
  return(pp)
}
sig_boxplot_w<-function(dat,leg,ylab,xlab='',palette=ggsci::pal_lancet()(10)[3:4]){
  dat=na.omit(dat)
  colnames(dat)=c('group','gene')
  dat=dat[order(dat$group),]
  all.combn=combn(as.character(unique(dat$group)),2)
  my_comparisons=lapply(seq_len(ncol(all.combn)), function(i) all.combn[,i])
  pp=ggboxplot(dat, 
               x='group', y='gene', color = 'group',
               palette =  palette,
               short.panel.labs = T,outlier.shape = NA)+
    stat_compare_means(comparisons=my_comparisons,method="wilcox.test",label = "p.signif")+
    ylab(ylab)+xlab(xlab)+labs(color=leg)
  return(pp)
}
mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  if(jitter){
    if(is.null(point_size)){
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
    }else{
      p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
    }
  }
  
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  #最佳截断
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
bioForest=function(rt=null,col){
  #读取输入文件
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #输出图形
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #绘制右边森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 未过滤任何基因
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 降序排列
  return(geneList)
}
mg_nomogram=function(clinical_riskscore,
                     os,
                     status,
                     title='Nomogram',
                     quick=T,
                     mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3#对观测2的六个指标在列线图上进行计分展示
  #,observation=pbc[2,] #也可以不展示
  #预测3年和5年的死亡风险，此处单位是day
  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE #cox回归中需要TRUE
  #              ,showP = T #是否展示统计学差异
  #              ,droplines = F#观测2示例计分是否画线
  #,colors = mg_colors[1:3] #用前面自己定义的颜色
  #,rank="decreasing") #根据统计学差异的显著性进行变量的排序
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) #展示观测的可信区间
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=crbind2DataFrame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

#######00.数据准备###############
genecode=read.delim('GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]

#TCGA####
#临床数据
tcga_cli<-read.delim('origin_datas/TCGA/Merge_STAD_clinical.txt',sep='\t',header = T)
colnames(tcga_cli)[1:20]
tcga_cli=tcga_cli[,c("A0_Samples","A17_Age","A18_Sex","A3_T","A4_N","A5_M","A6_Stage","A7_Grade","A2_Event","A1_OS")]
colnames(tcga_cli)<-c('Samples','Age','Gender','T.stage','N.stage','M.stage','Stage','Grade','OS','OS.time')
tcga_cli$Samples=paste0(tcga_cli$Samples,"-01")
rownames(tcga_cli)<-tcga_cli$Samples
tcga_cli$OS.time
tcga_cli=tcga_cli %>% drop_na(OS.time)
tcga_cli=tcga_cli[tcga_cli$OS.time>0,]
table(tcga_cli$OS)
tcga_cli$OS[tcga_cli$OS=='Alive']<-0
tcga_cli$OS[tcga_cli$OS=='Dead']<-1
tcga_cli$Status<-ifelse(tcga_cli$OS==0,'Alive','Dead')	
tcga_cli$Age[tcga_cli$Age=='Not Available']<-NA
fivenum(as.numeric(tcga_cli$Age),na.rm = T)
tcga_cli$Age1<-ifelse(tcga_cli$Age>67,'>67','<=67')	
table(tcga_cli$Age1)
tcga_cli=crbind2DataFrame(tcga_cli)
head(tcga_cli)
dim(tcga_cli)
#414  12
table(tcga_cli$Stage)
tcga_cli$Stage <- gsub('[ABC]', '', tcga_cli$Stage)
tcga_cli$Stage <- gsub('Stage ', '', tcga_cli$Stage)
tcga_cli$Stage[tcga_cli$Stage=='']<-NA
table(tcga_cli$Grade)
tcga_cli$Grade[tcga_cli$Grade=='GX']<-NA
table(tcga_cli$T.stage)
tcga_cli$T.stage <- gsub('[abc]', '', tcga_cli$T.stage)
tcga_cli$T.stage[tcga_cli$T.stage=='TX']<-NA
table(tcga_cli$N.stage)
tcga_cli$N.stage <- gsub('[abc]', '', tcga_cli$N.stage)
tcga_cli$N.stage[tcga_cli$N.stage=='NX' | tcga_cli$N.stage=='']<-NA
table(tcga_cli$M.stage)
tcga_cli$M.stage[tcga_cli$M.stage=='MX']<-NA
head(tcga_cli)

#表达谱
tcga_data<-read.delim('origin_datas/TCGA/Merge_TCGA-STAD_TPM.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga_data[1:4,1:4]
table(substr(colnames(tcga_data),14,15))
sample_T=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==1)]#肿瘤样本
sample_T=intersect(sample_T,tcga_cli$Samples)
length(sample_T)#350
sample_N=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==11)]#正常样本
length(sample_N)#32
range(tcga_data)
tcga_tpm_log=log2(tcga_data[intersect(rownames(tcga_data),mrna_genecode$SYMBOL),c(sample_T,sample_N)]+1)
tcga_tpm_log_T=tcga_tpm_log[,sample_T]
dim(tcga_tpm_log_T)
# 19160   350
table(substr(colnames(tcga_tpm_log_T),14,15))
tcga_cli=tcga_cli[intersect(colnames(tcga_tpm_log_T),tcga_cli$Samples),]
dim(tcga_cli)
#350   12


# GSE66229##########
load('origin_datas/GEO/GSE66229_cli_exp.RData')

####01.聚类####################
dir.create('results/01.Cluster')
meta.gene=readMatrix('origin_datas/Metabolism_associated_genes.txt',row = F,header = T)
meta.gene=as.character(meta.gene$`Gene Symbol`)
length(meta.gene)
#2752

meta.exp=tcga_tpm_log_T[intersect(meta.gene,rownames(tcga_tpm_log_T)),]
meta.exp.fit=meta.exp[rowSums(meta.exp>1)>= ncol(tcga_tpm_log_T)/2,]
dim(meta.exp.fit)
#1807  350


tcga.meta.cox=cox_batch(dat=tcga_tpm_log_T[intersect(rownames(meta.exp.fit),rownames(tcga_tpm_log_T)),tcga_cli$Samples],
                        time=tcga_cli$OS.time,event=tcga_cli$OS)
table(tcga.meta.cox$p.value<0.05)

tcga.meta.cox.fit=tcga.meta.cox[tcga.meta.cox$p.value<0.05,]
tcga.meta.cox.fit$type=ifelse(tcga.meta.cox.fit$HR>1,'Risk','Protective')
table(tcga.meta.cox.fit$type)
writeMatrix(tcga.meta.cox.fit,outpath = 'results/01.Cluster/tcga.meta.cox.txt')
###########差异基因
tcga_type=data.frame(Samples=c(sample_T,sample_N),Type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples

tcga.limma=mg_limma_DEG(exp = tcga_tpm_log[,tcga_type$Samples],group = tcga_type$Type,ulab = 'Tumor',dlab = 'Normal')
tcga.limma$Summary
tcga.degs=rownames(tcga.limma$DEG[which(abs(tcga.limma$DEG$logFC)>log2(1.5) & tcga.limma$DEG$adj.P.Val<0.05),])
length(tcga.degs)
writeMatrix(dat =round(tcga.limma$DEG[which(abs(tcga.limma$DEG$logFC)>log2(1.5) & tcga.limma$DEG$adj.P.Val<0.05),],3),
            outpath = 'results/01.Cluster/tcga.degs.txt')


##############预后显著的差异表达代谢相关基因
pdf('results/01.Cluster/Fig1a.pdf',height = 5,width = 7,onefile = F)
mg_venn_plot(list(Metabolic_gene_cox=rownames(tcga.meta.cox.fit),TCGA_DEGs=tcga.degs))
dev.off()


#############一致性聚类
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[1]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[2]
consen_gene=intersect(rownames(tcga.meta.cox.fit),tcga.degs)
length(consen_gene)
########TCGA###########
tcga_consen_data=as.matrix(tcga_tpm_log_T[intersect(consen_gene,rownames(tcga_tpm_log_T)),])
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = T
                                           , seed = 123456)
k=2
tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Cluster=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
writeMatrix(tcga.subtype,row=F,header=T,outpath = 'results/01.Cluster/tcga.subtype.txt')

cluster.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,
                        data = data.frame(OS.time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                                          , OS = tcga_cli[rownames(tcga.subtype),]$OS
                                          , Cluster=tcga.subtype$Cluster)),
           data=data.frame(OS.time = tcga_cli[rownames(tcga.subtype),]$OS.time/365
                           , OS = tcga_cli[rownames(tcga.subtype),]$OS
                           , Cluster=tcga.subtype$Cluster),
           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
           title='TCGA-LIHC',ggtheme=custom_theme(),
           linetype = c("solid", "dashed","strata")[1],
           palette = cluster.color,
           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
           #legend = c(0.8,0.75), # 指定图例位置
           legend.title = "",
           legend.labs = c("C1","C2"))
fig1a=mg_merge_plot(cluster.km$plot,cluster.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')



######GSE66229(**)###########
GSE66229_consen_data=as.matrix(GSE66229_exp[intersect(consen_gene,rownames(GSE66229_exp)),])
GSE66229_consen_data=sweep(GSE66229_consen_data,1,apply(GSE66229_consen_data, 1, median))
GSE66229_consen_data=as.matrix(GSE66229_consen_data)
dim(GSE66229_consen_data)
GSE66229_clust_subtype <- ConsensusClusterPlus(GSE66229_consen_data
                                               , maxK = 10, reps = 500, pItem = 0.8
                                               , pFeature = 1
                                               , title = "GSE66229_subtype"
                                               , clusterAlg = clusterAlg_name
                                               , distance = distance_name
                                               , plot = "pdf"
                                               , writeTable = T
                                               , seed = 123456)
k=2
GSE66229.subtype <- data.frame(Samples = names(GSE66229_clust_subtype[[k]]$consensusClass),
                               Cluster=GSE66229_clust_subtype[[k]]$consensusClass)
GSE66229.subtype$Cluster=paste0('C',GSE66229.subtype$Cluster)
table(GSE66229.subtype$Cluster)
writeMatrix(GSE66229.subtype,row=F,header=T,outpath = 'results/01.Cluster/GSE66229.subtype.txt')
cluster.km1=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,
                                   data = data.frame(OS.time = GSE66229_cli[rownames(GSE66229.subtype),]$OS.time/365
                                                     , OS = GSE66229_cli[rownames(GSE66229.subtype),]$OS
                                                     , Cluster=GSE66229.subtype$Cluster)),
                      data=data.frame(OS.time = GSE66229_cli[rownames(GSE66229.subtype),]$OS.time/365
                                      , OS = GSE66229_cli[rownames(GSE66229.subtype),]$OS
                                      , Cluster=GSE66229.subtype$Cluster),
                      conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                      title='GSE66229',ggtheme=custom_theme(),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = cluster.color,
                      legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                      #legend = c(0.8,0.75), # 指定图例位置
                      legend.title = "",
                      legend.labs = c("C1","C2"))
fig1b=mg_merge_plot(cluster.km1$plot,cluster.km1$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')

fig1ab=mg_merge_plot(fig1a,fig1b,labels = c('D','E'))
savePDF('results/01.Cluster/Fig1b.pdf',fig1ab,height  = 6,width  = 12)


####02.亚型间临床特征################
dir.create('results/02.subtype.cli')
tcga.subtype.cli=merge(tcga_cli,tcga.subtype,by='Samples')
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples
head(tcga.subtype.cli)

cli_anno=tcga.subtype.cli[order(tcga.subtype.cli$Cluster),c(13,11,3:8,12)]
cli_anno$N.stage[which(cli_anno$N.stage=='N1'|cli_anno$N.stage=='N2'|cli_anno$N.stage=='N3')]<-'N1+N2+N3'

chisq1=chisq.test(table(tcga.subtype.cli$Cluster,tcga.subtype.cli$Age1))
chisq2=chisq.test(table(tcga.subtype.cli$Cluster,tcga.subtype.cli$Gender))
chisq3=chisq.test(table(tcga.subtype.cli$Cluster,tcga.subtype.cli$T.stage))
chisq4=chisq.test(table(cli_anno$Cluster,cli_anno$N.stage))
chisq5=chisq.test(table(tcga.subtype.cli$Cluster,tcga.subtype.cli$M.stage))
chisq6=chisq.test(table(tcga.subtype.cli$Cluster,tcga.subtype.cli$Stage))
chisq7=chisq.test(table(tcga.subtype.cli$Cluster,tcga.subtype.cli$Grade))
chisq8=chisq.test(table(tcga.subtype.cli$Cluster,tcga.subtype.cli$Status))

cli.chisq.test1=data.frame(Features=c('Age1','Gender','T.stage','N.stage','M.stage','Stage','Grade','Status'),
           pval=round(c(chisq1$p.value,chisq2$p.value,chisq3$p.value,chisq4$p.value,chisq5$p.value,chisq6$p.value,chisq7$p.value,chisq8$p.value),3))
cli.chisq.test1
writeMatrix(cli.chisq.test1,'results/02.subtype.cli/subtype.cli.chisq.test.txt',row = F)
#   Features  pval
# 1     Age1 0.400
# 2   Gender 0.365
# 3  T.stage 0.237
# 4  N.stage 0.024
# 5  M.stage 0.013
# 6    Stage 0.332
# 7    Grade 0.000
# 8   Status 0.006


consen_gene_cox=tcga.meta.cox.fit[consen_gene,]
gene_anno=data.frame(type=consen_gene_cox$type[order(consen_gene_cox$type)])
rownames(gene_anno)=rownames(consen_gene_cox[order(consen_gene_cox$type),])

cluster.color=pal_nejm()(10)[3:4]
names(cluster.color)=c('C1','C2')

age.col=c("#BEBADA","#FB8072")
names(age.col)=c('<=67','>67')

sex.col=c("#BEBADA","#FB8072")
names(sex.col)=c('MALE','FEMALE')

tstage.col=brewer.pal(n = 4,"Set2")
names(tstage.col)=c('T1','T2','T3','T4')

nstage.col=c("#8DD3C7","#FFFFB3")#brewer.pal(n = 4,"Set3")
names(nstage.col)=c('N0','N1+N2+N3')

mstage.col=c("#8DA0CB","#E78AC3")
names(mstage.col)=c('M0','M1')

stage.col=brewer.pal(n = 4,"Accent")
names(stage.col)=c('I','II','III','IV')

grade.col=brewer.pal(n = 3,"Pastel2")
names(grade.col)=c('G1','G2','G3')

status.col=c("#E41A1C","#377EB8")
names(status.col)=c('Alive','Dead')

colnames(cli_anno)
color_anno=list(Cluster=cluster.color,Age1=age.col,Gender=sex.col,T.stage=tstage.col,N.stage=nstage.col,M.stage=mstage.col,
                Stage=stage.col,Grade=grade.col,Status=status.col)
pdf('results/02.subtype.cli/Fig2.pdf',height = 12,width = 10)
pheatmap(tcga_tpm_log_T[rownames(gene_anno),rownames(cli_anno)],
         scale = 'row',
         color =  colorRampPalette(c("navy", "white", "red"))(200),
         main="Heatmap", # 设置图形标题
         annotation_col = cli_anno,
         annotation_row = gene_anno,
         annotation_colors = color_anno,
         display_numbers = F, # 热图上显示数值
         cluster_cols = F, # 去掉横向、纵向聚类
         cluster_rows = F,
         show_rownames = T, #去掉横、纵坐标id
         show_colnames = F,
         fontsize_row = 12, # 分别设置横向和纵向字体大小
         fontsize_col = 16) 
dev.off()



#######03.代谢特征############
dir.create('results/03.meta.sign')
load('results/tcga.metabo.ssgsea.RData')
tcga.metabo.ssgsea[1:5,1:5]
writeMatrix(tcga.metabo.ssgsea,'results/03.meta.sign/tcga.metabo.ssgsea.txt')
colnames(tcga.metabo.ssgsea)=gsub(' ','_',colnames(tcga.metabo.ssgsea))
colnames(tcga.metabo.ssgsea)=gsub(',','_',colnames(tcga.metabo.ssgsea))
colnames(tcga.metabo.ssgsea)=gsub('-','_',colnames(tcga.metabo.ssgsea))
tcga.pathway.score=t(tcga.metabo.ssgsea)
tcga.pathway.score=crbind2DataFrame(tcga.pathway.score)
diff_pathway<-function(dat,group){
  dat=data.frame(cluster=group,t(dat))
  gr=c('C1','C2')
  dat1=dat[dat$cluster==gr[1],-1]
  dat2=dat[dat$cluster==gr[2],-1]
  pathway=unique(colnames(dat)[-1])
  p_vale=data.frame()
  for (i in pathway){
    dd1=wilcox.test(dat1[,i],dat2[,i])$p.value
    p_vale=rbind(p_vale,data.frame(pathway=i,p.value=dd1))
  }
  return(p_vale)
}
tcga.pathway.score=crbind2DataFrame(tcga.pathway.score)
tcga.pahtway.diff<-diff_pathway(dat=tcga.pathway.score[,tcga.subtype.cli$Samples],group=tcga.subtype.cli$Cluster)
tcga.pahtway.diff$lab<-ifelse(tcga.pahtway.diff$p.value<0.0001,'****',
                              ifelse(tcga.pahtway.diff$p.value<0.001,'***',
                                     ifelse(tcga.pahtway.diff$p.value<0.01,'**',
                                            ifelse(tcga.pahtway.diff$p.value<0.05,'*',''))))
head(tcga.pahtway.diff)
cut_off=0.05
table(tcga.pahtway.diff$p.value<cut_off)
# FALSE  TRUE 
# 24    89 
diff.pathway=tcga.pahtway.diff$pathway[which(tcga.pahtway.diff$p.value<cut_off)]


library(pheatmap)
anno_col=data.frame(Cluster=tcga.subtype.cli$Cluster)
rownames(anno_col)=tcga.subtype.cli$Samples
anno_col=anno_col[order(anno_col$Cluster),,drop=F]

pdf('results/03.meta.sign/Fig3.pdf',height =12,width =12)
#pal_nejm()(10)[3:4]
pheatmap(as.matrix(tcga.pathway.score[diff.pathway,rownames(anno_col)]),
         name = 'ssGSEA score',scale = 'row',
         clustering_method = c('ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average','mcquitty', 'median', 'centroid')[7],
         clustering_distance_rows = "correlation",
         show_rownames = T,cluster_rows = T,
         cutree_rows = 2,
         #gaps_row = c(80, 9),
         #border_color = "black",
         show_colnames = F,cluster_cols = F,
         #cutree_cols = 2,
         gaps_col  = 185,
         annotation_col = anno_col,breaks = unique(c(seq(-2,2,length=100))),
         annotation_colors = list(Cluster=c('C1'='#E18727FF',
                                           'C2'='#20854EFF')),
         color = colorRampPalette(c('#003399', "white", '#CC0033'))(100))
dev.off()

#####04.亚型间免疫/通路#################
dir.create('results/04.subtype.immu.path')
load('results/tcga.ciber.RData')
writeMatrix(tcga.ciber,'results/04.subtype.immu.path/tcga.cibersort.txt')
mg_PlotMutiBoxplot(data = tcga.ciber[tcga.subtype.cli$Samples,1:22],
                   group = tcga.subtype.cli$Cluster,
                   test_method = c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[2],
                   ylab = 'Score',group_cols = ggsci::pal_nejm()(10)[3:4],
                   legend.pos = 'top',
                   add = 'boxplot')
p.dat1=cbind(tcga.ciber[tcga.subtype.cli$Samples,1:22],Cluster=tcga.subtype$Cluster)
p.dat1=melt(p.dat1)
head(p.dat1)
p1=p.dat1 %>% 
  ggplot(aes(x=variable, y=value,fill=Cluster)) +
  geom_boxplot()+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values =pal_nejm()(10)[3:4])+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Cluster") +
  #theme_light()+
  theme_classic()+
  theme(legend.position = 'top',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 将图表标题居中


load('results/tcga.immu.est.RData')
writeMatrix(tcga.est,'results/04.subtype.immu.path/tcga.est.txt')
mg_PlotMutiBoxplot(data = tcga.est[tcga.subtype.cli$Samples,],
                   group = tcga.subtype.cli$Cluster,
                   test_method = c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[2],
                   ylab = 'Score',group_cols = ggsci::pal_nejm()(10)[3:4],
                   legend.pos = 'top',
                   add = 'boxplot')
p.dat2=cbind(tcga.est[tcga.subtype.cli$Samples,],tcga.subtype)
head(p.dat2)
p2=p.dat2 %>% 
  ggplot(aes(x=Cluster, y=StromalScore,fill=Cluster)) +
  geom_violin()+  #根据Ancestry的不同因子使用不同颜色，其实用R默认颜色也不错，这里只是展示一下如何提取喜欢的图片颜色。
  scale_fill_manual(values = pal_nejm()(10)[3:4])+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'wilcox.test')+
  theme_classic()+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 15))

p3=p.dat2 %>% 
  ggplot(aes(x=Cluster, y=ImmuneScore,fill=Cluster)) +
  geom_violin()+  #根据Ancestry的不同因子使用不同颜色，其实用R默认颜色也不错，这里只是展示一下如何提取喜欢的图片颜色。
  scale_fill_manual(values = pal_nejm()(10)[3:4])+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  theme_classic(base_size = 20)+
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'wilcox.test')+
  theme(legend.position = 'none',axis.text = element_text(color = 'black'),
        axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
        axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 15))

fig4a=mg_merge_plot(mg_merge_plot(p2,p3,ncol = 1,nrow = 2,common.legend = T,labels = c('A','B')),
                    p1,ncol=2,nrow=1,labels = c('','C'),widths = c(1,3),common.legend = T)
savePDF('results/04.subtype.immu.path/Fig4a.pdf',fig4a,height = 5,width = 10)


load('results/tcga.immu.ssgsea.RData')
writeMatrix(tcga.immu.ssgsea,'results/04.subtype.immu.path/tcga.immu.ssgsea.txt')
p.dat4=cbind(tcga.immu.ssgsea[tcga.subtype.cli$Samples,],Cluster=tcga.subtype$Cluster)
p.dat4=crbind2DataFrame(p.dat4)
p.dat4=melt(p.dat4)
head(p.dat4)
fig4b=p.dat4 %>% 
  ggplot(aes(x=variable, y=value,fill=Cluster)) +
  geom_boxplot()+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values =pal_nejm()(10)[3:4])+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Cluster",title = 'ssGSEA') +
  #theme_light()+
  theme_classic()+
  theme(legend.position = 'top',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 将图表标题居中
fig4b

#########MCPcounter
writeMatrix(tcga.mcp,'results/04.subtype.immu.path/tcga.mcp.txt')
p.dat5=cbind(tcga.mcp[tcga.subtype.cli$Samples,],Cluster=tcga.subtype$Cluster)
p.dat5=crbind2DataFrame(p.dat5)
p.dat5=melt(p.dat5)
head(p.dat5)
fig4c=p.dat5 %>% 
  ggplot(aes(x=variable, y=value,fill=Cluster)) +
  geom_boxplot()+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values =pal_nejm()(10)[3:4])+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Cluster",title = 'MCP-Count') +
  #theme_light()+
  theme_classic()+
  theme(legend.position = 'top',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 将图表标题居中
fig4c

fig4bc=mg_merge_plot(fig4b,fig4c,ncol = 1,nrow = 2,common.legend = T,labels = c('A','B'))
savePDF('results/04.subtype.immu.path/FigS1.pdf',fig4bc,height = 12,width = 10)



#############mdp-count and ssgsea heatmap
cluster.anno=data.frame(Cluster=tcga.subtype$Cluster[order(tcga.subtype$Cluster)])
rownames(cluster.anno)=tcga.subtype$Samples[order(tcga.subtype$Cluster)]

tcga.immu.dat.sub=t(cbind(tcga.mcp[tcga.subtype$Samples,],tcga.immu.ssgsea[tcga.subtype$Samples,]))
immu.type=data.frame(type=rep(c('MCP-Count','ssGSEA'),c(ncol(tcga.mcp),ncol(tcga.immu.ssgsea))))
rownames(immu.type)=rownames(tcga.immu.dat.sub)

pdf('results/04.subtype.immu.path/Fig4b.pdf',height = 7,width = 11,onefile = F)
pheatmap(t(tcga.immu.dat.sub[rownames(immu.type),rownames(cluster.anno)]),
         scale = 'column',
         #border="white", # 设置边框为白色
         color =  colorRampPalette(c("navy", "white", "red"))(200),
         main="Heatmap", # 设置图形标题
         annotation_col = immu.type,
         annotation_row = cluster.anno,
         annotation_colors = list(Cluster=c('C1'='#E18727FF',
                                            'C2'='#20854EFF')),
         #cutree_cols = 2, #列划为2块，
         #cutree_rows =2, # 行为2块
         #cellwidth = 6,cellheight = 5, # 设置热图方块
         display_numbers = F, # 热图上显示数值
         cluster_cols = F, # 去掉横向、纵向聚类
         cluster_rows = F,
         show_rownames = F, #去掉横、纵坐标id
         show_colnames = T,
         #legend_breaks=c(-5,0,5),
         fontsize_row = 12, # 分别设置横向和纵向字体大小
         fontsize_col = 16) 
dev.off()


####################TIDE###########
# tcga_tide_dat <- t(scale(t(tcga_tpm_log_T),scale = F))
# dim(tcga_tide_dat)
# write.table(tcga_tide_dat,file = 'results/04.subtype.immu.path/tcga_tide_dat.txt',quote = F, sep = '\t')

tcga_tide_res<-read.csv('results/04.subtype.immu.path/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)


tide_sel=c('TIDE','IFNG','Exclusion','Dysfunction','MDSC')
tcga.tide.dat.sub=cbind.data.frame(tcga.subtype.cli[,'Cluster'],tcga_tide_res[tcga.subtype.cli$Samples,tide_sel])
colnames(tcga.tide.dat.sub)[1]='Cluster'
tcga.tide.dat.sub=melt(tcga.tide.dat.sub)
head(tcga.tide.dat.sub)
colnames(tcga.tide.dat.sub)=c('Cluster','type','value')
head(tcga.tide.dat.sub)
p.tide1=tcga.tide.dat.sub %>% 
  ggplot(aes(x=Cluster, y=value,fill=Cluster)) +
  geom_boxplot()+ facet_wrap(~type,scales = 'free',nrow = 1,ncol = 5)+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values = ggsci::pal_nejm()(10)[3:4])+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Cluster") +
  theme_light()+
  theme(legend.position = 'none',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, hjust = 0.5)) # 将图表标题居中
p.tide1
savePDF('results/04.subtype.immu.path/Fig4d.pdf',p.tide1,height = 4,width = 10)



##############通路GSEA
mg_RunGSEA_use=function(mod=c('exp_group','exp_gene','rank')[1],exp_Path=NULL, sample_group_path=NULL, outFolder=NULL,gene=NULL, column=NULL,lower=50,upper=50, gmt_Path=c("KEGG",'GO_BP','GO_CC','GO_MF','reactome','HALLMARK','TF')[1],plot_svg=FALSE,top=10,min=5,max=5000,outLog=T){
  
  if(is.null(exp_Path)|is.null(mod)|is.null(outFolder)|is.null(gmt_Path)){
    return(NULL)
  }
  if(plot_svg){
    svg='true'
  }else{
    svg='false'
  }
  if(gmt_Path=='KEGG'){
    #gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c2.cp.kegg.v7.0.symbols.gmt')
    
  }
  else if(gmt_Path=='GO_BP'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.bp.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='GO_CC'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.cc.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='GO_MF'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c5.mf.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='reactome'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c2.cp.reactome.v7.0.symbols.gmt')
  }
  else if(gmt_Path=='HALLMARK'){
    #gmt_Path=paste0(MG_Grobal_baseFolder,'/source/h.all.v7.0.symbols.gmt')
    
  }
  else if(gmt_Path=='TF'){
    gmt_Path=paste0(MG_Grobal_baseFolder,'/source/c3.tft.v7.0.symbols.gmt')
  }
  if(file.exists(paste0(getwd(),'/',outFolder))){
    outFolder=paste0(getwd(),'/',outFolder)
  }else if(!file.exists(outFolder)){
    dir.create(outFolder)
    if(file.exists(paste0(getwd(),'/',outFolder))){
      outFolder=paste0(getwd(),'/',outFolder)
    }
  }
  
  if(file.exists(paste0(getwd(),'/',gmt_Path))){
    gmt_Path=paste0(getwd(),'/',gmt_Path)
  }
  if(file.exists(paste0(getwd(),'/',exp_Path))){
    exp_Path=paste0(getwd(),'/',exp_Path)
  }
  
  command=NULL
  if(mod=='exp_group'){
    if(!is.null(exp_Path)&!is.null(sample_group_path)&!is.null(outFolder)){
      if(file.exists(paste0(getwd(),'/',sample_group_path))){
        sample_group_path=paste0(getwd(),'/',sample_group_path)
      }
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar exp_group '
                     ,exp_Path,' ',sample_group_path,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max) 
    }
  }else if(mod=='exp_gene'){
    if(!is.null(exp_Path)&!is.null(gene)&!is.null(outFolder)){
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar exp_gene '
                     ,exp_Path,' ',gene,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max,' ',lower,' ',upper) 
    }
  }else if(mod=='rank'){
    if(!is.null(exp_Path)&!is.null(column)&!is.null(outFolder)){
      command=paste0(MG_Grobal_baseFolder,'/jre/GSEA/MG_GSEA.jar rank '
                     ,exp_Path,' ',column-1,' ',outFolder,' ',gmt_Path,' ',svg,' ',top,' ',min,' ',max) 
    }
  }
  
  if(!is.null(command)){
    if(MG_Grobal_System=='win'){
      command=paste0(MG_Grobal_baseFolder,'/jre/bin/java -jar ',command)
    }else{
      command=paste0('java -jar ',command)
    }
    print(paste0('RunGSEA CMD:',command))
    logs=system(command, intern = !outLog, 
                ignore.stdout = FALSE, ignore.stderr = FALSE, 
                wait = TRUE, input = NULL, show.output.on.console = TRUE, 
                minimized = FALSE, invisible = TRUE)
    if(outLog){
      if(logs==0){
        print('Run GSEA succ')
      }else{
        print('Run GSEA error')
      }
    }else{
      print('Runed GSEA')
      print(logs)
      logs=logs[grep('######/',logs)]
      if(length(logs)==1){
        logs=unlist(strsplit(logs[1],'/'))
        if(length(logs)>1){
          return(logs[2:length(logs)])
        }
      }
    }
  }
  return(NULL)
}
dir.create('results/04.subtype.immu.path/GSEA')
writeMatrix(dat = tcga_tpm_log_T[,tcga.subtype$Samples],outpath = 'results/04.subtype.immu.path/tcga_tpm_log_T.txt',row = T,header = T)
writeMatrix(dat = tcga.subtype,outpath = 'results/04.subtype.immu.path/tcga.subtype.use.txt',row = F,header = T)
# mg_RunGSEA_use(mod = 'exp_group',exp_Path = 'results/04.subtype.immu.path/tcga_tpm_log_T.txt'
#            ,sample_group_path = 'results/04.subtype.immu.path/tcga.subtype.use.txt'
#            ,outFolder = 'results/04.subtype.immu.path/GSEA'
#            ,gmt_Path = 'HALLMARK',outLog=F)
tcga_GSEA=parseGSEAResult('results/04.subtype.immu.path/GSEA/my_analysis.Gsea.1666661900437/')#######HALLMARK
head(tcga_GSEA$EnrichTable)
tcga_GSEA.res<-tcga_GSEA$EnrichTable
table(tcga_GSEA.res$NES>0)
tcga_GSEA.res=tcga_GSEA.res[,c("Term","NES","NP","FDR")]
tcga_GSEA.res$group<-ifelse(tcga_GSEA.res$NES>0,'C1','C2')
tcga_GSEA.res$NES=abs(tcga_GSEA.res$NES)
tcga_GSEA.res=tcga_GSEA.res[order(tcga_GSEA.res$group,-tcga_GSEA.res$NES,decreasing = T),]
tcga_GSEA.res$Term=gsub('HALLMARK_','',tcga_GSEA.res$Term)
tcga_GSEA.res$Term=factor(tcga_GSEA.res$Term,levels = tcga_GSEA.res$Term)
tcga_GSEA.res=tcga_GSEA.res[tcga_GSEA.res$FDR<0.25 ,]
dim(tcga_GSEA.res)
writeMatrix(tcga_GSEA.res,outpath = 'results/04.subtype.immu.path/tcga_GSEA.res.txt')
library(ggplot2)
pdf('results/04.subtype.immu.path/Fig4e.pdf',height = 15,width = 10)
ggplot(tcga_GSEA.res,aes(x=Term,y=group,size=tcga_GSEA.res$NES,color=tcga_GSEA.res$FDR))+
  geom_point()+coord_flip()+xlab('')+ylab('')+
  scale_color_gradient(low = 'red',high = 'blue')+
  theme_bw()+theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))
dev.off()

####05.亚型间突变特征################
dir.create('results/05.subtype.mut')
# #不同分子亚型之间基因突变的差异
tcga.maf=getTCGAMAFByCode('STAD')
tcga.subtype.use=tcga.subtype
colnames(tcga.subtype.use)[1]='Tumor_Sample_Barcode'
tcga.subtype.use$Tumor_Sample_Barcode=substr(tcga.subtype.use$Tumor_Sample_Barcode,1,12)
tcga.subtype.use.C1=tcga.subtype.use[which(tcga.subtype.use$Cluster=='C1'),]
tcga.subtype.use.C2=tcga.subtype.use[which(tcga.subtype.use$Cluster=='C2'),]
tcga.subtype.use.all=tcga.subtype.use[order(tcga.subtype.use$Cluster),]
write.table(tcga.subtype.use.C1,file='results/05.subtype.mut/tcga.subtype.c1.txt')
write.table(tcga.subtype.use.C2,file='results/05.subtype.mut/tcga.subtype.c2.txt')
# write.table(tcga.subtype.use.all,file='results/05.subtype.mut/tcga.subtype1.txt')
tcga.maf1=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.subtype.use.C1$Tumor_Sample_Barcode))
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,clinicalData = 'results/05.subtype.mut/tcga.subtype.c1.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.subtype.use.C2$Tumor_Sample_Barcode))
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,clinicalData = 'results/05.subtype.mut/tcga.subtype.c2.txt')
tcga.maf2@clinical.data


tcga.mut.dat <- tcga.maf
tcga.mut.dat <- as.data.frame(tcga.mut.dat@data)
tcga.mut.dat <- tcga.mut.dat[, c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]

tcga.mut.dat$Variant_Classification <- 1
tcga.mut.dat <- reshape2::dcast(data = tcga.mut.dat, Hugo_Symbol ~ Tumor_Sample_Barcode)
class(tcga.mut.dat)
rownames(tcga.mut.dat) <- tcga.mut.dat$Hugo_Symbol
tcga.mut.dat <- tcga.mut.dat[, -1]

colnames(tcga.mut.dat) <- paste0(colnames(tcga.mut.dat), '-01')
mut.samples <- intersect(colnames(tcga.mut.dat), tcga.subtype.cli$Samples)


tcga.mut.dat <- tcga.mut.dat[, mut.samples]
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples
tcga_mut_cli <- tcga.subtype.cli[mut.samples, ]

tcga.mut.dat.freq <- as.data.frame(rowSums(tcga.mut.dat))
colnames(tcga.mut.dat.freq) <- 'Freq'
tcga.mut.dat.freq$Genes <- rownames(tcga.mut.dat.freq)
library(dplyr)
str(tcga.mut.dat.freq)
head(tcga.mut.dat.freq)
tcga.mut.dat.freq <- dplyr::arrange(tcga.mut.dat.freq, desc(Freq))
head(tcga.mut.dat.freq)
dim(tcga.mut.dat.freq)

mut.genes <- rownames(tcga.mut.dat.freq)[tcga.mut.dat.freq$Freq > 3]
length(mut.genes)

tcga.mut.dat <- ifelse(tcga.mut.dat > 0, 'Mutant', 'WildType')
dim(tcga.mut.dat)

mut.res <- data.frame(C1 = NA,
                      C2 = NA)
mut.p <- c()
for (ge in mut.genes) {
  print(ge)
  tmp <- table(tcga.mut.dat[ge, ], tcga_mut_cli$Cluster)
  pvalue <- fisher.test(tmp)
  mut.p <- c(mut.p, pvalue$p.value)
  mut.res <- rbind(mut.res, tmp[1, ])
}
mut.res <- na.omit(mut.res)
rownames(mut.res) <- mut.genes
class(mut.res)
mut.res$P.value <- mut.p

table(mut.res$P.value < 0.05)
# FALSE  TRUE 
# 9059  1126
mut.res.filtered <- mut.res[which(mut.res$P.value < 0.05), ]
mut.res.filtered
dim(mut.res.filtered)
writeMatrix(mut.res.filtered,'results/05.subtype.mut/subtype.mut.gene.txt')
writeMatrix(mut.res.filtered,'P_20221018_STAD_metabolic/Files/subtype.mut.gene.txt')
ggsci::pal_nejm()(10)[3:4]
pdf('results/05.subtype.mut/Fig5a.pdf',height = 6,width = 6)
oncoplot(maf=tcga.maf1,clinicalFeatures = 'Cluster',
         genes = rownames(mut.res.filtered)[1:15],
         sortByAnnotation = T,
         annotationColor = list(Cluster=c(C1='#E18727FF')))
dev.off()
pdf('results/05.subtype.mut/Fig5b.pdf',height = 6,width = 6)
oncoplot(maf=tcga.maf2,clinicalFeatures = 'Cluster',
         genes = rownames(mut.res.filtered)[1:15],
         sortByAnnotation = T,
         annotationColor = list(Cluster=c(C2='#20854EFF')))
dev.off()
# mg_group_oncoplot_byMAF(maf=tcga.maf
#                         ,group=tcga.subtype.use.all
#                         ,genes = rownames(mut.res.filtered)[1:15]
#                         ,annotationColor='nejm')


##############分子特征
tcga.sub.all=readMatrix(paste0(MG_Grobal_baseFolder,'/source/PMC5982584_supplement_2.txt'))
stad.sub.all<-tcga.sub.all[tcga.sub.all$`TCGA Study`=='STAD',]
head(stad.sub.all)
rownames(stad.sub.all)=paste0(rownames(stad.sub.all),'-01')
stad.sub.all=stad.sub.all[intersect(tcga.subtype$Samples,row.names(stad.sub.all)),]
dim(stad.sub.all)
#349  63
head(stad.sub.all)
stad.sub.all$Samples=rownames(stad.sub.all)
table(stad.sub.all$`TCGA Subtype`)
#GI.CIN      GI.EBV       GI.GS GI.HM-indel   GI.HM-SNV 
#195          26         42          55           7 
table(stad.sub.all$`Immune Subtype`)
######## C1（伤口愈合），C2（INF-r 占优势），C3（炎症），C4（淋巴细胞耗竭），C5（在免疫学上沉默）和 C6（TGF-beta 占优势）
stad.sub.tcga=merge(stad.sub.all,tcga.subtype.cli,by='Samples')
stad.sub.tcga=merge(stad.sub.tcga,stad.purity[,c('Samples','Cluster','TMB')],by='Samples')

col.selected=c('Aneuploidy Score','Homologous Recombination Defects','Fraction Altered','Number of Segments')
rownames(stad.sub.tcga)=stad.sub.tcga$Samples
tcga.mut.dat.sub=cbind.data.frame(tcga.subtype[,'Cluster'],stad.sub.tcga[tcga.subtype$Samples,c(col.selected,'TMB')])
colnames(tcga.mut.dat.sub)[1]='Cluster'
tcga.mut.dat.sub=melt(tcga.mut.dat.sub)
head(tcga.mut.dat.sub)
colnames(tcga.mut.dat.sub)=c('Cluster','type','value')
head(tcga.mut.dat.sub)

pdf('results/05.subtype.mut/Fig5c.pdf',height = 4,width = 10)
tcga.mut.dat.sub %>% 
  ggplot(aes(x=Cluster, y=value,fill=Cluster)) +
  geom_boxplot()+ facet_wrap(~type,scales = 'free',nrow = 1,ncol = 5)+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values = pal_nejm()(10)[3:4])+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Cluster), label = "p.format", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Cluster") +
  theme_light()+
  theme(legend.position = 'none',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 1)) # 将图表标题居中
dev.off()


###########06.WGCNA##########################
dir.create('results/06.WGCNA')
# cluster gene
library(WGCNA)
allowWGCNAThreads(nThreads = 36)#允许R语言程序最大线程运行
enableWGCNAThreads(nThreads = 36)# 打开多线程

tcga.mads=apply(tcga_tpm_log_T, 1, mad)
tpm_T2=tcga_tpm_log_T[which(tcga.mads>0.5),]
dim(tpm_T2)
# tpm_T2=(2^tpm_T2-1)
tpm_T2=t(tpm_T2)
dim(tpm_T2)
range(tpm_T2)

pdf('results/06.WGCNA/1.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()


tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.25,
                                 minModuleSize=50)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('results/06.WGCNA/2.pdf',height = 6,width = 10)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# 
 writeMatrix(tpm_T2.module$Modules,outpath = 'results/06.WGCNA/tcga.wgcna.module.genes.txt')
pdf('results/06.WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()

#### 模块特征向量聚类
# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('results/06.WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()


############### 模块特征向量与代谢亚型的关系
c1.inds=which(tcga.subtype$Cluster=='C1')
c2.inds=which(tcga.subtype$Cluster=='C2')
gtype.group=rep('C',length(tcga.subtype$Cluster))
gtype.group[c1.inds]<-'C1'
gtype.group[c2.inds]<-'C2'

design <- model.matrix(~0+factor(gtype.group))
colnames(design)=c('C1','C2')

spms=design

MEs_col<-tpm_T2.module$MEs
dim(MEs_col)

modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('results/06.WGCNA/5.pdf',width = 12,height = 6)
labeledHeatmap(Matrix = data.frame(t(modTraitCor)),
               xLabels = colnames(t(modTraitCor)),
               yLabels = rownames(t(modTraitCor)),
               cex.lab = 1,
               ySymbols = colnames(modTraitCor), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(t(textMatrix)),
               setStdMargins = FALSE,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#基因的模块成员度（module membership）计算
#即各基因表达值与相应模块特征基因的相关性，其衡量了基因在全局网络中的位置
## Notes: signedKME函数的列名是根据输入数据的列名称从第3个字母开始截取，所以如果前面的模块
## 的名称前面的ME如果已经去除，记得在前面加上两个字母
geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))
#各基因表达值与临床表型的相关性分析
geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "turquoise"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
turquoise.gene=names(which(moduleGenes))
length(turquoise.gene)
table(tpm_T2.module$Modules[,'mergedColors'])

pdf('results/06.WGCNA/6.pdf',height = 6,width = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 'C2']),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for C2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
dev.off()

#############模块基因GO KEGG
turquoise.gene.kegg=enrichmentORA(turquoise.gene,maxNum = length(turquoise.gene),
                                  mp_dbs=c('pathway_KEGG',
                                           'geneontology_Biological_Process',
                                           'geneontology_Cellular_Component',
                                           'geneontology_Molecular_Function'))
turquoise.gene.kegg.fit=turquoise.gene.kegg[turquoise.gene.kegg$FDR<0.05,]
dim(turquoise.gene.kegg.fit)
writeMatrix(turquoise.gene.kegg.fit,'results/06.WGCNA/turquoise.module.gene.go.kegg.txt')


pdf('results/06.WGCNA/7.pdf',height = 8,width = 12)
dotplot_batch(enrichmentORA = turquoise.gene.kegg,top = 10,FDR = T,
              dbs = c('pathway_KEGG',
                      'geneontology_Biological_Process',
                      'geneontology_Cellular_Component',
                      'geneontology_Molecular_Function'))
dev.off()

########07.module+nomo###############
dir.create('results/07.module')
pre.genes=turquoise.gene

tcga_model_data=t(tcga_tpm_log_T[pre.genes,tcga_cli$Samples])
colnames(tcga_model_data)=gsub('-','__',colnames(tcga_model_data))
tcga_model_data=merge(data.frame(Samples=tcga_cli$Samples,OS=tcga_cli$OS,OS.time=tcga_cli$OS.time),
                      data.frame(Samples=rownames(tcga_model_data),tcga_model_data),by='Samples')
rownames(tcga_model_data)=tcga_model_data$Samples
tcga_model_data=tcga_model_data[,-1]
tcga_model_data=crbind2DataFrame(tcga_model_data)
dim(tcga_model_data)

########module_select#############
dir.create('files/model_select')

select_gene_zscore <- function(dat1, dat2 = NULL, dat3 = NULL, 
                               a = 1, n = 100, 
                               ratio = 0.5, cut_p = 0.05, 
                               years = c(1,3,5)){
  library(timeROC)
  library(survival)
  library(glmnet)
  library(mosaic)
  sample.index <- data.frame()
  SampleingTime <- c()
  SamplingSeed <- c()
  tra.auc <- c()
  test.auc <- c()
  geo1.auc <- c()
  geo2.auc <- c()
  all.auc <- c()
  tra.p <- c()
  test.p <- c()
  geo1.p <- c()
  geo2.p <- c()
  all.p <- c()
  tcga.dat <- dat1
  geo1.dat <- dat2
  geo2.dat <- dat3
  gene.list <- c()
  GeneNums <- c()
  for (i in a:n) {
    set.seed(i)
    par(mfrow = c(2, 3))
    myd.index <- seq(1,nrow(tcga.dat),1)
    tra.index <- sample(myd.index,size = round(nrow(tcga.dat)*ratio))
    tra.dat <- tcga.dat[tra.index,]
    test.dat <- tcga.dat[-tra.index,]
    write.table(rownames(tra.dat), file = paste0('files/model_select/tra.dat_zscore_', i, '.txt'), sep = '\t', quote = F, row.names = F)
    write.table(rownames(test.dat), file = paste0('files/model_select/test.dat_zscore_', i, '.txt'), sep = '\t', quote = F, row.names = F)
    tra.cox <- t(apply(tra.dat[,3:c(ncol(tra.dat))],2,function(x){
      vl=as.numeric(x)
      tm=tra.dat$OS.time
      ev=tra.dat$OS
      #ev=ifelse(ev=='Alive',0,1)
      dat=data.frame(tm,ev,vl)[which(tm > 0 & !is.na(vl)),]
      return(coxFun(dat))
    }))
    colnames(tra.cox) <- c('p.value','HR','Low 95%CI','High 95%CI')
    tra.cox <- as.data.frame(tra.cox)
    write.table(tra.cox, file = paste0('files/model_select/sig_cox_zscore_', i, '.txt'), sep = '\t', quote = F)
    cut.cox=tra.cox[which(tra.cox[,1]<cut_p),] 
    if (nrow(cut.cox) > 1) {
      print(paste("Processing....: ",i," resampling",sep=""))
      flush.console()
      geneid=rownames(cut.cox)
      geneid=geneid[which(!is.na(geneid))]
      if (length(geneid)>1) {
        sample.index<-c(sample.index,data.frame(tra.index))
        set.seed(i)
        cv.fit <- cv.glmnet(as.matrix(tra.dat[,geneid]), cbind(time=tra.dat$OS.time, 
                                                               status=tra.dat$OS)
                            ,family="cox", nlambda=100, alpha=1) 
        sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
        write.table(sig.coef, file = paste0('files/model_select/sig.coef_zscore_', i, '.txt'), sep = '\t', quote = F)
        if (length(names(sig.coef)>1)) {
          dat1 <- cbind(time=tra.dat$OS.time,
                        status=tra.dat$OS,
                        tra.dat[,match(names(sig.coef),
                                       colnames(tra.dat))])
          fmla <- as.formula(paste0("Surv(time, status) ~"
                                    ,paste0(names(sig.coef),collapse = '+')))
          cox <- coxph(fmla, data =as.data.frame(dat1))
          # lan=coef(cox)
          # print(lan)
          cox1=step(cox, trace = 0)
          lan=coef(cox1)
          write.table(lan, file = paste0('files/model_select/lan_zscore_', i, '.txt'), sep = '\t', quote = F)
          # lan=sig.coef
          final_gene=names(cox1$coefficients)
          # final_gene=names(cox$coefficients)
          GeneNums <-c(GeneNums, length(final_gene))
          
          risk.tra=as.numeric(lan%*%as.matrix(t(tra.dat[,final_gene])))
          risk.tra=zscore(risk.tra)
          ROC.DSST <- timeROC(T=tra.dat$OS.time,
                              delta=tra.dat$OS,
                              marker=risk.tra,
                              cause=1,
                              weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),
                              iid=TRUE)
          tra.auc=c(tra.auc,max(ROC.DSST$AUC))
          trap=plotKMCox(data.frame(tra.dat$OS.time,tra.dat$OS,ifelse(risk.tra>=0,'H','L')))
          tra.p=c(tra.p,trap)
          risk.test=as.numeric(lan%*%as.matrix(t(test.dat[,final_gene])))
          risk.test=zscore(risk.test)
          ROC.DSST1=timeROC(T=test.dat$OS.time,delta=test.dat$OS,marker=risk.test,cause=1,weighting="marginal",
                            times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
          test.auc=c(test.auc,max(ROC.DSST1$AUC))
          testp=plotKMCox(data.frame(test.dat$OS.time,test.dat$OS,ifelse(risk.test>=0,'H','L')))
          test.p=c(test.p,testp)
          risk.all=as.numeric(lan%*%as.matrix(t(tcga.dat[,final_gene])))
          risk.all=zscore(risk.all)
          ROC.DSST3=timeROC(T=tcga.dat$OS.time,delta=tcga.dat$OS,marker=risk.all,cause=1,weighting="marginal",
                            times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
          all.auc=c(all.auc,max(ROC.DSST3$AUC))
          allp=plotKMCox(data.frame(tcga.dat$OS.time,tcga.dat$OS,ifelse(risk.all>=0,'H','L')))
          all.p=c(all.p,allp)
          # final_gene1=as.character(gene.type[final_gene,1])
          
          if (length(intersect(final_gene,colnames(geo1.dat)))==length(final_gene)) {
            risk.geo1=as.numeric(lan%*%as.matrix(t(geo1.dat[,intersect(final_gene,colnames(geo1.dat))])))
            risk.geo1=zscore(risk.geo1)
            ROC.DSST2=timeROC(T=geo1.dat$OS.time,delta=geo1.dat$OS,marker=risk.geo1,cause=1,weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
            geo1.auc=c(geo1.auc,max(ROC.DSST2$AUC))
            geop=plotKMCox(data.frame(geo1.dat$OS.time,geo1.dat$OS,ifelse(risk.geo1>=0,'H','L')))
            geo1.p=c(geo1.p,geop)
          }else{
            geo1.auc=c(geo1.auc,NA)
            geo1.p=c(geo1.p,NA)
          }
          
          if (length(intersect(final_gene,colnames(geo2.dat)))==length(final_gene)) {
            risk.geo2=as.numeric(lan%*%as.matrix(t(geo2.dat[,intersect(final_gene,colnames(geo2.dat))])))
            risk.geo2=zscore(risk.geo2)
            ROC.DSST4=timeROC(T=geo2.dat$OS.time,delta=geo2.dat$OS,marker=risk.geo2,cause=1,weighting="marginal",
                              times=c(365*years[1],365*years[2],365*years[3]),iid=TRUE)
            geo2.auc=c(geo2.auc,max(ROC.DSST4$AUC))
            geop=plotKMCox(data.frame(geo2.dat$OS.time,geo2.dat$OS,ifelse(risk.geo2>=0,'H','L')))
            geo2.p=c(geo2.p,geop)
          }else{
            geo2.auc=c(geo2.auc,NA)
            geo2.p=c(geo2.p,NA)
          }
          print(c('AUC',max(ROC.DSST$AUC),max(ROC.DSST1$AUC),max(ROC.DSST3$AUC)))
          print(c('P',trap,testp,allp))
          print(c('num',length(final_gene)))
          print("........")
          SampleingTime=c(SampleingTime,i)
          gene.list=c(gene.list,as.character(final_gene))
          SamplingSeed=c(SamplingSeed,rep(i,length(final_gene)))
        }
      }
    }
  }
  myd.clustering.df=data.frame("SamplingTime"=SampleingTime,
                               "TrainRiskP"=tra.p,"TrainRiskAUC"=tra.auc,
                               "TestRiskP"=test.p,"TestRiskAUC"=test.auc,
                               "TCGARiskP"=all.p,"TCGARiskAUC"=all.auc,
                               "GEO1RiskP"=geo1.p,"GEO1RiskAUC"=geo1.auc,
                               "GEO2RiskP"=geo2.p,"GEO2RiskAUC"=geo2.auc,
                               "GeneNums"=GeneNums)
  sample.index=as.data.frame(sample.index)
  colnames(sample.index)=paste("seed",SampleingTime,sep="")
  myd.clustering.genes=data.frame("SamplingSeed"=SamplingSeed,"Genes"=gene.list)
  return(list(myd.clustering.df,sample.index,myd.clustering.genes))
}


# myd_exp_resampling <- select_gene_zscore(tcga_model_data,
#                                          #gse21257_model_data,
#                                          a = 1,
#                                          n = 500,
#                                          ratio = 0.7,
#                                          cut_p = 0.01, 
#                                          years = c(1,2,3,4,5))
# myd_exp_resampling[[1]]
write.csv(myd_exp_resampling[[1]],
          file = 'files/tcga-OS-1-500-0.7-0.01.csv',
          quote = F, row.names = F)
num <- 432
tra.samples <- rownames(read.delim(paste0('files/model_select/tra.dat_zscore_',num,'.txt'), 
                                   header = T, row.names = 1, stringsAsFactors = F))
test.samples <- rownames(read.delim(paste0('files/model_select/test.dat_zscore_',num,'.txt'),
                                    header = T, row.names = 1, stringsAsFactors = F))
tra.data <- tcga_model_data[tra.samples, ]
dim(tra.data)
writeMatrix(tra.data,'results/07.module/tcga.train.dat.txt')
test.data <- tcga_model_data[test.samples, ]
dim(test.data)
writeMatrix(test.data,'results/07.module/tcga.test.dat.txt')

########随机分组卡方检验
colnames(tcga_cli)
tcga.cli.chisq=data.frame(cohort=rep(c('Train','Test'),c(length(tra.samples),length(test.samples))),
                          tcga_cli[c(tra.samples,test.samples),c(12,3:8,11)])
writeMatrix(tcga.cli.chisq,'results/07.module/tcga.cli.chisq.txt',row = F)


#单因素
tra.cox <- t(apply(tra.data[,3:c(ncol(tra.data))],2,function(x){
  vl=as.numeric(x)
  tm=tra.data$OS.time
  ev=tra.data$OS
  #ev=ifelse(ev=='Alive',0,1)
  dat=data.frame(tm,ev,vl)[which(tm > 0 & !is.na(vl)),]
  return(coxFun(dat))
}))
colnames(tra.cox)=c('p.value','HR','Low 95%CI','High 95%CI')
length(which(tra.cox[,1]<0.01))
tra.cox <- na.omit(tra.cox)
writeMatrix(round(tra.cox,3),'results/07.module/tra.cox.txt')

filter_genes <- rownames(tra.cox[tra.cox[,1]<0.01, ])
filter_genes_dat=data.frame(Names=rownames(tra.cox[tra.cox[,1]<0.01,]),tra.cox[tra.cox[,1]<0.01,])
length(filter_genes)#30

tra.cox=as.data.frame(tra.cox)
tra.cox$coef=log(tra.cox[,2])
tra.cox$type=rep('None',nrow(tra.cox))
tra.cox$type[which(tra.cox$p.value<0.01 & tra.cox$coef>0)]='Risk'
tra.cox$type[which(tra.cox$p.value<0.01 & tra.cox$coef<0)]='Protective'


#lasso
library(glmnet)
set.seed(num)
fit1=glmnet(as.matrix(tra.data[,filter_genes])
            #,factor(samps)
            ,cbind(time=tra.data$OS.time,
                   status=tra.data$OS)
            ,family="cox"
            #,family="binomial"
            #,type.measure="deviance"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix(tra.data[,filter_genes])
                  #,factor(samps)
                  ,cbind(time=tra.data$OS.time,
                         status=tra.data$OS)
                  ,family="cox"
                  #,family="binomial"
                  #,type.measure="deviance"
                  ,nlambda=100
                  , alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
pdf('results/07.module/lasso.pdf',height = 10,width = 5)
mg_plot_lasso(fit1,
              cv.fit,
              # lambda = cv.fit$lambda.min,
              show_text=T,
              figLabels=c('A','B'))
dev.off()
# 基因的多因素
tcga_dat1 <- cbind(time=tra.data$OS.time,
                   status=tra.data$OS,
                   tra.data[,names(sig.coef)])

fmla <- as.formula(paste0("Surv(time, status) ~"
                          ,paste0(names(sig.coef),collapse = '+')))

length(names(sig.coef))
cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
cox1=step(cox)
pdf('results/07.module/Fig7a.pdf',height = 5,width = 8,onefile = F)
ggforest(cox1, data =tcga_dat1, 
         main = "Hazardratio", 
         fontsize =1.0, 
         noDigits = 2)#画森林图
dev.off()
lan <- coef(cox1)
round(lan, 3)
genes <- names(cox1$coefficients)
length(genes)#
paste0(round(lan, 3), '*', names(lan),collapse = '+')
tra.cox[genes,]

# 训练集
risk.tr <- as.numeric(lan%*%as.matrix(t(tra.data[,genes])))
risk.trz <- mosaic::zscore(risk.tr)

tcga.tra.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ ifelse(risk.tr>median(risk.tr),'H','L'),
                                data = tra.data),
                   data=tcga_cli[tra.samples,],
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='TCGA-train',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   palette = ggsci::pal_lancet()(8),
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   #legend = c(0.8,0.75), # 指定图例位置
                   legend.title = "",
                   legend.labs = c("High","Low"))
fig7a1=mg_merge_plot(tcga.tra.km$plot,tcga.tra.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig7a1
fig7a2=ggplotTimeROC(tra.data$OS.time,
             tra.data$OS,
             risk.trz,mks = c(1,2,3,4,5))

# 测试集
fmla.test <- as.formula(paste0("Surv(OS.time, OS) ~"
                               ,paste0(genes,collapse = '+')))
cox.test <- coxph(fmla.test, data =as.data.frame(test.data))
test_lan <- coef(cox.test)
risk.te=as.numeric(test_lan%*%as.matrix(t(test.data[,genes])))

risk.tez=mosaic::zscore(risk.te)
table(test.data$OS.time > 365 * 3)

tcga.te.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ ifelse(risk.te>median(risk.te),'H','L'),
                                    data = test.data),
                       data=tcga_cli[test.samples,],
                       conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                       title='TCGA-test',ggtheme=custom_theme(),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = ggsci::pal_lancet()(8),
                       legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                       #legend = c(0.8,0.75), # 指定图例位置
                       legend.title = "",
                       legend.labs = c("High","Low"))
fig7b1=mg_merge_plot(tcga.te.km$plot,tcga.te.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig7b1

fig7b2=ggplotTimeROC(test.data$OS.time,
                    test.data$OS,
                    risk.tez,mks = c(1,2,3,4,5))

# tcga
fmla.tcga <- as.formula(paste0("Surv(OS.time, OS) ~"
                               ,paste0(genes,collapse = '+')))
cox.tcga <- coxph(fmla.tcga, data =as.data.frame(tcga_model_data))
tcga_lan <- coef(cox.tcga)
risk.tcga=as.numeric(tcga_lan%*%as.matrix(t(tcga_model_data[tcga_cli$Samples,genes])))

risk.tcgaz=mosaic::zscore(risk.tcga)

tcga.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ ifelse(risk.tcga>median(risk.tcga),'H','L'),
                                   data = tcga_model_data),
                      data=tcga_cli,
                      conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                      title='TCGA-STAD',ggtheme=custom_theme(),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = ggsci::pal_lancet()(8),
                      legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                      #legend = c(0.8,0.75), # 指定图例位置
                      legend.title = "",
                      legend.labs = c("High","Low"))
fig7c1=mg_merge_plot(tcga.km$plot,tcga.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig7c1
fig7c2=ggplotTimeROC(tcga_model_data$OS.time,
                    tcga_model_data$OS,
                    risk.tcgaz,mks = c(1,2,3,4,5))


########### GSE66229 数据集##############
GSE66229_model_data <- t(GSE66229_exp[intersect(pre.genes,rownames(GSE66229_exp)),GSE66229_cli$Samples])
GSE66229_model_data=merge(data.frame(Samples=GSE66229_cli$Samples,OS.time=GSE66229_cli$OS.time,OS=GSE66229_cli$OS),
                          data.frame(Samples=rownames(GSE66229_model_data),GSE66229_model_data),
                          by='Samples')
rownames(GSE66229_model_data)=GSE66229_model_data$Samples
GSE66229_model_data=GSE66229_model_data[,-1]
GSE66229_model_data=crbind2DataFrame(GSE66229_model_data)
dim(GSE66229_model_data)

fmla.GSE66229 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                   ,paste0(intersect(colnames(GSE66229_model_data),genes),collapse = '+')))
cox.GSE66229 <- coxph(fmla.GSE66229, data =as.data.frame(GSE66229_model_data))
GSE66229_lan <- coef(cox.GSE66229)
risk.GSE66229=as.numeric(GSE66229_lan%*%as.matrix(t(GSE66229_model_data[,intersect(colnames(GSE66229_model_data),genes)])))

# risk.GSE66229=as.numeric(lan%*%as.matrix(t(GSE66229_model_data[,genes])))
risk.GSE66229z=mosaic::zscore(risk.GSE66229)

GSE66229.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ ifelse(risk.GSE66229>median(risk.GSE66229),'H','L'),
                                data = GSE66229_model_data),
                   data=GSE66229_cli,
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='GSE66229',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   palette = ggsci::pal_lancet()(8),
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   #legend = c(0.8,0.75), # 指定图例位置
                   legend.title = "",
                   legend.labs = c("High","Low"))
fig7d1=mg_merge_plot(GSE66229.km$plot,GSE66229.km$table,ncol = 1,nrow = 2,heights = c(3.5,1),align = 'v')
fig7d1

fig7d2=ggplotTimeROC(GSE66229_model_data$OS.time,
                    GSE66229_model_data$OS,
                    risk.GSE66229z,mks = c(1,2,3,4,5))



#####单因素+多因素+nomogram########
fig7b=mg_merge_plot(fig7a1,fig7a2,fig7b1,fig7b2,fig7c1,fig7c2,fig7d1,fig7d2,ncol = 4,nrow = 2,
                    labels = c('A','','B','','C','','D',''))
savePDF('results/07.module/Fig7b.pdf',fig7b,height = 10,width = 16)
tcga.risktype.cli=data.frame(tcga_cli,
                             Cluster=tcga.subtype[tcga_cli$Samples,]$Cluster,
                             Riskscore=risk.tcga,
                             Risktype=ifelse(risk.tcga>median(risk.tcga),'High','Low'))

head(tcga.risktype.cli)
writeMatrix(tcga.risktype.cli,outpath = 'results/08.risktype.cli/tcga.risktype.cli.txt',row = F)

tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)
table(tcga_cox_datas$T.stage)
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T1'|tcga_cox_datas$T.stage=='T2']<-'T1+T2'
tcga_cox_datas$T.stage[tcga_cox_datas$T.stage=='T3'|tcga_cox_datas$T.stage=='T4']<-'T3+T4'
table(tcga_cox_datas$N.stage)
tcga_cox_datas$N.stage[tcga_cox_datas$N.stage=='N1'|tcga_cox_datas$N.stage=='N2'|tcga_cox_datas$N.stage=='N3']<-'N1+N2+N3'
table(tcga_cox_datas$M.stage)
table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'
table(tcga_cox_datas$Grade)
tcga_cox_datas$Grade[which(tcga_cox_datas$Grade=='G1')]=NA
#Age
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~T.stage,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~N.stage,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

#M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~M.stage,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat

#Stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

#Grade
Grade_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Grade,
                               data=tcga_cox_datas))
Grade_sig_cox_dat <- data.frame(Names=rownames(Grade_sig_cox[[8]]),
                                HR = round(Grade_sig_cox[[7]][,2],3),
                                lower.95 = round(Grade_sig_cox[[8]][,3],3),
                                upper.95 = round(Grade_sig_cox[[8]][,4],3),
                                p.value=round(Grade_sig_cox[[7]][,5],3))
Grade_sig_cox_dat

#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     T.stage_sig_cox_dat,
                     N.stage_sig_cox_dat,
                     M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        "T.stage",
                        "N.stage",
                        "M.stage",
                        "Stage",
                        "Grade",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
writeMatrix(data.sig,'results/07.module/data.sig.txt',row = F)
pdf('results/07.module/Fig7c.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap = 8,lineheight = 10)
dev.off()

#多因素
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~ Age +Gender+Stage+Grade+Risktype, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Names=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Age",
                         "Gender",
                         "Stage",
                         "Grade",
                         "RiskScore")
data.muti$Names <- rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/07.module/Fig7d.pdf', width = 8, height = 5,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap = 8,lineheight = 10)
dev.off()

####################列线图

pdf('results/07.module/Fig7_Nomogram.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age,
                                Stage=tcga_cox_datas$Stage),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,
                     mks = c(1,3,5))
dev.off()
mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))



##############08.临床特征########
dir.create('results/08.risktype.cli')
# tcga.risktype.cli=data.frame(tcga_cli,
#                              Riskscore=tcga.module.risk$result$riskscore,
#                              Risktype=ifelse(tcga.module.risk$result$riskscore>median(tcga.module.risk$result$riskscore),'High','Low'))
# 
head(tcga.risktype.cli)
fig8a=mg_violin_1(data = tcga.risktype.cli[,c('Stage',"Riskscore")],
            melt = T,leg.title = 'Stage',ylab = 'Riskscore',
            group_col = ggsci::pal_lancet()(8),jitter = F)

fig8b=mg_violin_1(data = tcga.risktype.cli[,c('Grade',"Riskscore")],
            melt = T,leg.title = 'Grade',ylab = 'Riskscore',
            group_col = ggsci::pal_lancet()(8),jitter = F)

fig8c=mg_violin_1(data = tcga.risktype.cli[,c('Cluster',"Riskscore")],
            melt = T,leg.title = 'Cluster',ylab = 'Riskscore',
            group_col = ggsci::pal_lancet()(8),jitter = F,
            test_method = 'wilcox.test',cmp_test_method = 'wilcox.test')
fig8=mg_merge_plot(fig8a,fig8b,fig8c,labels = LETTERS[1:3],ncol = 3,nrow = 1,widths = c(1,1,0.75))
savePDF('results/08.risktype.cli/Fig8a.pdf',fig8,height = 4,width = 10)

chisq9=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$Age1))
chisq10=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$Gender))
chisq11=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$T.stage))
chisq12=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$N.stage))
chisq13=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$M.stage))
chisq14=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$Stage))
chisq15=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$Grade))
chisq16=chisq.test(table(tcga.risktype.cli$Risktype,tcga.risktype.cli$Status))

cli.chisq.test2=data.frame(Features=c('Age1','Gender','T.stage','N.stage','M.stage','Stage','Grade','Status'),
           pval=round(c(chisq9$p.value,chisq10$p.value,chisq11$p.value,chisq12$p.value,
                        chisq13$p.value,chisq14$p.value,chisq15$p.value,chisq15=chisq16$p.value),3))
cli.chisq.test2
writeMatrix(cli.chisq.test2,'results/08.risktype.cli/risktype.cli.chisq.test.txt',row = F)
# Features  pval
# 1     Age1 0.362
# 2   Gender 1.000
# 3  T.stage 0.020
# 4  N.stage 0.452
# 5  M.stage 0.363
# 6    Stage 0.130
# 7    Grade 0.224
# 8   Status 0.000

head(tcga.risktype.cli)
cli_anno1=tcga.risktype.cli[order(tcga.risktype.cli$Riskscore),c(14,15,13,11,12,3:8)]
cli_anno1$N.stage[which(cli_anno1$N.stage=='N1'|cli_anno1$N.stage=='N2'|cli_anno1$N.stage=='N3')]<-'N1+N2+N3'
table(cli_anno1$N.stage)


risktype.color=pal_lancet()(8)[1:2]
names(risktype.color)=c('High','Low')

# age.col=c("#BEBADA","#FB8072")
# names(age.col)=c('<=65','>65')
# 
# sex.col=c("#BEBADA","#FB8072")
# names(sex.col)=c('MALE','FEMALE')
# 
# tstage.col=brewer.pal(n = 4,"Set2")
# names(tstage.col)=c('T1','T2','T3','T4')
# 
# nstage.col=c("#8DD3C7","#FFFFB3")#brewer.pal(n = 4,"Set3")
# names(nstage.col)=c('N0','N1+N2+N3')
# 
# mstage.col=c("#8DA0CB","#E78AC3")
# names(mstage.col)=c('M0','M1')
# 
# stage.col=brewer.pal(n = 4,"Accent")
# names(stage.col)=c('I','II','III','IV')
# 
# grade.col=brewer.pal(n = 3,"Pastel2")
# names(grade.col)=c('G1','G2','G3')
# 
# status.col=c("#E41A1C","#377EB8")
# names(status.col)=c('Alive','Dead')

colnames(cli_anno1)
color_anno1=list(Risktype=risktype.color,Age1=age.col,Gender=sex.col,T.stage=tstage.col,N.stage=nstage.col,M.stage=mstage.col,
                Stage=stage.col,Grade=grade.col,Status=status.col,Riskscore=brewer.pal(7,"YlOrRd"),Cluster=cluster.color)
pdf('results/08.risktype.cli/Fig8b.pdf',height = 8,width = 10)
pheatmap(tcga_tpm_log_T[genes,rownames(cli_anno1)],
         scale = 'row',
         #border="white", # 设置边框为白色
         color =  colorRampPalette(c("navy", "white", "red"))(200),
         main="Heatmap", # 设置图形标题
         annotation_col = cli_anno1,
         #annotation_row = gene_anno,
         annotation_colors = color_anno1,
         #cutree_cols = 2, #列划为2块，
         #cutree_rows =2, # 行为2块
         #cellwidth = 6,cellheight = 5, # 设置热图方块
         display_numbers = F, # 热图上显示数值
         cluster_cols = F, # 去掉横向、纵向聚类
         cluster_rows = F,
         show_rownames = T, #去掉横、纵坐标id
         show_colnames = F,
         #legend_breaks=c(-5,0,5),
         fontsize_row = 12, # 分别设置横向和纵向字体大小
         fontsize_col = 16) 
dev.off()


#####09.分组间免疫相关性#############
dir.create('results/09.risktype.immu')
##########tcga.ciber基本无意义
mg_PlotMutiBoxplot(data = tcga.ciber[tcga.risktype.cli$Samples,1:22],
                   group =  tcga.risktype.cli$Risktype,
                   test_method = c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[2],
                   ylab = 'Score',group_cols = ggsci::pal_nejm()(10)[3:4],
                   legend.pos = 'top',
                   add = 'boxplot')

p.dat6=cbind(Risktype=tcga.risktype.cli$Risktype,tcga.immu.ssgsea[tcga.risktype.cli$Samples,])
p.dat6=crbind2DataFrame(p.dat6)
p.dat6=melt(p.dat6)
head(p.dat6)

p.immu1=
  p.dat6 %>% 
  ggplot(aes(x=variable, y=value,fill=Risktype)) +
  geom_boxplot()+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values =risktype.color)+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Risktype",title = 'ssGSEA') +
  #theme_light()+
  theme_classic()+
  theme(legend.position = 'top',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 将图表标题居中


p.dat7=cbind(Risktype=tcga.risktype.cli$Risktype,tcga.mcp[tcga.risktype.cli$Samples,])
p.dat7=crbind2DataFrame(p.dat7)
p.dat7=melt(p.dat7)
head(p.dat7)
p.immu2=
  p.dat7 %>% 
  ggplot(aes(x=variable, y=value,fill=Risktype)) +
  geom_boxplot()+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values =risktype.color)+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Risktype",title = 'MCP-Count') +
  #theme_light()+
  theme_classic()+
  theme(legend.position = 'top',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) # 将图表标题居中
p.immu2
figs2=mg_merge_plot(p.immu1,p.immu2,ncol=1,nrow=2,labels = c('A','B'),common.legend = T)
savePDF('results/09.risktype.immu/FigS2.pdf',figs2,height = 9,width = 10)

mg_PlotMutiBoxplot(data = tcga.est[tcga.risktype.cli$Samples,],
                   group = tcga.risktype.cli$Risktype,
                   test_method = c('t.test','wilcox.test','paired_t.test','paired_wilcox.test','anova','kruskal.test')[2],
                   ylab = 'Score',group_cols =risktype.color,
                   legend.pos = 'top',
                   add = 'boxplot')

tcga.est.use=tcga.est[tcga.risktype.cli$Samples,]
colnames(tcga.est.use)=paste0(colnames(tcga.est.use),'_EST')

tcga.immu.ssgsea.use=tcga.immu.ssgsea[tcga.risktype.cli$Samples,]
colnames(tcga.immu.ssgsea.use)=paste0(colnames(tcga.immu.ssgsea.use),'_ssGSEA')

tcga.mcp.use=tcga.mcp[tcga.risktype.cli$Samples,]
colnames(tcga.mcp.use)=paste0(colnames(tcga.mcp.use),'_MCP-COUNT')


immu.dat=cbind(tcga.est.use,tcga.immu.ssgsea.use,tcga.mcp.use)
immu.heatmap.dat=t(immu.dat)

immu_cell_type=data.frame(immu_cell=rownames(immu.heatmap.dat))
immu_cell_type$type=rep('none',nrow(immu_cell_type))
immu_cell_type$type[grep('EST',immu_cell_type$immu_cell)]<-'EST'
immu_cell_type$type[grep('ssGSEA',immu_cell_type$immu_cell)]<-'ssGSEA'
immu_cell_type$type[grep('MCP-COUNT',immu_cell_type$immu_cell)]<-'MCP-COUNT'

head(immu_cell_type)

###############相关性气泡图
immu.RS.dat=cbind(tcga.risktype.cli[,'Riskscore'],tcga.est.use,tcga.immu.ssgsea.use,tcga.mcp.use)
colnames(immu.RS.dat)[1]='Riskscore'
head(immu.RS.dat)
immu.cor.RS=Hmisc::rcorr(as.matrix(immu.RS.dat))
immu.cor.RS.res=data.frame(Names=names(immu.cor.RS$r['Riskscore',]),
                         cor=as.numeric(immu.cor.RS$r['Riskscore',]),
                         pval=as.numeric(immu.cor.RS$P['Riskscore',]))
immu.cor.RS.res$type=rep('none',nrow(immu.cor.RS.res))
immu.cor.RS.res$type[grep('EST',immu.cor.RS.res$Names)]<-'EST'
immu.cor.RS.res$type[grep('ssGSEA',immu.cor.RS.res$Names)]<-'ssGSEA'
immu.cor.RS.res$type[grep('MCP-COUNT',immu.cor.RS.res$Names)]<-'MCP-COUNT'
immu.cor.RS.res=immu.cor.RS.res[-1,]
immu.cor.RS.res$pval=ifelse(immu.cor.RS.res$pval<0.0001,'****',
                            ifelse(immu.cor.RS.res$pval<0.001,'***',
                                   ifelse(immu.cor.RS.res$pval<0.01,'**',
                                          ifelse(immu.cor.RS.res$pval<0.05,'*',''))))

immu.cor.RS.res$Names <- factor(immu.cor.RS.res$Names
                              ,levels=immu.cor.RS.res$Names[order(immu.cor.RS.res$type,decreasing = T)], ordered=TRUE)

head(immu.cor.RS.res)
writeMatrix(immu.cor.RS.res,'results/09.risktype.immu/immu.cor.RS.res.txt',row = F)

pdf('results/09.risktype.immu/Fig9a.pdf',height = 8,width = 12)
ggplot(immu.cor.RS.res,aes(x=cor,y=Names,color=type))+
  geom_point(aes(size=pval))+xlab('Correlation coefficient')+ylab('Immune Cell')+
  #scale_color_gradient(low = 'red',high = 'blue')+
  theme_bw()+theme(axis.text.y = element_text(size =10))
dev.off()

################TIDE+免疫评分
tcga.immu.dat.ri=cbind.data.frame(Risktype=tcga.risktype.cli[,'Risktype'],
                                  TIDE=tcga_tide_res[tcga.risktype.cli$Samples,'TIDE'],
                                  tcga.est[tcga.risktype.cli$Samples,])
head(tcga.immu.dat.ri)
tcga.immu.dat.ri=melt(tcga.immu.dat.ri)
head(tcga.immu.dat.ri)
colnames(tcga.immu.dat.ri)=c('Risktype','type','value')
head(tcga.immu.dat.ri)
p.tide2=tcga.immu.dat.ri %>% 
  ggplot(aes(x=Risktype, y=value,fill=Risktype)) +
  geom_boxplot()+ facet_wrap(~type,scales = 'free',nrow = 1,ncol = 4)+
  #scale_color_manual(values = sub.col) +  #箱线图颜色
  scale_fill_manual(values = risktype.color)+   #箱线图填充颜色
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  labs(x="", y = "Score", fill = "Risktype") +
  theme_light()+
  theme(legend.position = 'none',                 #图例位置
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 8, hjust = 0.5)) # 将图表标题居中
p.tide2
savePDF('results/09.risktype.immu/Fig9b.pdf',p.tide2,height = 4,width = 12)

save.image(file = '20221018_STAD_metabolic.RData')
