

rm(list=ls())


library(R.matlab)
library(RColorBrewer)

# summary of q
ncase <- 2
dtab <- matrix(NA,ncase,50)
tmp <- list(ncase)
for(case in 1:ncase){
  tmp[[case]] <- readMat(paste('paras',case,'.mat',sep=''))
  tmpd <- round(table(tmp[[case]]$ps)/sum(table(tmp[[case]]$ps)),4)
  dtab[case,as.integer(names(tmpd))+1] <- tmpd
}


ds <- which(apply(dtab,2,function(x)return(any(!is.na(x)))))
dtab <- dtab[,ds]
dtab1 <- t(rbind(ds-1, dtab))
findmax <- function(x) return(which(x==max(x,na.rm=T)))
dtab1 <- rbind(dtab1, dtab1[apply(dtab1,2,findmax),1])
dtab1 <- as.data.frame(dtab1)
names(dtab1) <- c('number of selected variables')
dtab1[nrow(dtab1),1] <- 'posterior mode:'
# write.csv(dtab1,file='dtab.csv',na='.',row.names=F)
qmods <- as.integer(dtab1[nrow(dtab1),2:ncol(dtab1)])

# figure version
fig1 <- as.data.frame((dtab)); names(fig1) <- ds-1; 
row.names(fig1) <- c('GDD hypothesis 1', 'GDD hypothesis 2')



dirs <- './' #./tex/'
#dirs <- 'C:\\Users\\Zhen\\Google Drive\\Yang\\chapter1\\tex\\'


# postscript(file=paste(dirs,'histq.eps',sep=''),pointsize=11,width=5.5,height=2.5,horizontal=F,paper='special')

cols0 = brewer.pal(11,'RdYlGn')[c(6,9)] #6
# palette(rainbow(K))
par(mfrow=c(1,1),pch=18,xaxt='s',yaxt='s',mar=c(1.5,2,1,0)+.1,cex.main=.8, font.main=1, cex.axis=.8, mgp=c(1.1,.2,0), tck=-0.01)
barplot(as.matrix(fig1),beside=TRUE, col=cols0, legend = rownames(fig1), args.legend=list(cex=.8))
#abline(h=0, col='red',lty=2)

dev.off()




## MCMC traceplot of the qs
#lattice.options(default.theme = canonical.theme(color = FALSE))
#my.padding <- list(strip = 1, top.padding = 0, main.key.padding = 0, key.axis.padding = 0, axis.ylab.padding = 0, ylab.key.padding  = 0, right.padding = 0, left.padding = 0, key.sub.padding = 0, bottom.padding = 0, axis.right = 0, key.right = 0)
#
#
#
#nam_case <- c('Distance')
#ncase <- length(nam_case)
#nch <- 5; nsample <- 1e3; tot <- nch*nsample
#K <- length(nam_case)
#
#fig2 <- data.frame()
#for(case in 1:K) fig2 <- rbind(fig2, cbind(as.vector(tmp[[case]]$allps), rep(1:nsample,nch)))
#names(fig2) <- c('y','x')
#fig2$grp <- as.factor(rep(rep(paste('chain =',1:nch), each=nsample), K))
#fig2$case <- factor(rep(nam_case, each=tot), levels=nam_case)
#
#cols <- brewer.pal(nch,'Paired') #rainbow(nch)
#myfig <- xyplot(y~x|case, groups=grp,data=fig2, type='l',xlab='',ylab=list('Number of selected variables', cex=.8),layout=c(1,1),col=cols[1:nch], lty=1)
#myfig <- update(myfig, par.settings = list(layout.heights = my.padding, layout.widths = my.padding), 
# strip = strip.custom(par.strip.text = list(cex=.8, lines=2, font=1)), scales = list(alternating = T, cex = .7, abbreviate = F, 
# y=list(cex=.7,rot = c(0,0), tck=c(0.4,0.4)), 
# x=list(cex=.7,rot = c(0,0), tck=c(0.4,0.4))) 
#)
#
##postscript(file=paste('./tex/trace2.eps',sep=''),pointsize=11,width=5.5,height=2.5,horizontal=F,paper='special')
## print(myfig)
##dev.off()

# rename for the figures
load(file='Xnam')
Xnam[Xnam=='tag_length'] <- 'length'
Xnam[Xnam=='years_diff'] <- 'years_lag'
# Xnam[Xnam=='d5_tag'] <- 'GDD'
Xnam[Xnam=='Dip'] <- 'Diporeia'
Xnam[Xnam=='sexFemale'] <- 'sex: Female'
Xnam[Xnam=='GDD1'] <- 'GDD_Diff'
Xnam[Xnam=='GDD2_M9'] <- 'Sep*GDDlake'
Xnam[Xnam=='GDD2_M10'] <- 'Oct*GDDlake'
Xnam[Xnam=='GDD2_M11'] <- 'Nov*GDDlake'
Xnam[Xnam=='GDD2_M12'] <- 'Dec*GDDlake'
Xnam <- gsub('rec_Y','rec_Y: ',Xnam)
Xnam <- gsub('tag_Y','tag_Y: ',Xnam)
Xnam <- gsub('rec_M','rec_M: ',Xnam)
Xnam <- gsub('tag_site','tag_site: ',Xnam)
Xnam <- paste(Xnam, " ") # add some gap from axis

# variable wise summary

nam_ind <- c(rep(1,2),rep(2,6))
var_all <- list(Xnam)
# xadjs <- c(.15,.55,.5,.37,.64,.27,.31,.31)
xadjs <- c(1,.55,.5,.37,.64,.27,.31,.31)



# par(mfrow=c(2,2))

#========================================  loop starts
for(case in 1:ncase){
  
  #case <- 1
  # nams0 <- var_all[[nam_ind[case]]]
  
  # be careful
  if(case == 1) nams0 <- Xnam[1:33]
  if(case == 2) nams0 <- Xnam[-33]
  
  
  
  
  
  a1 <- tmp[[case]]$mat.var
  a0 <- cbind(nams0[a1[,1]], round(a1[,-1],4))
  a0 <- as.data.frame(a0)
  for(j in 2:ncol(a0)) a0[,j] <- as.numeric(as.character(a0[,j]))
  #names(a0) <- c('Variable', 'Posterior Mean', '2.5%Q', '97.5%Q', 'Selectivity')
  # write.csv(a0[1:nmax,], file='varmat.csv',row.names=F)
  names(a0) <- c('Variable', 'Mean', 'Lower', 'Upper', 'Selectivity')
  
  fig3 <- as.data.frame(nams0)
  fig3 <- merge(x=fig3, y=a0, by.x='nams0',by.y='Variable')
  #h0 <- 1e5; if(qmods[case]>0) h0 <- a0$Selectivity[qmods[case]]
  #fig3$sig <- as.integer(fig3$Selectivity >= h0)
  fig3$sig <- (fig3$Lower*fig3$Upper > 0)
  ra <- range(as.vector(fig3[,2:4]), na.rm=T)+.05*c(-1,1)
  
  # be careful here
  nams <- fig3$nams0
  #nams[6:16] <- paste('rec_M',c(11,12,1:10),sep='')
  #nams[15:22] <- paste('rec_Y',c(2003:2005,2007:2011),sep='')
  #all(sort(levels(fig3$nams0)) == sort(nams))    # must be TRUE
  fig3$nams0 <- factor(fig3$nams0, levels=Xnam)
  # fig3$nams0[order(fig3$nams0)]
  fig3 <- fig3[order(fig3$nams0, decreasing =T),]
  
  
  
  
  
  # postscript(file=paste(dirs,'vs',case,'.eps',sep=''),pointsize=10,width=6.4,height=7.5,horizontal=F,paper='special')
  # pdf(file=paste('./tex/kF1.pdf',sep=''),pointsize=9,width=6,height=4.2,paper='special')
  
  jpeg(paste(dirs,'vs',case,'.jpg',sep=''), quality=100, height=7.5,width=6.4, units='in',pointsize=11, res=300) 
  
  collab <- c('black','red') #'gray60'
  cols <- collab[fig3$sig+1]
  # change mar() for best margin
  par(yaxt='n',xaxt='s', mar=c(1, 9.3, .3, 0)+.3, mgp=c(1.2,0.1,0), tck=-0.01, cex.axis=.8, cex.main=.8)
  qs <- nrow(fig3); vec <- 1:nrow(fig3); #vec[vec>qs] <- NA
  plotCI(y=vec,x=fig3$Mean,li=fig3$Lower,ui=fig3$Upper, col='black', barcol=cols, gap=0, sfrac=0.004, lty=1, main='',xlab="",ylab="", pch=16, cex=.5, err='x')   #nam_case[case]
  abline(h=1:qs, col='gray90', lty=2)
  abline(v=0, col='blue')
  plotCI(y=1:qs,x=fig3$Mean,li=fig3$Lower,ui=fig3$Upper, col='black', barcol=cols, gap=0, sfrac=0.004, lty=1, main="",xlab="",ylab="", pch=16, cex=.8, add=T, lwd=2,err='x')
  #text(y=1:qs,x=rep(ra[1]-xadjs[case],qs), labels = fig3$nams0, srt = 0, adj = c(0,0),xpd = TRUE, cex=.85, col=collab[fig3$sig+1], font=1)#-1.35 for nvar=1
  par(yaxt='s', las=1, mgp=c(1.4,0.3,0),tck=-0.005)
  xlabs <- fig3$nams0 #gsub(patter='tag_site',x=fig3$nams0,replacement='tag_site: ')
  axis(side=2, at=1:qs,labels = xlabs,font=1, cex.axis=.9)
  text(y=1:qs+.2,x=fig3$Mean-.04, labels = fig3$Selectivity, srt = 0, adj = c(0,0),xpd = TRUE, cex=.75, col=collab[fig3$sig+1], font=1)
  
  dev.off()
  
  
  
  
  
  # model wise summary
  nvar <- case
  nmax <- 12
  modtab <- tmp[[nvar]]$tab
  bmat <- round(tmp[[nvar]]$bmat,4); lb <- round(tmp[[nvar]]$lb,4); ub <- round(tmp[[nvar]]$ub,4)
  bmat0 <- round(tmp[[nvar]]$bmat0,4);
  # bmat[1] is the beta for the best model
  p <- length(nams)
  modmat <- matrix(NA, nrow(modtab), p)
  mat_m <- mat_l <- mat_u <- matrix(NA, nrow(modtab), p)
  for(i in 1:p){
    inds <- which(modtab[,i]==1)
    modmat[inds,i] <- paste(bmat[inds,i], ' (', lb[inds,i], ', ', ub[inds,i], ')', sep='')
    mat_m[inds,i] <- bmat[inds,i]
    mat_l[inds,i] <- lb[inds,i]
    mat_u[inds,i] <- ub[inds,i]
  }
  modmat <- as.data.frame(modmat)
  names(modmat) <- nams0
  modmat$Selectivity <- round(modtab[,ncol(modtab)],4)
  del.inds <- integer()
  for(j in 1:ncol(modmat)) if(all(is.na(modmat[1:nmax,j]))) del.inds <- c(del.inds, j)
  modmat <- modmat[,-del.inds]
  sav.mod <- modmat[1:nmax,]
  # write.csv(sav.mod, file='modmat.csv', na='.', row.names=F)
  
  # summary of q for top models
  myqs <- apply(sav.mod, 1, function(x) sum(!is.na(x)))[1:12] -1 #exclude the last col: selectivity
  mys <- round(modmat$Selectivity[1:nmax], 4)
  
  nams4 <- names(sav.mod)[-ncol(sav.mod)]
  mat_m <- mat_m[1:nmax, -del.inds] # nmodel by nvariable
  mat_l <- mat_l[1:nmax, -del.inds]
  mat_u <- mat_u[1:nmax, -del.inds]
  
  
  q <- length(nams4)
  nm <- nrow(mat_m)
  fig3 <- as.data.frame(cbind(as.numeric(mat_m),as.numeric(mat_l),as.numeric(mat_u)))
  names(fig3) <- c('means','lbs','ubs')
  fig3$mods <- rep(1:nrow(mat_m), ncol(mat_m)) # model indicator
  fig3$vars <- rep(1:q, each=nrow(mat_m))
  inds <- order(aggregate(fig3$lbs, list(fig3$vars), min, na.rm=T)$x)
  fig3$vars <- factor(fig3$vars, levels=inds)
  fig3 <- fig3[order(fig3$vars), ]
  fig3$vars <- as.integer(fig3$vars)
  x0 <- seq(-1,1,len=nm)/2.1
  fig3$x <- fig3$vars - x0[fig3$mods]
  
  
  
  
  # postscript(file=paste(dirs,'ms',case,'.eps',sep=''),pointsize=9,width=6.4,height=8.4,horizontal=F,paper='special')
  
  jpeg(paste(dirs,'ms',case,'.jpg',sep=''), quality=100, height=8.4,width=6.4, units='in',pointsize=11, res=300) 
  
  cols1 <- c(brewer.pal(9,'Set1'),brewer.pal(12,'Paired'),brewer.pal(8,'Set2'))[1:q] # indicate the variables
  cols1[6] <- '#01665e' 
  # cols1 <- cols1[inds]
  cols1[9] <- '#3f007d'; cols1[10] <- 'black' # replace light colors with darker one.
  cols1 <- rev(cols1)
  #plot(1:length(cols1),col=cols1,cex=3,pch=16)
  #plot(1:length(cols1), col=cols1)
  pchs <- 1:25; pchs <- pchs[-c(1,3,9:16)]  # indicate the models  # rep(21, nm); #
  pchs[10] <- 15 
  pchs <- pchs[1:nm]  #rep(NA, nm)#
  ltys <- 1
  par(yaxt='n',xaxt='s', mar=c(1,10,0,0)+.5, mgp=c(1.2,0.4,0), tck=-0.01, cex.axis=1)
  plotCI(y=fig3$x,x=fig3$means,li=fig3$lbs,ui=fig3$ubs, pch=pchs[fig3$mods], barcol=cols1[fig3$vars], gap=0, sfrac=0.002, lty=ltys,lwd=1, main='',xlab="",ylab="", err='x',col=cols1[fig3$vars], cex=1,type='p',pt.bg='red')
  abline(v=0,lty=1,col='blue') 
  abline(h=1:q,lty=2,col='gray80') 
  plotCI(y=fig3$x,x=fig3$means,li=fig3$lbs,ui=fig3$ubs, pch=pchs[fig3$mods], barcol=cols1[fig3$vars], gap=0, sfrac=0.002, lty=ltys,lwd=1, main="",xlab="",ylab="", err='x',col=cols1[fig3$vars], cex=1, add=T,type='p',pt.bg='white')
  pchs2 <- c(1:nm)
  pos <- (fig3$lbs-.03)*(fig3$mods%%2) + (fig3$ubs+.03)*(!fig3$mods%%2)
  text(y=fig3$x,x=pos, labels=pchs2[fig3$mods], cex=.8)
  par(yaxt='s', las=1, mgp=c(1.2,0.24,0))
  xlabs <- nams4[inds]
  # xlabs <- gsub(patter='tag_site',x=nams4[inds],replacement='tag_site: ')
  axis(side=2, at=1:q,labels = xlabs,font=1, cex.axis=1)
  legend('bottomright', legend=paste('Model ',1:nm,' (',myqs,', ',mys,')', sep=''), inset=.02, col='black', bg='white', cex=1, lty=ltys,pch=pchs,pt.bg='white') # 
  
  dev.off()
  
  
  
  # # get formula for the top model
  # tmpa <- fig3[fig3$mods == 1, ]; row.names(tmpa) <- gsub(x=nams4[inds],pattern='_',replacement='\\_')
  # myeqn <- paste('{\\tt logdist} &=& ', bmat0[1])
  # for(j in 1:nrow(tmpa)){
  #   if(!is.na(tmpa$means[j])){
  #     tmpchar <- '+'; if(tmpa$means[j]<0) tmpchar <- '-'
  #     myeqn <- paste(myeqn, tmpchar, '{\\tt',row.names(tmpa)[j],'}\\times', abs(round(tmpa$means[j],3)))
  #    }
  # }
  # 
  # cat(myeqn, '\n')
  
  # if(case==1) tab1 <- tmpa  
  # if(case==2) tab2 <- tmpa  
  
}
#========================================  loop ends




# # be careful
# tab1$xnam <- row.names(tab1)
# tab2$xnam <- row.names(tab2)
# tab1 <- na.omit(tab1)
# tab2 <- na.omit(tab2)
# 
# row.names(tab1) <- NULL
# row.names(tab2) <- NULL
# 
# mat <- cbind(tab1[,-c(4:6)], tab2[,-c(4:6)])
# mat <- mat[,c(4,1:3, 8,5:7)]
# 
# # create the table for the manuscript
#   round1 <- function(x,d=3){
#     y <- x
#     #if(is.character(x)) if(any(grepl('_',x))) x <- gsub('_','\\\\_',x)
#     if(is.numeric(x)){
#     y <- as.character(round(x,d))
#     for(i in 1:length(y)){
#      k <- as.integer(gregexpr('.',y[i],fixed=T))
#      if(k<0) {y[i] <- paste(y[i],'.',sep=''); for(j in 1:d) y[i] <- paste(y[i],'0',sep='')}
#      if(k>=0){ b <- nchar(y[i])-as.integer(gregexpr('.',y[i],fixed=T))
#         if(b<d) for(j in 1:(d-b))y[i] <- paste(y[i],'0',sep='')
#      }
#     }
#    }
#     return(y)
#   }
# 
# mat <- apply(mat,2,round1)
# 
# for(i in 1:nrow(mat)){
#   #cat(paste(mat[i,],'&'),'\n')
#   for(j in 1:ncol(mat)){
#     x <- mat[i,j]
#     if(j %in% c(1,5)) {
#       if(any(grepl('_',x))) cat(gsub('_','\\\\_',x))
#       else cat(x)
#      }
#     else cat('$',x,'$', sep='')
#     if(j<ncol(mat)) cat('&')
#   }
#   cat('\\\\\n')
# }














