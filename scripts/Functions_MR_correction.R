library(data.table)

### Function to check if alleles are alligned or flipped

flipper=function(x){
  
  if(x[1]==x[4] & x[2]==x[3]){
    return(as.numeric(x[5])*-1)
  }else{
    return(NA)
  }
  
}


gwma.data.frame=function(bGWAS.obj=res.bg           # bgwas output object
                         ,correct.se=T              # bgwas infletes the se shuould it be corrected to be the original se (this option makes no difference if no.se==T)
                         ,chi.df=1                  # Degrees of freedom for test statistic of corrected pvalues
                         ,z_filter=1                # pvalue threshold below which the expected estimate is considered equal to 0, default means no filter is applied
                         ,no.se=T                   # Should the correction take into account the uncertainty in the expected effect estimation, 
                                                    #default TRUE which means that the point estimate is used without taking into account the associated SE.
                         ,results.data.frame=res.t  # data frame with the original results, it should coincide with the one used as input for the bGWAS step
                         ,rs.lab="rsid"             # Label for the SNP identifier column
                         ){
  
  res=results.data.frame
  
  
  names(res)[grep(rs.lab,names(res))]="rs"
  prior.res=data.frame(extract_results_bGWAS(bGWAS.obj,"all"))
  res=merge(x=prior.res,y=res[,c("rs","beta1","se")],by="rs",all.x=T)
  #### alt=a1 ref=a0
  ### prepare for flipping
  
  to.flip=which(sign(res$beta1)!=sign(res$observed_Z) )
  res$beta1[to.flip]=res$beta1[to.flip]*-1
  
  if(correct.se==T) res$prior_std_error=sqrt(res$prior_std_error^2-1)
  
  res$expected_p=pchisq(res$prior_estimate^2,df=1,lower=F)
  
  if(no.se==FALSE){
    res$z_expected=res$prior_estimate/res$prior_std_error
    res$expected_p=pchisq(res$z_expected^2,df=1,lower=F)
  }
  res$prior_estimate[which(res$expected_p>z_filter)]=0
  if(no.se==T){
    res$z_diff=res$observed_Z-res$prior_estimate
    
  }else{
    res$z_diff=(res$observed_Z-res$prior_estimate)/sqrt(1+(res$prior_std_error^2))
  }
  res$corr.beta=res$z_diff*res$se
  res$corr.p=pchisq(res$z_diff^2,df=chi.df,lower=F)
  res$p=pchisq(res$observed_Z^2,df=1,lower=F)
  res$corr2raw_ratio=res$corr.beta/res$beta1
  return(res)
  
  
}


#----
# Define miami function
#----

miami <- function (x, y=NULL, chr = "chr", bp = "pos", p = "p", snp = "rsid", col1 = c("#e34a33", "#fdbb84"), col2 = c("#2b8cbe","#a6bddb"), chrlabs = NULL, suggestiveline = -log10(1e-05),
                   genomewideline = -log10(1e-08), ymin = -15, ymax = 15, x.name = "x", y.name = "y", highlight = NULL, highlight_col="red", label = NULL, logp = TRUE, ymid = 0, ...)
{
  CHR = BP = P = index = NULL
  ymid <- abs(ymid)
  ymin = ymin + ymid
  ymax = ymax - ymid
  
  #----
  # Check variables
  #----
  
  column_check <- function(column, df) {
    if (!(column %in% names(df)))
      stop(paste0("Column ", column, " not found in ",deparse(substitute(df)),"!"))
    if (!(column == snp) & !is.numeric(df[[column]]))
      stop(paste(column, "column in",deparse(substitute(df)),"should be numeric."))
  }
  
  for (column in c(bp,p,chr,snp)) column_check(column, x)
  if(!is.null(y)) for (column in c(bp,p,chr,snp)) column_check(column, y)
  
  #----
  # Create data frame
  #----
  print("##LOAD D1")
  
  d1 = data.table(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) d1$SNP <- x[[snp]]
  if (!is.null(highlight)) d1$highlight <- x[[highlight]]
  d1 <- subset(d1, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d1 <- d1[order(d1$CHR, d1$BP), ]
  if (logp) {
    d1$logp <- -log10(d1$P) - ymid
    d1 <- d1[d1$logp > 0,]
  } else {
    d1$logp <- d1$P - ymid
    d1 <- d1[d1$logp > 0,]
  }
  
  
  if(!is.null(y)) {
    print("##LOAD D2")
    d2 = data.table(CHR = y[[chr]], BP = y[[bp]], P = y[[p]])
    if (!is.null(y[[snp]])) d2$SNP <- y[[snp]]
    if (!is.null(highlight)) d2$highlight <- y[[highlight]]
    d2 <- subset(d2, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d2 <- d2[order(d2$CHR, d2$BP), ]
    if (logp) {
      d2$logp <- log10(d2$P) + ymid
      d2 <- d2[d2$logp < 0,]
    } else {
      d2$logp <- -d2$P + ymid
      d2 <- d2[d2$logp < 0,]
    }
    
    d <- data.frame(rbind(d1,d2))
    
  } else {
    d <- data.frame(d1)
  }
  
  #d <- d[abs(d$logp) > ymid,]
  
  #----
  # Define plotting parameters
  #----
  
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  } else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index ==
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP +
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index ==
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  
  print("##CALL PLOT")
  
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", yaxt = "n", cex.axis = 1.5, cex.lab = 1.5,
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(ymin-1,ymax+1), xlab = xlabel, ylab = expression(-log[10](italic(p))),mar=5*c(5.1,6.1,4.1,2.1))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  yticks <- seq(-5*ceiling(abs(ymin/5)), 5*ceiling(abs(ymax/5)), by = 5)
  yticks <- sort(c(yticks+sign(yticks)*(5-ymid),5-ymid,-5+ymid))
  if(is.null(y)) yticks <- yticks[yticks>=0]
  ylabels <- abs(yticks)+ymid
  
  if ((abs(suggestiveline) - ymid) > 0) {
    abline(h = suggestiveline - ymid, col = "blue")
    abline(h = -suggestiveline + ymid, col = "blue")
  }
  if ((abs(genomewideline) - ymid) > 0) {
    abline(h = genomewideline - ymid, col = "red")
    abline(h = -genomewideline + ymid, col = "red")
  }
  abline(h = 0, col = "black")
  if(!is.null(y)) legend(0, y = ymax, c(x.name,y.name), fill = c(col1[1],col2[1]), bty="n", cex = 1.5)
  
  
  if (nchr == 1) {
    axis(1, cex.axis = 1.5, cex.lab = 1.5, las=1, ...)
    axis(2, at = yticks, labels = ylabels, las=1, cex.axis=1.5, cex.lab = 1.5, ...)
  } else {
    axis(1, at = ticks, labels = labs, las=1, cex.axis = 1.5, cex.lab = 1.5, ...)
    axis(2, at = yticks, labels = ylabels, las=1, cex.axis = 1.5, cex.lab = 1.5, ...)
  }
  col1 = rep(col1, max(d$CHR))
  col2 = rep(col2, max(d$CHR))
  
  d1 <- d[d$logp>0,]
  d2 <- d[d$logp<0,]
  
  if (nchr == 1) {
    with(d1[d1$logp<ymax,], points(pos, logp, pch = 20, col = col1[1], ...))
    with(d1[d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), pch = 24, col = col1[1], ...))
    with(d2[d$logp>ymin,], points(pos, logp, pch = 20, col = col2[1], ...))
    with(d2[d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), pch = 25, col = col2[1], ...))
  } else {
    print("##PLOT X")
    icol = 1
    for (i in unique(d1$index)) {
      
      if(!is.null(highlight)) {
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 0, ], points(pos, logp, col = col1[icol], pch = 20, ...))  
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex= 1.5, lwd= 3, ...))
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 2, ], points(pos, logp, col = highlight_col, bg=highlight_col, pch = 24, cex= 1.5, lwd= 2, ...))
        with(d1[d1$index == unique(d1$index)[i] & d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), col = highlight_col, bg=highlight_col, pch = 24, cex=1.5, lwd = 2, ...))
      } else {
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax, ], points(pos, logp, col = col1[icol], pch = 20,cex=0.8, ...))  
        with(d1[d1$index == unique(d1$index)[i] & d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), pch = 24, col = col1[icol], cex=1, lwd = 1, ...))
      }
      
      
      icol = icol + 1
    }
    
    print("##PLOT Y")
    if(!is.null(y)) {
      
      icol = 1
      for (i in unique(d2$index)) {
        
        if(!is.null(highlight)) {
          with(d2[d2$index == unique(d2$index)[i] & d2$logp>ymin & d2$highlight == 0, ], points(pos, logp, col = col2[icol], pch = 20, ...))
          with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymax & d2$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex = 1.5, lwd=3, ...))
          with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymax & d2$highlight == 2, ], points(pos, logp, col = highlight_col, bg=highlight_col, pch = 25, cex = 1.5, lwd=2, ...))
          with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), col = highlight_col, bg=highlight_col, pch = 25, cex=1.5, lwd=2, ...))
        } else {
          with(d2[d2$index == unique(d2$index)[i] & d2$logp>ymin, ], points(pos, logp, col = col2[icol], pch = 20,cex=0.8, ...))
          with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), pch = 25, col = col2[icol], cex=1, lwd=1, ...))
        }
        
        
        icol = icol + 1
      }
    }
  }
  
  
  
  #----
  # Plot names
  #----
  
  if (!is.null(label)) {
    required_cols <- c("LABEL","CHR","POS","P")
    if (is.null(colnames(label))) warning(paste("label must be a data.frame with col.names:",paste0(required_cols,collapse=" ")))
    label <- data.frame(label)
    
    if (all(required_cols %in% colnames(label))) {  
      label$pos <- label$POS
      for (i in unique(label$CHR)[!unique(label$CHR)==1]) {
        lastbase <- tail(d[d$index==(i-1),"pos"],n=1)
        label[label$CHR==i,"pos"] <- label[label$CHR==i,"pos"]+lastbase
      }
      
      #with(label, points(pos, P, col = "#974d26", pch = 4, cex=2, ...))
      
      label$P <- ifelse(label$P>0, label$P-ymid+2*strheight(label$LABEL), label$P+ymid-2*strheight(label$LABEL))
      label$P <- ifelse(label$P>ymax,ymax-strheight(label$LABEL)/2,label$P)
      label$P <- ifelse(label$P<ymin,ymin+strheight(label$LABEL)/2,label$P)
      label <- label[order(abs(label$P),label$pos),]
      
      moved <- label[0,]
      if (nrow(label)>1) {
        
        n <- nrow(label)
        repeat {
          for(i in 1:(nrow(label)-1)) {
            for(j in (i+1):nrow(label)) {
              if(abs(label[i,"pos"]-label[j,"pos"]) < 1.2*(strwidth(label[i,"LABEL"]) + strwidth(label[j,"LABEL"])) & 
                 abs(label[i,"P"]-label[j,"P"]) < 2.5*strheight(label[i,"LABEL"]) &
                 label[i,"P"]<ymax & label[i,"P"] > ymin ) {
                if(abs(label[i,"P"]) < abs(label[j,"P"])) { 
                  moved <- rbind(moved, label[i,])
                  label[i,"P"] <- label[j,"P"] + sign(label[i,"P"])*2.5*strheight(label[i,"LABEL"])
                  if(!is.null(highlight)) {
                    if( any(abs(d[d$highlight,"logp"]-label[i,"P"]) < 2*strheight(label[i,"LABEL"]) & 
                            abs(d[d$highlight,"pos"] - label[i, "pos"]) < 2*strwidth(label[i,"LABEL"])) ) {
                      label[i,"P"] <- label[i,"P"] + sign(label[i,"P"])*2*strheight(label[i,"LABEL"])
                    }
                  }
                  moved <- rbind(moved, label[i,])
                } else {
                  moved <- rbind(moved, label[j,])
                  label[j,"P"] <- label[i,"P"] + sign(label[j,"P"])*2.5*strheight(label[i,"LABEL"])
                  if( any(abs(d[d$highlight,"logp"]-label[j,"P"]) < 2*strheight(label[j,"LABEL"]) & 
                          abs(d[d$highlight,"pos"] - label[j, "pos"]) < 2*strwidth(label[j,"LABEL"])) ) {
                    label[j,"P"] <- label[j,"P"] + sign(label[j,"P"])*2*strheight(label[j,"LABEL"])
                  }
                  moved <- rbind(moved, label[j,])
                }
              }
            }
          }
          if(nrow(moved)==n) break
          n <- nrow(moved)
        }
        
      }
      label$P <- ifelse(label$P>ymax,ymax-strheight(label$LABEL)/6,label$P)
      label$P <- ifelse(label$P<ymin,ymin+strheight(label$LABEL)/6,label$P)
      
      moved <- data.table(moved)
      moved <- moved[order(pos,P),]
      m1 <- moved[P>0,.(x=pos,y1=min(P),y2=max(P)),by=c("LABEL")]
      m2 <- moved[P<0,.(x=pos,y1=min(P),y2=max(P)),by=c("LABEL")]
      
      if(nrow(m1)>0) {
        for(i in 1:nrow(m1)) lines(x=m1[i,x] + c(0,0), y=c(m1[i,y1],m1[i,y2]) - strheight(m1[i,LABEL]))
      }
      if(nrow(m2)>0) {
        for(i in 1:nrow(m2)) lines(x=m2[i,x] + c(0,0), y=c(m2[i,y1],m2[i,y2]) + strheight(m2[i,LABEL]))
      }
      
      shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                             theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
        xy <- xy.coords(x,y)
        xo <- r*strwidth('A')
        yo <- r*strheight('A')
        # draw background text with small shift in x and y in background colour
        for (i in theta) {
          text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
        }
        # draw actual text in exact xy position in foreground colour
        text(xy$x, xy$y, labels, col=col, ... )
      }
      
      with(label, shadowtext(x=pos,y=P,labels=LABEL,cex=1.5,col="black",bg="white",font=4,r=0.2))
      
    } else {
      warning(paste("label data.frame is missing columns:",paste0(required_cols[!required_cols%in%colnames(label)],collapse=" ")))
    }
  }
  
  par(xpd = FALSE)
}

gwma.miami=function(gwma.obj=corrected.gwas, output.file="Cheese_miamiplot.bmp",main="Cheese consumption"){

ps=gwma.obj[,c("chrm","pos","rs","p")]
cps=gwma.obj[,c("chrm","pos","rs","corr.p")]
names(cps)[4]="p"
bmp(paste0(output.file),height=960, width=2400)
miami (x=ps, y=cps, chr = "chrm", bp = "pos", p = "p", snp = "rs"
       , col1 = c("#e34a33", "#fdbb84"), col2 = c("#2b8cbe","#a6bddb")
       , chrlabs = NULL, suggestiveline = -log10(1e-05),
       genomewideline = -log10(1e-08), ymin = -30
       , ymax = 30, x.name = "Raw results", y.name = "Corrected results", highlight = NULL
       , highlight_col="red", label = NULL, logp = TRUE, ymid = 0,main=main)

dev.off()

}








