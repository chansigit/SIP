library( ggplot2 )
library(mgcv)

args<-commandArgs(T)
dataInput<-args[1]
outFileDir<-args[2]

data <- read.table( dataInput,header=TRUE,row.names=1 )

#data <- read.table( dataInput,header=TRUE,row.names=1 )

message( "Basic information:" )
print( dim(data) )
data.matrix <- as.matrix( data )
print(dim(data.matrix))


################################################################################
idx.matrix <- data.matrix > 0
num.genes.detected.per.cell <- apply( idx.matrix,MARGIN=2,FUN=sum )
#print(length(num.genes.detected.per.cell))
num.genes.detected.median <- median( num.genes.detected.per.cell )
print(num.genes.detected.median)
setEPS()
postscript(paste(outFileDir,sep="/","genes.detected.per.cell.eps"),width=10,height=8)
data.plot <- data.frame(
  count = num.genes.detected.per.cell
)
p <- ggplot( data.plot,aes(x=count) )

p +geom_histogram(binwidth=100,fill='lightblue1',color='gray',aes(y=..count..))+labs( x="Genes detected/Cell",y="Counts" )+theme_bw()+
   annotate( "text",label=paste0("Meidan: ",num.genes.detected.median),size=8,colour='black',x=2100,y=110 )+
   geom_density( aes(y=100*..count..) )+
   theme(plot.title=element_text(size=20),
         axis.title.y=element_text(size = 16, vjust=+0.2),
         axis.title.x=element_text(size = 16, vjust=-0.2),
         axis.text.y=element_text(size = 14),
         axis.text.x=element_text(size = 14),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank())
dev.off()

################################################################################
## filter genes and cells
cell.thr <- 5
cell.less.than.thr.genes <- which( num.genes.detected.per.cell < cell.thr )
if (length(cell.less.than.thr.genes) != 0){
    data.matrix.filtered <- data.matrix[,-cell.less.than.thr.genes ]
}else{
    data.matrix.filtered <- data.matrix
}

idx.matrix.cell = data.matrix.filtered > 0
num.genes.appearance <- apply( idx.matrix.cell,MARGIN=1,FUN=sum )
genes.no.appearance <- which( num.genes.appearance == 0 )
if (length(genes.no.appearance) != 0){
    data.matrix.filtered <- data.matrix.filtered[-genes.no.appearance,]
}

print(dim(data.matrix.filtered))

################################################################################
## Expression.matrix
cell.counts <- apply( data.matrix.filtered,MARGIN=2,FUN=sum )
normalization.data <- t( t(data.matrix.filtered) / cell.counts )
tpm.data <- normalization.data * 10000
log.tpm.data <- log2( tpm.data + 1 )

################################################################################
## feature extraction
mu <- apply( log.tpm.data,MARGIN=1,FUN=mean )
sigma <- apply( log.tpm.data,MARGIN=1,FUN=sd )

log.mu <- log10(mu)
log.cv.square <- 2 * log10( ( (sigma) / mu )  )

## fit a GAM model
log.cv2.mu <- data.frame(
  cv2 = log.cv.square,
  mu = log.mu
)

fit.gam.cv2.mu <- gam( cv2~s(mu),data=log.cv2.mu )
print( summary(fit.gam.cv2.mu) )
fit.gam <- fit.gam.cv2.mu$fitted.values

num.of.cells <- dim(data.matrix.filtered)[2]
chi2.stat <- num.of.cells * 10^(log.cv.square - fit.gam)
p.chi2.stat <- pchisq( chi2.stat,df=num.of.cells-1 )
pvalue.chi2 <- 1 - p.chi2.stat

sig.thr <- 0.05
sig.genes._ <- pvalue.chi2 < sig.thr
sig.genes.idx <- which( pvalue.chi2 < sig.thr )

log.tpm.with.sig.genes <- log.tpm.data[sig.genes.idx,]
print(dim(log.tpm.with.sig.genes))
#fit.lowess.cv2.mu <- lowess( x = log.mu,y = log.cv.square )
#idx.lowess <- match( log.mu,fit.lowess.cv2.mu$x )
#print(head(idx.lowess))

setEPS()
postscript( paste(outFileDir,sep='/', "genes.cv.mu.eps"),width=8,height=6 )

data.plot <- data.frame(
  x = log.mu,
  y = log.cv.square,
  y_gam = fit.gam,
  sig_genes = factor(sig.genes._)
)
p <- ggplot( data.plot,aes(x=x,y=y,colour=sig_genes) )
p+geom_point(size=0.1)+theme_bw()+
  geom_line( aes(x,y_gam),colour='red' )+
  annotate( 'text',label=("Fitted with 'gam' from 'mgcv' package"),size=4,x=-1.5,y=-1 )+
  labs(x="log2(mu)",y="log2(cv^2)")+
  theme(plot.title=element_text(size=20),
        axis.title.y=element_text(size = 16, vjust=+0.2),
        axis.title.x=element_text(size = 16, vjust=-0.2),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.off()


expr <- t(log.tpm.with.sig.genes)
expr.pca<-prcomp(expr, retx = TRUE, center = TRUE, scale = FALSE)
expr.pc<-data.frame(expr.pca$x[,1:10])
kmeans<-kmeans(expr.pc, 2)
expr.visual<-data.frame(expr.pc[,1:2])
expr.visual$cluster<-as.factor(kmeans$cluster)

setEPS()
postscript(paste(outFileDir,sep='/','expr.pca.kmeans.eps'), width=8, height=6)
p <- ggplot( expr.visual,aes(x=PC1,y=PC2,colour=cluster) )
p+geom_point(size=0.5)+theme_bw()+
  labs(x="PC1",y="PC2")+
  theme(plot.title=element_text(size=20),
        axis.title.y=element_text(size = 16, vjust=+0.2),
        axis.title.x=element_text(size = 16, vjust=-0.2),
        axis.text.y=element_text(size = 14),
        axis.text.x=element_text(size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dev.off()

write.table(log.tpm.with.sig.genes, paste(outFileDir,sep='/', "log.tpm.with.sig.genes.txt"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep='\t')
