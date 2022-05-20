library(mvabund)
x <- as.matrix(read.csv('./genus80_top30.csv',row.names=1))
x <- as.matrix(read.csv('../feature_tables/52N_genus_t80_selbal.csv',row.names=1))
meta <- as.data.frame(read.csv('../meta/meta52_current.csv',row.names=1))

x <- as.mvabund(x)
plot(x)
meanvar.plot(x)

mod1 <- manyglm(x ~ Zosyn, data=meta, cor.type="I",)
plot(mod1)
anov <- anova.manyglm(mod1,p.uni="unadjusted",cor.type='I')
anov <- summary.manyglm(mod2,p.uni="adjusted")
a <- colnames(x)
summary.manyglm(mod1)
names(anov)
anov$uni.p
anov$table
write.csv(anov$table,'piptazo.csv')
write.csv(anov$uni.p,'piptazouni.csv')
