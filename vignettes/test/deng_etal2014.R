

################   Deng et al 2014  data analysis  ####################################

library(devtools)

read.data1 = function() {
  x = tempfile()
  download.file('https://cdn.rawgit.com/kkdey/singleCellRNASeqMouseDeng2014/master/data/Deng2014MouseEsc.rda', destfile=x, quiet=TRUE)
  z = get(load((x)))
  return(z)
}

Deng2014MouseESC <- read.data1()

deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

voom_deng_counts <- limma::voom(deng.counts)$E

dat <- t(voom_deng_counts)

pr <- prcomp(dat);

par(mfrow=c(1,1))
plot(pr$x[,1], pr$x[,2])

cor_sample <- cov2cor(cov(voom_deng_counts))
nsamples <- 10;

library(vegan)
screeplot(pr, bstick=TRUE, npcs=20)
bstick_var <- bstick(pr)
countPC <- numPC_generate(pr, npcs=20)

fields::image.plot(cor_sample)

nsamples <- 10

ash_cor_sample <- ash_cor(cor_sample, nsamples)$ash_cor_PD;
fields::image.plot(as.matrix(ash_cor_sample))

ash_cor_sample2 <- ash_cor2(cor_sample, nsamples)$ash_cor_PD;
fields::image.plot(as.matrix(ash_cor_sample2))

ash_cor_sample3 <- ash_cor3(cor_sample, nsamples)$ash_cor_PD;
fields::image.plot(as.matrix(ash_cor_sample3))

library(corpcor)
strimmer_cor <- cor.shrink(t(dat))



diags <- eigen(cor_sample)$values
ash_cov_sample3 <- diag(sqrt(diags))%*%ash_cor_sample3%*%diag(sqrt(diags))
ash_cov_sample2 <- diag(sqrt(diags))%*%ash_cor_sample2%*%diag(sqrt(diags))
ash_cov_sample <- diag(sqrt(diags))%*%ash_cor_sample%*%diag(sqrt(diags))
cov_sample <- diag(sqrt(diags))%*%cor_sample%*%diag(sqrt(diags))
strimmer_cov_sample <- diag(sqrt(diags))%*%strimmer_cor%*%diag(sqrt(diags))

library(glasso)
glasso_cov_sample <- glasso::glasso(cov_sample, rho = 0.7)

eigen_ash_sample3 <- eigen(ash_cov_sample3, only.values = TRUE)
eigen_ash_sample2 <- eigen(ash_cov_sample2, only.values = TRUE)
eigen_ash_sample <- eigen(ash_cov_sample, only.values = TRUE)
eigen_sample <- eigen(cov_sample, only.values = TRUE)
eigen_strimmer_sample <- eigen(strimmer_cov_sample, only.values = TRUE)
eigen_glasso_sample <- eigen(glasso_cov_sample$w, only.values = TRUE)

library(ggplot2)

dim <- dim(cor_sample)[1]
dim <- 25
eigendata <- data.frame(
  eigenorder = 1:dim,
  sample_cov = sort(log(as.numeric(eigen_sample$values)+1),  decreasing=TRUE)[1:dim],
  ash_cov = sort(log(as.numeric(eigen_ash_sample$values)+1),  decreasing=TRUE)[1:dim],
  ash_cov2 = sort(log(as.numeric(eigen_ash_sample2$values)+1),  decreasing=TRUE)[1:dim],
  ash_cov3 = sort(log(as.numeric(eigen_ash_sample3$values)+1),  decreasing=TRUE)[1:dim],
  strimmer_cov = sort(log(as.numeric(eigen_strimmer_sample$values)+1),  decreasing=TRUE)[1:dim],
  glasso_cov = sort(log(as.numeric(eigen_glasso_sample$values)+1),  decreasing=TRUE)[1:dim]
)

colnames(eigendata) <- c("eigenorder",
                         "sample_cov",
                         "ash_cov",
                         "ash_cov2",
                         "ash_cov3",
                         "strimmer_cov",
                         "glasso_cov")

library(ggplot2)
ggplot(eigendata, aes(eigenorder)) +
  geom_line(aes(y = sample_cov, colour = "sample cov")) +
  geom_line(aes(y = ash_cov, colour = "ash cov")) +
  geom_line(aes(y = ash_cov2, colour = "ash cov2"))+
  geom_line(aes(y = ash_cov3, colour = "ash cov 3"))+
  geom_line(aes(y = strimmer_cov, colour = "strimmer cov"))+
  geom_line(aes(y = glasso_cov, colour = "glasso cov"))+
  xlab("Order of eigenvalues (sorted)")+
  ylab("log(Eigenvalues + 1)")+
  scale_colour_manual(values=c("blue", "red", "green", "orange", "pink", "brown"))+
  ggtitle(paste0("Eigenvalues distribution n/p=", round(dim(cor_sample)[1]/nsamples, 4),"\n for different shrinkage methods"))+
  theme(plot.title = element_text(lineheight=.8, face="bold"))


#################  image plots for the correlation matrices  ##########################

library(fields)
set.seed(1)
par(mfrow=c(2,2))
cols = gray.colors(100)
image.plot(cov2cor(pop_cov$Sigma), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample2), col=cols, nlevel=100)
image.plot(as.matrix(ash_cor_sample3), col=cols, nlevel=100)

###############   PCA  plot  #######################################

eigen_ash <- eigen(ash_cov_sample)
plot(pr$x[,1], pr$x[,2])

######################   eigenvalues vs prcomp  ############################

pr <- prcomp(dat)
color=c("red","blue","cornflowerblue","black","cyan","darkblue",
        "brown4","burlywood","darkgoldenrod1","darkgray","deepskyblue","darkkhaki",
        "firebrick","darkorchid","hotpink","green","magenta","yellow", "azure1","azure4");

plot(pr$x[,1], pr$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))

mean_dat <- (rep(1,dim(dat)[1]))%*% t(colMeans(dat))
dat_adjusted <- dat - mean_dat
pr1 <- eigen(dat_adjusted %*% t(dat_adjusted))
eigs <- eigen(cov_sample)$values

library(expm)
projected_dat <- sqrtm(ash_cov_sample)%*% sqrtm(solve(cov_sample)) %*% dat
pr_projected <- prcomp(projected_dat)
plot(pr_projected$x[,1], pr_projected$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))

projected_dat <- sqrtm(ash_cov_sample2)%*% sqrtm(solve(cov_sample)) %*% dat
pr_projected <- prcomp(projected_dat)
plot(pr_projected$x[,1], pr_projected$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))

projected_dat <- sqrtm(ash_cov_sample3)%*% sqrtm(solve(cov_sample)) %*% dat
pr_projected <- prcomp(projected_dat)
plot(pr_projected$x[,1], pr_projected$x[,2], col=color[factor(deng.meta_data$cell_type)],
     pch=20, cex=1, ylim=c(-200,200))

###############   Minimal Spanning Tree  #################################

library(igraph)
rownames(projected_dat) <- rownames(dat)
dist_tree <- as.matrix(stats::dist(projected_dat, method = "euclidean"))
system.time(mst2 <- ape::mst(dist_tree));

indices_vertex <- unique(factor(deng.meta_data$cell_type))
mst2.graph <- igraph::graph.adjacency(as.matrix(mst2));
plot(mst2.graph,
     edge.arrow.mode=0,
     vertex.size = 5,
     vertex.shape="square",
     vertex.label=NA,
     vertex.color=color[factor(deng.meta_data$cell_type)],
     layout=igraph::layout_as_tree(mst2.graph,11)
     #    layout= layout.kamada.kawai
     # layout= igraph::layout_with_kk
)
legend("topleft", legend=indices_vertex,  fill=color[factor(indices_vertex)])


dist_tree2 <- as.matrix(stats::dist(t(deng.counts), method = "euclidean"))
system.time(mst22 <- ape::mst(dist_tree2));

indices_vertex <- unique(factor(deng.meta_data$cell_type))
mst22.graph <- igraph::graph.adjacency(as.matrix(mst22));
plot(mst22.graph,
     edge.arrow.mode=0,
     vertex.size = 5,
     vertex.shape="square",
     vertex.label=NA,
     vertex.color=color[factor(deng.meta_data$cell_type)],
     layout=igraph::layout_as_tree(mst22.graph,11)
     #    layout= layout.kamada.kawai
     # layout= igraph::layout_with_kk
)
legend("topleft", legend=indices_vertex,  fill=color[factor(indices_vertex)])

###############   Checking for grouping    ##############################

library(igraph)
gg <- get.edgelist(mst22.graph)
gg_indices <- t(sapply(1:dim(gg)[1], function(x) return(c(which(rownames(dat) == gg[x,1]), which(rownames(dat) == gg[x,2])))))
gg_type <- t(apply(gg_indices, 1, function(x) return(deng.meta_data$cell_type[x])))

wrong_matches <- sum(sapply(1:dim(gg_type)[1], function(x)
                          {
                              if (gg_type[x,1] != gg_type[x,2]){
                                return(1)
                              }else{
                                return(0)
                              }}))/dim(gg_type)[1]


###################   the inverse correlation matrix  ############################

invcor_shrink_strimmer <- invcor.shrink(t(dat));

invcor_shrink_strimmer[which(abs(invcor_shrink_strimmer) < 0.3)]=0;
row.names(invcor_shrink_strimmer) <- deng.meta_data$cell_type;

grids.net <- network::network(as.matrix(abs(invcor_shrink_strimmer)),   vertex.attr=NULL, vertex.attrnames=row.names(glasso_out$wi.cor.shrunk),
                              directed=F, hyper=FALSE, loops=FALSE, multiple=FALSE, bipartite = F)

adjacency.matrix <- as.sociomatrix(grids.net)
all(grids.net[,]==adjacency.matrix)
edgelist.matrix <- as.matrix(grids.net,matrix.type="edgelist")

grids.net %v% "vertex.names" <- row.names(invcor_shrink_strimmer)
#get.vertex.attribute(grids.net,"vertex.names")

signed.weight.vec <- array(0,dim(edgelist.matrix)[1])
for(m in 1:dim(edgelist.matrix)[1])
{
  signed.weight.vec[m] <- invcor_shrink_strimmer[edgelist.matrix[m,1], edgelist.matrix[m,2]];
}
signed.vec <- array(0,dim(edgelist.matrix)[1])
signed.vec[which(signed.weight.vec>0)]="Pos"
signed.vec[which(signed.weight.vec<0)]="Neg"

color.vec <- plyr::mapvalues(signed.vec, from = c("Pos", "Neg"), to = c("blue", "green"))

grids.net %e% "color.vec" <- color.vec

list.edge.attributes(grids.net)

plot(grids.net,
     label=grids.net %v% "vertex.names",
     thresh=0.1,
     displaylabels=TRUE,
     label.cex=rep(0.4, dim(invcor_shrink_strimmer)[1]),
     edge.col=grids.net %e% "color.vec",
     edge.label.cex=1,
     vertex.cex= rep(0.5, dim(invcor_shrink_strimmer)[1]),
     vertex.col= rep(2, dim(invcor_shrink_strimmer)[1]),
     vertex.lwd= rep(1, dim(invcor_shrink_strimmer)[1]),
    # boxed.labels=TRUE,
     mode="kamadakawai")

##################   inverse ash correlation matrix   ###########################

invcor_shrink_ash <- solve(ash_cor_sample);

invcor_shrink_ash[which(abs(invcor_shrink_ash) < 0.3)]=0;
row.names(invcor_shrink_ash) <- deng.meta_data$cell_type;

grids.net <- network::network(as.matrix(abs(invcor_shrink_ash)),   vertex.attr=NULL, vertex.attrnames=row.names(glasso_out$wi.cor.shrunk),
                              directed=F, hyper=FALSE, loops=FALSE, multiple=FALSE, bipartite = F)

adjacency.matrix <- as.sociomatrix(grids.net)
all(grids.net[,]==adjacency.matrix)
edgelist.matrix <- as.matrix(grids.net,matrix.type="edgelist")

grids.net %v% "vertex.names" <- row.names(invcor_shrink_ash)
#get.vertex.attribute(grids.net,"vertex.names")

signed.weight.vec <- array(0,dim(edgelist.matrix)[1])
for(m in 1:dim(edgelist.matrix)[1])
{
  signed.weight.vec[m] <- invcor_shrink_ash[edgelist.matrix[m,1], edgelist.matrix[m,2]];
}
signed.vec <- array(0,dim(edgelist.matrix)[1])
signed.vec[which(signed.weight.vec>0)]="Pos"
signed.vec[which(signed.weight.vec<0)]="Neg"

color.vec <- plyr::mapvalues(signed.vec, from = c("Pos", "Neg"), to = c("blue", "green"))

grids.net %e% "color.vec" <- color.vec

list.edge.attributes(grids.net)

plot(grids.net,
     label=grids.net %v% "vertex.names",
     thresh=0.1,
     displaylabels=TRUE,
     label.cex=rep(0.4, dim(invcor_shrink_strimmer)[1]),
     edge.col=grids.net %e% "color.vec",
     edge.label.cex=1,
     vertex.cex= rep(0.5, dim(invcor_shrink_strimmer)[1]),
     vertex.col= rep(2, dim(invcor_shrink_strimmer)[1]),
     vertex.lwd= rep(1, dim(invcor_shrink_strimmer)[1]),
     # boxed.labels=TRUE,
     mode="kamadakawai")

#######################   GLASSO  estimates   ####################################

system.time(glasso_out <-glasso::glasso(cov_sample, rho=0.7))

invcor_shrink_glasso <-glasso_out$wi;

#invcor_shrink_ash[which(abs(invcor_shrink_ash) < 0.3)]=0;
row.names(invcor_shrink_glasso) <- deng.meta_data$cell_type;

grids.net <- network::network(as.matrix(abs(invcor_shrink_glasso)),   vertex.attr=NULL, vertex.attrnames=row.names(glasso_out$wi.cor.shrunk),
                              directed=F, hyper=FALSE, loops=FALSE, multiple=FALSE, bipartite = F)

adjacency.matrix <- as.sociomatrix(grids.net)
all(grids.net[,]==adjacency.matrix)
edgelist.matrix <- as.matrix(grids.net,matrix.type="edgelist")

grids.net %v% "vertex.names" <- row.names(invcor_shrink_glasso)
#get.vertex.attribute(grids.net,"vertex.names")

signed.weight.vec <- array(0,dim(edgelist.matrix)[1])
for(m in 1:dim(edgelist.matrix)[1])
{
  signed.weight.vec[m] <- invcor_shrink_glasso[edgelist.matrix[m,1], edgelist.matrix[m,2]];
}
signed.vec <- array(0,dim(edgelist.matrix)[1])
signed.vec[which(signed.weight.vec>0)]="Pos"
signed.vec[which(signed.weight.vec<0)]="Neg"

color.vec <- plyr::mapvalues(signed.vec, from = c("Pos", "Neg"), to = c("blue", "green"))

grids.net %e% "color.vec" <- color.vec

plot(grids.net,
     label=grids.net %v% "vertex.names",
     thresh=0.1,
     displaylabels=TRUE,
     label.cex=rep(0.4, dim(invcor_shrink_glasso)[1]),
     edge.col=grids.net %e% "color.vec",
     edge.label.cex=1,
     vertex.cex= rep(0.5, dim(invcor_shrink_glasso)[1]),
     vertex.col= rep(2, dim(invcor_shrink_glasso)[1]),
     vertex.lwd= rep(1, dim(invcor_shrink_glasso)[1]),
     # boxed.labels=TRUE,
     mode="circle")
