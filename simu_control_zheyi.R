
library(stats)
library(brainwaver)
library(igraph)
library(tnet)
library(geonames)
library(rjson)
library(fields)

library(fossil)
library(plotrix)
library(animation)
library(geonames)
library(rjson)
library(fields)
library(matrixcalc)
library(stats)
library(sna)
library(igraph)

setwd('C:/Sync/temp/Code/')

#############################
###       load data       ###
#############################
cities <- read.csv("0india_node.csv", header=TRUE)
link<-read.csv("0india_edge.csv", header=TRUE)
filled<-read.csv("0filled_in.csv", header=TRUE)
link_org<-read.csv("0india_edge_org3.csv", header=TRUE)

### from edgelist(dataframe) to adjacency matrix ###
adj_org <-matrix(0,nrow=67,ncol=67)
colnames(adj_org) <- cities[,2]
rownames(adj_org) <- cities[,2]

count_cnt = dim(link_org)[1]

for(cnt in 1:count_cnt)
{
	rowo = link_org[cnt,4]
	colo = link_org[cnt,5]
	
	adj_org[rowo,colo]=adj_org[rowo,colo]+link_org[cnt,9]
	adj_org[colo,rowo]=adj_org[colo,rowo]+link_org[cnt,9]
}
### from adjacency matrix to graph ###
g_org <-graph.adjacency(adj_org,mode="undirected",weighted=TRUE)

####################################################
### calculate network index:                     ###
###    modularity, clustering(transitivity),     ###
###    strength(degree), clustering? path length ###
###    efficiency???
####################################################
fc1<-fastgreedy.community(g_org)
mod_org<-modularity(g_org,membership(fc1),E(g_org)$weight)

clu_org1<-transitivity(g_org,type="globalundirected",weights=NULL)

deg_org<-graph.strength(g_org)

net_st<-cbind(get.edgelist(g_org,names=FALSE),E(g_org)$weight)
# #For undirected networks, the symmetrise function must be run
#if(!is.directed(g_st))
net_st<-as.tnet(net_st,type="weighted one-mode tnet")
net_st<-symmetrise_w(net_st)
# #Ensure that it conforms to the tnet standard

# # g2<-tnet_igraph(net,type="weighted one-mode tnet")
##result.clu<-clustering_w(net, measure=c("bi", "am", "gm", "ma", "mi"))
clu_org2<-clustering_w(net_st)


path_org<-average.path.length(g_org, directed=FALSE, unconnected=TRUE)

### 不熟悉global.efficiency, 注意这里用的是0/1数据
b_adj_org <- 1*as.matrix(adj_org!=0)
adjc<-dim(adj_org)[1]
eff<-global.efficiency(as.matrix(b_adj_org), matrix(1,adjc,adjc))
eff_org<-mean(eff$nodal.eff)






#######################################
### modeling & parameter estimation ###
### method: brute force             ###
#######################################

totLink = sum(adj_org)/2

totPop = sum(cities[,7])

pop = log(cities[,7]+1)
nat = cities[,8]
x = cities[,3]
y = cities[,4]

# factors
popF = pop%*%t(pop) 
distF = rdist.earth(cbind(x,y),cbind(x,y),miles=FALSE)
distF = log(distF+1)
### 矩阵乘积是啥...
### 为什么都取了log呢？这是在干啥

### optinum value:
# # eps=1 # disturbance
# alpha = 1.5 #distance
# beta = 0.3 #population
# theta = 2 #border same country
# gamma = 1  #degree

nocity = dim(cities)[1]
nolink = dim(link)[1]

brutalsize = 5;
count2 = 1
mode = 2

result = matrix(0,nrow = brutalsize^4, ncol = 17)
colnames(result)<-c('alpha','beta','theta','gamma','qap','qap-sig','mod','clu_u','clu_w','ave.path','eff','pE','pC','pM','pP','pDeg','hamming')
ptm=proc.time()


	alpha = 4
	beta = 2
	theta = 1
	eps = 1
	gamma = 1
	
	ctrF = matrix(1,nrow=nocity,ncol=nocity)
	
	for (ti in 1:nocity)
	{
		for (tj in 1:nocity)
		{
			if(nat[ti]!=nat[tj])
			{
				
			ctrF[ti,tj]=theta
			ctrF[tj,ti]=theta
			} 
			
			
		}
	}	
	
	
tempmod=matrix(0,100,1)
tempeff=matrix(0,100,1)
tempclu=matrix(0,100,1)
temppath=matrix(0,100,1)

tempdeg =matrix(0,67,1)

adj_rese <-matrix(0,nrow=nocity,ncol=nocity)
	colnames(adj_rese) <- cities[,2]
	rownames(adj_rese) <- cities[,2]
	
	indF = matrix(1,67,67)

########################
### simulation start ###	
for (rep in 1:100)
{
	distF1 = distF^alpha
	popF1 = popF^beta
	ctrF1 = ctrF^theta
	

	adj_simu <-matrix(0,nrow=nocity,ncol=nocity)
	colnames(adj_simu) <- cities[,2]
	rownames(adj_simu) <- cities[,2]
	


	
for (j in 1:totLink)
{	
	
	
	if(mode==1){
	degr <- as.matrix(rowSums(adj_simu))+eps
	tempD <- degr%*%t(degr)
	degrF1<-tempD^gamma
	probp <-indF*degrF1*popF1/distF1/ctrF1
	diag(probp)<-0
	} else if(mode==2){
	
	tempD<-(adj_simu)%*%(adj_simu)+eps
	degrF1<-tempD^gamma
	probp <-degrF1*popF1/distF1/ctrF1
	diag(probp)<-0	
	} else if (mode ==7){
	degr <- as.matrix(rowSums(adj_simu))	
	
	degrF1<-matrix(0,nrow=nocity,ncol=nocity)
	
	for (tpi in 1:nocity)
	{
		for(tpj in 1:nocity)
		{
			degrF1[tpi,tpj] = degr[tpi]+degr[tpj]+eps
		}
	}
	degrF1<-degrF1^gamma
	probp <-degrF1*popF1/distF1/ctrF1
	diag(probp)<-0
		
	}
	
		
	probp2<-vec(probp)
	sumProb<-sum(probp2)
	probp2<-probp2/sumProb
	
	a<-sample.int(n=nocity*nocity,size=1,prob=probp2)
		
	city1 = (a-1)%%nocity+1
	city2 = (a-1)%/%nocity+1
	
	
	adj_simu[city1,city2]<-adj_simu[city1,city2]+1
	adj_simu[city2,city1]<-adj_simu[city2,city1]+1
	
	# if(adj_simu[city1,city2]==5)
	# {
	# indF[city1,city2]=0
		# indF[city2,city1]=0
	# }
		
}


count_cntf = dim(filled)[1]

for(cnt in 1:count_cntf)
{
	rowf = filled[cnt,4]
	colf = filled[cnt,5]
	
	if(adj_simu[rowf,colf]==0){
	
	adj_simu[rowf,colf]=adj_simu[rowf,colf]+filled[cnt,9]
	adj_simu[colf,rowf]=adj_simu[colf,rowf]+filled[cnt,9]
	}
}


g_simu <-graph.adjacency(adj_simu,mode="undirected",weighted=TRUE)

# if(is.connected(g_simu))
# {
	# g_simu = g_simu
# }else{
	# mcg<-clusters(g_simu)
	# comp1 = which(mcg$membership == which.max(mcg$csize))
	# indtemp = which.max(mcg$csize)

		# for(h in 1:mcg$no)
	# {
		# if(h!=indtemp)
		# {comp2 = which(mcg$membership == h)
		
		# indt<-which(probp == max(probp[comp1,comp2]), arr.ind=TRUE)
		# ra_c1 = indt[1,1]
		# ra_c2 = indt[1,2]
		
		# adj_simu[ra_c1,ra_c2]=1
		# adj_simu[ra_c2,ra_c1]=1

		# }
	# }
	# g_simu <-graph.adjacency(adj_simu,mode="undirected",weighted=TRUE)
	
	# }
 
 adj_rese = adj_rese+adj_simu
 
qapcorr<-array(dim = c(2,67,67))
qapcorr[1,,] = adj_org
qapcorr[2,,] = adj_simu
qapresult<-qaptest(qapcorr,gcor,g1=1,g2=2)
corrcoef=qapresult$testval
corrsig=qapresult$pgreq




result[count2,1]=alpha
result[count2,2]=beta
result[count2,3]=theta
result[count2,4]=gamma


result[count2,5]=result[count2,5]+corrcoef
result[count2,6]=result[count2,6]+corrsig

fc1<-fastgreedy.community(g_simu)
modi<-modularity(g_simu,membership(fc1),E(g_simu)$weight)
result[count2,7]<-result[count2,7]+modi

clui1<-transitivity(g_simu,type="globalundirected",weights=NULL)

result[count2,8]<-result[count2,8]+clui1

net_st<-cbind(get.edgelist(g_simu,names=FALSE),E(g_simu)$weight)
# #For undirected networks, the symmetrise function must be run
if(!is.directed(g_simu))
{net_st<-as.tnet(net_st,type="weighted one-mode tnet")
	}
net_st<-symmetrise_w(net_st)
# #Ensure that it conforms to the tnet standard

# # g2<-tnet_igraph(net,type="weighted one-mode tnet")
##result.clu<-clustering_w(net, measure=c("bi", "am", "gm", "ma", "mi"))
clui2<-clustering_w(net_st)
result[count2,9]<-result[count2,9]+clui2

pathi2<-average.path.length(g_simu, directed=FALSE, unconnected=TRUE)
result[count2,10]<-result[count2,10]+pathi2

b_adj_simu <- 1*as.matrix(adj_simu!=0)
adjc<-dim(adj_simu)[1]
eff<-global.efficiency(as.matrix(b_adj_simu), matrix(1,adjc,adjc))
effi<-mean(eff$nodal.eff)
result[count2,11]<-result[count2,11]+effi


result[count2,17]<-result[count2,17]+sum(abs(adj_simu-adj_org))

tempmod[rep] = modi
tempeff[rep] = effi
tempclu[rep] = clui2
temppath[rep] = pathi2


tempdeg<-tempdeg+ as.matrix(graph.strength(g_simu))


}





#######################
### hypothesis test ###
#######################
	
pe<-t.test(tempeff,mu=eff_org)
pE<-pe$p.value
	
pc<-t.test(tempclu,mu=clu_org2)
pC<-pc$p.value
	
pm<-t.test(tempmod,mu=mod_org)
pM<-pm$p.value
	
pp<-t.test(temppath,mu=path_org)
pP<-pp$p.value	

deg_stat = round(as.matrix(tempdeg/100),0)
	
# pdeg<-chisq.test(deg_stat,deg_org)
deg_s = density(as.matrix(deg_stat))
org_s = density(as.matrix(deg_org))
pdeg<-ks.test(org_s$y,deg_s$y)
pDeg<-pdeg$p.value


result[count2,5]=result[count2,5]/100
result[count2,6]=result[count2,6]/100
result[count2,7]=result[count2,7]/100
result[count2,8]=result[count2,8]/100
result[count2,9]=result[count2,9]/100
result[count2,10]=result[count2,10]/100
result[count2,11]=result[count2,11]/100
result[count2,17]=result[count2,17]/100

result[count2,12] = pE
result[count2,13] = pC
result[count2,14] = pM
result[count2,15] = pP
result[count2,16] = pDeg

 adj_rese  =  adj_rese/100

count2 = count2+1

dt = proc.time() - ptm


######################
###  save results ####
######################
write.csv(result, file = "/data/india/result_control3_4_2_1_1.csv")

write.csv( adj_rese,file="/data/india/accum_matrix_4_2_1_1.csv")


adj_rese2 <-matrix(0,nrow=nocity,ncol=nocity)
	colnames(adj_rese2) <- cities[,2]
	rownames(adj_rese2) <- cities[,2]
	adj_rese2 = adj_rese
	
	diff2= adj_rese-adj_rese2
	diffv2<-as.matrix(diff2[lower.tri(diff2)])
	
	write.csv( diffv2,file="/data/india/diffv2.csv")
	
	