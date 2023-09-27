
library(stringr)
library(paleobioDB)
library(mclust)
library(msir)
library(stats4)
library(iNEXT)
library(doParallel)

# bootstrap re-sampling paleobioDB for cycad + angio abundance/diversity

angio_data <- pbdb_occurrences(limit="all", base_name="Angiospermae", show=c("phylo", "ident", "time"))
#write.csv(angio_data,"angio_paleoDB_data.csv") # write out as csv to avoid repeat downloads (saves time)
#angio_data=read.csv("angio_paleoDB_data.csv")
mid_age=(as.numeric(angio_data $eag)+as.numeric(angio_data $lag))/2
angio_data =cbind(angio_data,mid_age)

cycad_data <- pbdb_occurrences(limit="all", base_name="Cycadales", show=c("phylo", "ident", "time"))
#write.csv(cycad_data,"cycad_paleoDB_data.csv") # write out as csv to avoid repeat downloads (saves time)
#cycad_data =read.csv("cycad_paleoDB_data.csv")
mid_age=(as.numeric(cycad_data $eag)+as.numeric(cycad_data $lag))/2
cycad_data =cbind(cycad_data,mid_age)


# removing non-cycad taxa and non-foliage organs

cycad_data= cycad_data[!cycad_data$gnl=="Phasmatocycas",]
cycad_data= cycad_data[!cycad_data$gnl=="Cycadospadix",]
cycad_data= cycad_data[!cycad_data$gnl=="Behuninia",]
cycad_data= cycad_data[!cycad_data$gnl=="Dichotozamites",]
cycad_data= cycad_data[!cycad_data$gnl=="Anthrophyopsis",]
cycad_data= cycad_data[!cycad_data$gnl=="Ctenophyllum",]
cycad_data= cycad_data[!cycad_data$gnl=="Sphenozamites",]
cycad_data= cycad_data[!cycad_data$gnl=="Stangerites",]
cycad_data= cycad_data[!cycad_data$gnl=="Zamiostrobus",]
cycad_data= cycad_data[!cycad_data$gnl=="Bucklandia",]
cycad_data= cycad_data[!cycad_data$gnl=="Ctenozamites",]
cycad_data= cycad_data[!cycad_data$gnl=="Leptocycas",]
cycad_data= cycad_data[!cycad_data$gnl=="Dioonites",]
cycad_data= cycad_data[!cycad_data$gnl=="Haitingeria",]
cycad_data= cycad_data[!cycad_data$gnl=="Glossozamites",]
cycad_data= cycad_data[!cycad_data$gnl=="Heilungia",]
cycad_data= cycad_data[!cycad_data$gnl=="Bureja",]
cycad_data= cycad_data[!cycad_data$gnl=="Sarmatiella",]
cycad_data= cycad_data[!cycad_data$gnl=="Eucommiidites",]
cycad_data= cycad_data[!cycad_data$gnl=="Archaeocycas",]
cycad_data= cycad_data[!cycad_data$gnl=="Brunoa",]
cycad_data= cycad_data[!cycad_data$gnl=="Worsdellia",]
cycad_data= cycad_data[!cycad_data$gnl=="Menucoa",]
cycad_data= cycad_data[!cycad_data$gnl=="Bororoa",]
cycad_data= cycad_data[!cycad_data$gnl=="Chamberlania",]
cycad_data= cycad_data[!cycad_data$gnl=="Beaniopsis",]


# setting up parallel resampling

cores=detectCores()
cluster=makeCluster(cores[1]-1)
registerDoParallel(cluster)

start_time2=Sys.time()

boot_data=  cycad_data # run separately for cycad and angio data (update exported file names accordingly below)

n_samples=1e3 # number of samples in re-sampled dataset

n_iterations=1e1 # number of bootstrap iterations (1e3 for final run)

# run resampling

z_out=foreach(i=1: n_iterations,.combine=rbind) %dopar% {

	library(stringr)
	#library(paleobioDB)
	library(mclust)
	library(msir)
	library(stats4)
	library(iNEXT)

	int=25
	intervals=seq(0,350,by=int)

	i=1
	out=NULL
	for (i in 1:(length(intervals)-1)){

	k=1
	k_out=NULL
	for(k in 1:n_samples){
	
	cur_num=round(runif(1,1,dim(boot_data)[1]),digits=0)
	cur_data= boot_data[cur_num,]
	k_out=rbind(k_out,cur_data)
	
	}

	samples=k_out
	cur_interval_start=intervals[i]
	cur_interval_end=intervals[i]+int
	cur_data= samples[samples $mid_age>cur_interval_start,]
	cur_data=cur_data[cur_data$mid_age<cur_interval_end,]
	cur_genera=unique(na.omit(cur_data$gnl))
	cur_raw_genus_richness=length(cur_genera)
	cur_occurrences=length(na.omit(cur_data$gnl))
	cur_sites=length(unique(na.omit(cur_data$cid)))

	j=1
	j_out=NULL
	for (j in 1: cur_raw_genus_richness){
		
		cur_genus=cur_genera[j]
		cur_genus_data=cur_data[cur_data$gnl==cur_genus,]
		cur_genus_counts=length(na.omit(cur_genus_data$gnl))
		j_out=c(j_out, cur_genus_counts)
		
	}
	
	if(max(j_out)>1){TRIPS = doTRiPS_abs(j_out,t=int)}
	if(max(j_out)<2){TRIPS=cbind(c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA))}
	if(max(j_out>0)){chao_richness=ChaoRichness(j_out,conf=0.95)}
	if(max(j_out<1)){chao_richness=c(NA,NA,NA,NA,NA)}
	results=c(cur_interval_start,cur_interval_end,TRIPS[3,1],TRIPS[3,2],TRIPS[3,3], cur_raw_genus_richness,(cur_interval_start+cur_interval_end)/2,paste(sort(cur_genera),collapse=","), cur_occurrences, as.numeric(chao_richness[2]), as.numeric(chao_richness[3]), as.numeric(chao_richness[4]), as.numeric(chao_richness[5]),cur_sites)
	out=rbind(out,results)
	
	}

	
	w=1
	w_out=NULL
	
	for(w in 1:dim(out)[1]){
		
		cur_out=out[w,]
		cur_genera=unlist(strsplit(as.character(cur_out[8]),","))
		if(length(cur_genera)<1){cur_genera=NA}
		cur_results=cbind(w,cur_genera)
		w_out=rbind(w_out,cur_results)
		
	}
	
	unique_genera=na.omit(unique(w_out[,2]))
	
	v=1
	v_out=NULL
	
	for(v in 1:length(unique_genera)){
		
		cur_genus=unique_genera[v]
		cur_data=w_out[w_out[,2]==cur_genus,]
		cur_data=na.omit(cur_data)
		
		if(length(cur_data)>2){
		cur_range=seq(min(as.numeric(cur_data[,1])),max(as.numeric(cur_data[,1])),by=1)
		rangethrough=cbind(cur_range,cur_genus)
		}
		
		if(length(cur_data)<3){
		rangethrough=cur_data	
		}
		
		v_out=rbind(v_out,rangethrough)
		
	}
	
	u=1
	u_out=NULL
	bin_seq=as.numeric(unique(v_out[,1]))
	
	for(u in 1:length(bin_seq)){
		
		cur_bin=bin_seq[u]
		bin_data=v_out[v_out[,1]==cur_bin,]
		
		
		if(length(bin_data)>2){
		cur_genus_richness=length(bin_data[,2])		
		}
		
		if(length(bin_data)<3){
		cur_genus_richness=1
		}
		
		results=c(cur_bin,cur_genus_richness)
		u_out=rbind(u_out,results)
		
	}
	
	out=cbind(na.omit(out),u_out[,2])

	tempMatrix=out
	tempMatrix
}	
	
colnames(z_out)=c("min_age_Ma","max_age_Ma","TRIPS_1","TRIPS_2","TRIPS_3","genus_richness_raw","mid_age_Ma","genera","occurrences","genus_richness","chao_richness","chao_se","chao_lower","chao_upper","sites")
write.csv(z_out,"bootstrapped_angio_data.csv") # UPDATE FILE NAME cycad vs. angio
	
end_time2=Sys.time()
run_time2=end_time2-start_time2
run_time2
stopCluster(cluster)


# compile resampled data per time bin

z_out=read.csv("bootstrapped_angio_data_rangethrough_filtered_2023.03.16.csv") # UPDATE FILE NAME cycad vs. angio

bins=unique(z_out$mid_age_Ma)
y=1
y_out=NULL

for(y in 1:length(bins)){
	
	cur_bin=bins[y]
	cur_data=z_out[z_out$mid_age_Ma==cur_bin,]
	cur_richness_quant=as.numeric(quantile(cur_data$genus_richness,probs=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)))
cur_TRIPS_quant=as.numeric(quantile(cur_data$TRIPS_3,probs=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)))
	cur_Chao_quant=as.numeric(quantile(cur_data$chao_richness,probs=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)))
	cur_abundance_quant=as.numeric(quantile(cur_data$occurrences,probs=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)))
	
	cur_norm_richness= cur_data$genus_richness/cur_data$sites
	cur_norm_TRIPS= cur_data$TRIPS_3/cur_data$sites
	cur_norm_Chao= cur_data$chao_richness/cur_data$sites
	cur_norm_abundance= cur_data$occurrences/cur_data$sites
	
	
	#xbt= # number of taxa w/1st occurrence before interval and last occurrence after
	#xft= # number of taxa w/1st occurrence in this interval and who cross top boundary of interval
	#xbl= # number of taxa 1/last occurrence in this interval 
	
	
	v=1
	cur_bin_genera=unlist(strsplit(as.character(cur_data$genera),","))
	xbts=NULL
	xfts=NULL
	xbls=NULL
	
	for(v in 1:length(cur_bin_genera)){
		
		cur_genus=cur_bin_genera[v]
		prior_bin_data=z_out[z_out$mid_age_Ma>cur_bin,]
		next_bin_data=z_out[z_out$mid_age_Ma<cur_bin,]
		
		if(dim(prior_bin_data)[1]==0 & dim(next_bin_data)[1]>0){

			next_genera=unique(unlist(strsplit(as.character(next_bin_data$genera),",")))
			merged_next=length(next_genera[next_genera==cur_genus])
			if(merged_next>0){
				cur_xbt=0
				cur_xft=1
				cur_xbl=0
				}
			if(merged_next<1){
				cur_xbt=0
				cur_xft=0
				cur_xbl=0
				}


			}
		
		if(dim(next_bin_data)[1]==0){
			cur_xbt=NA
			cur_xft=NA
			cur_xbl=NA
		}
		
		if(dim(prior_bin_data)[1]>0 & dim(next_bin_data)[1]>0){
			
			prior_genera=unique(unlist(strsplit(as.character(prior_bin_data$genera),",")))
			next_genera=unique(unlist(strsplit(as.character(next_bin_data$genera),",")))
			merged_prior=length(prior_genera[prior_genera==cur_genus])
			merged_next=length(next_genera[next_genera==cur_genus])
			if(merged_prior>0 & merged_next>0){
				cur_xbt=1
				cur_xft=0
				cur_xbl=0
				}
			if(merged_prior>0 & merged_next<1){
				cur_xbt=0
				cur_xft=0
				cur_xbl=1
				}
			if(merged_prior<1 & merged_next>0){
				cur_xbt=0
				cur_xft=1
				cur_xbl=0
				}
				
			
		}
		
		xbts=c(xbts,cur_xbt)
		xfts=c(xfts,cur_xft)
		xbls=c(xbls,cur_xbl)
		
	}
	
	xbt=sum(xbts)
	xft=sum(xfts)
	xbl=sum(xbls)
	
	cur_origination=-log(xbt/(xbt+xft))/unique(cur_data$max_age_Ma-cur_data$min_age_Ma) # as in Foote (2000), cf. paleobioDB
	cur_extinction=-log(xbt/(xbt+xbl))/unique(cur_data$max_age_Ma-cur_data$min_age_Ma) # as in Foote (2000), cf. paleobioDB
	
	results=c(cur_bin,cur_richness_quant,cur_TRIPS_quant, cur_abundance_quant,mean(cur_data$genus_richness),2*sd(cur_data$genus_richness)/sqrt(length(cur_data$genus_richness)),mean(cur_data$occurrences),2*sd(cur_data$occurrences)/sqrt(length(cur_data$occurrences)),mean(cur_data$chao_richness),2*sd(cur_data$chao_richness)/sqrt(length(cur_data$chao_richness)),mean(cur_norm_richness),2*sd(cur_norm_richness)/sqrt(length(cur_norm_richness)),mean(cur_norm_TRIPS),2*sd(cur_norm_TRIPS)/sqrt(length(cur_norm_TRIPS)),mean(cur_norm_Chao),2*sd(cur_norm_Chao)/sqrt(length(cur_norm_Chao)),mean(cur_norm_abundance),2*sd(cur_norm_abundance)/sqrt(length(cur_norm_abundance)),mean(cur_data$sites), cur_origination, cur_extinction,length(cur_data$occurrences))
	y_out=rbind(y_out,results)

	
}

write.csv(y_out,"bin_out_angio.csv") # UPDATE FILE NAME cycad vs. angio


# plot data 

angio_bin_out=read.csv("bin_out_angio_mod.csv")
cycad_bin_out=read.csv("bin_out_cycad_mod.csv")

par(mfrow=c(2,1))
par(mar=c(0,0,0.5,0))
par(oma=c(4,5,1,1))
par(xpd=F)

plot(cycad_bin_out[,51]~ cycad_bin_out[,2],cex=0,ylim=c(0.0001,1),xlim=c(359,0),axes=F,log="y")
rect(359,0.0006,0,0.00015,col="grey85",border=F)
lines(cycad_bin_out[,51]~ cycad_bin_out[,2],lty=2,col="darkblue") # cycad origination rate
lines(cycad_bin_out[,52]~ cycad_bin_out[,2],col="darkblue") # cycad extinction rate
points(cycad_bin_out[,51]~ cycad_bin_out[,2],pch=21,bg="white",col="darkblue") # cycad origination rate
points(cycad_bin_out[,52]~ cycad_bin_out[,2],pch=21,bg="darkblue",col="darkblue") # cycad extinction rate
lines(angio_bin_out[,51]~ angio_bin_out[,2],lty=2,col="darkred") # angio origination rate
lines(angio_bin_out[,52]~ angio_bin_out[,2],col="darkred") # angio extinction rate
points(angio_bin_out[,51]~ angio_bin_out[,2],pch=21,bg="white",col="darkred") # angio origination rate
points(angio_bin_out[,52]~ angio_bin_out[,2],pch=21,bg="darkred",col="darkred") # angio extinction rate
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.001,0.01,0.1),labels=c("0.001","0.01","0.1"))
mtext(side=2,line=3,"origination / extinction")
text(20,0.0003,"zero",cex=0.9)

plot(angio_bin_out[,4]/max(angio_bin_out[,5])~ angio_bin_out[,2],xlim=c(359,0),cex=0,ylim=c(0.0045,1.75),log="y",axes=F,yaxs="i")
rect(cycad_bin_out[,2]+12.5,0.004,cycad_bin_out[,2]-12.5,cycad_bin_out[,30]/max(cycad_bin_out[,30]),border="darkblue",col="lightblue") # cycad relative abundance
arrows(cycad_bin_out[,2], cycad_bin_out[,25]/max(cycad_bin_out[,30]), cycad_bin_out[,2], cycad_bin_out[,35]/max(cycad_bin_out[,30]),length=0.02,angle=90,code=3,col="darkblue") # cycad relative abundance, 5-95% confidence interval error bars
rect(angio_bin_out[,2]+12.5, 0.004, angio_bin_out[,2]-12.5, angio_bin_out[,30]/max(angio_bin_out[,30]),col=rgb(1,0,0,0.2),border="darkred") # angio relative abundance
arrows(angio_bin_out[,2], angio_bin_out[,25]/max(angio_bin_out[,30]), angio_bin_out[,2], angio_bin_out[,35]/max(angio_bin_out[,30]),length=0.02,angle=90,code=3,col="darkred") # angio relative abundance, 5-95% confidence interval error bars
axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.01,0.1,1),labels=c("0.01","0.1","1"))
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
box(lwd=1)
mtext(side=1,line=2,"age (Myr ago)")
mtext(side=2,line=3,"relative occurrences")


# supplement

par(mfrow=c(1,1))
par(xpd=F)
plot(angio_bin_out[,4]/max(angio_bin_out[,5])~ angio_bin_out[,2],xlim=c(350,0),cex=0,ylim=c(0.005,70),log="y",axes=F,xlab=NA,ylab=NA)

arrows(65,1,65,0.001,length=0,lty=1)
points(65,0.01,pch=22,cex=4,bg="white",col="white")
text(65,0.01,"K-Pg")

lines(cycad_bin_out[,36]/max(cycad_bin_out[,36])~ cycad_bin_out[,2],col="darkblue",lty=1,lwd=1.5) # cycad resampled genus richness
lines(cycad_bin_out[,19]/max(cycad_bin_out[,19])~ cycad_bin_out[,2],col="darkblue",lty=2,lwd=1.5) # cycad TRIPS genus richness
lines(cycad_bin_out[,40]/max(cycad_bin_out[,40])~ cycad_bin_out[,2],col="darkblue",lty=3,lwd=1.5) # cycad Chao Richness

lines(angio_bin_out[,36]/max(angio_bin_out[,36])~ angio_bin_out[,2],col="darkred",lty=1,lwd=1.5) # angio resampled genus richness
lines(angio_bin_out[,19]/max(angio_bin_out[,19])~ angio_bin_out[,2],col="darkred",lty=2,lwd=1.5) # angio TRIPS genus richness
lines(angio_bin_out[,40]/max(angio_bin_out[,40])~ angio_bin_out[,2],col="darkred",lty=3,lwd=1.5) # angio Chao Richness

axis(side=2,las=1,tck=-0.02,mgp=c(3,0.5,0),at=c(0.01,0.1,1))
axis(side=1,las=1,tck=-0.02,mgp=c(3,0.5,0))
#box(lwd=1)
mtext(side=1,line=2,"age (Myr ago)")
mtext(side=2,line=2.5,"relative genus richness",at=0.1)


legend(275,30,c(expression(bold("Cycadales")),"rangethrough","TRIPS","Chao",expression(bold("Angiospermae")),"rangethrough","TRIPS","Chao"),lty=c(NA,1,2,3,NA,1,2,3),col=c(NA,"darkblue","darkblue","darkblue",NA,"darkred","darkred","darkred"),box.lwd=NA,ncol=2)

