#!/usr/bin/env Rscript

usage = "
--------------------------------------------------
Calculate sv_snp-snp_mid num and draw.
--------------------------------------------------
Usage:
	Rscript <Program> <mode> <vcftools.LD> <size> <outdir>
--------------------------------------------------
OPTIONs:
	mode:
		0 ---> only out put LD type
		1 ---> only draw heatmap and linesfig
		2 ---> out LD type and draw figs
	haploview.LD:
		LD file generated from the main program
	size:
		flanking ref size used in the main program
    outdir:
        path of the outputs
--------------------------------------------------
"

argv <- commandArgs(TRUE)

if(length(argv) != 4){stop(paste0("\nWrong num of options!\n",usage))}

mode=argv[1]
if(mode %in% c(0,1,2)){
    input=argv[2]
    size=as.numeric(argv[3])+1
    outdir=argv[4]
}else{
    stop("[ERROR]: Wrong mode option!")
}

# ? stringsAsFactors=F
# parse only siteA siteB R2 field
tt=read.table(input,header=T,stringsAsFactors=F)[,c(2,3,5)]
names(tt)=c("L1","L2","r.2")
tt$L1=as.numeric(tt$L1)
tt$L2=as.numeric(tt$L2)
tt$r.2=as.numeric(tt$r.2)
tt$r.2[is.na(tt$r.2)] = 0
# ! head(tt)
# ? creat a blank matrix
cc=(unique(c(tt$L1,tt$L2)))
# instead of check the input LD item, just fill the matrix and fill NA as 0
# cc_num=length(cc)
# if(cc_num != 2*size-1){stop(paste0("\nInput variant number ",cc_num," is not equal to ",2*size-1," !\n",usage))}
mm_size=2*size-1
mm.r2=matrix(ncol=mm_size,nrow=mm_size)
# ? fill the matrix based on the LD values (up/low triangular)
for (i in 1:nrow(tt)){mm.r2[tt$L1[i],tt$L2[i]]=tt$r.2[i];}
for (i in 1:nrow(tt)){mm.r2[tt$L2[i],tt$L1[i]]=tt$r.2[i];}
# ? fill diagonal and NA cell to 0
mm.r2[is.na(mm.r2)]=0
# ? get no-query-record matrix (query was in the middle of the matrix)
mm.snp=mm.r2[-size,-size]
mm.ss=mm.snp
# ? sort the matrix rows decreasingly
for (i in 1:nrow(mm.ss)){mm.ss[i,]=sort(mm.ss[i,],decreasing=T);}
# ? calc the median for each row
snp_mid=sapply(1:nrow(mm.ss),function(x){ median(sort(mm.ss[,x],decreasing=T))})
# ? extract the comparison between query and ref
sv_snp=sort(as.numeric(mm.r2[size,-size]),decreasing=T)
# ? get sv id from input name
sv=gsub("_flank.*$","",input)
sv=gsub(".*/","",sv)
# compare and get the type
if(mode != 1){
cmp=data.frame(snp_mid,sv_snp,cc=sv_snp-snp_mid)
cmp$cc[cmp$cc>0]=1
cmp$cc[cmp$cc<=0]=0
cmp.sum=sum(cmp$cc)
if(cmp.sum<nrow(cmp)/3){
	cmp.type="low"
}else if(cmp.sum>nrow(cmp)*2/3){
	cmp.type="high"
}else{
	cmp.type="mid"
}
cmp.out=data.frame(sv=sv,num=cmp.sum,total=nrow(cmp),type=cmp.type)
# output to file
ldtype_out=paste0(outdir,"/",sv,"_ldtype.txt")
write.table(cmp.out,file=ldtype_out,quote=F,col.names=F,row.names=F,sep="\t")
write.table(cmp.out,file=stdout(),quote=F,col.names=F,row.names=F,sep="\t")
}


if(mode != 0){
	# draw heatmap
	library(lattice)
    heatmap_pdf=paste0(outdir,"/",sv,"_heatmap.pdf")
	pdf(heatmap_pdf)
	print(levelplot(mm.r2,col.regions=colorRampPalette(c("grey90", "black"), space = "rgb")(100),scales=list(y=list(at=NULL),x=list(at=NULL)),main="", xlab="", ylab=""))
	dev.off()
	write(paste0("heatmap Fig: ",heatmap_pdf,"\n"),stderr())
	# draw line plot
	col.snp="grey"
	col.mid="black"
	col.sv="red"
	lwds=1
    line_pdf=paste0(outdir,"/",sv,"_lines.pdf")
	pdf(line_pdf)
	plot(c(snp_mid,1),type="n",lwd=lwds*2,xaxt = 'n', yaxt = 'n',main = "", sub = "", xlab = "",  ylab = "")
	for(i in 1:nrow(mm.ss)){
		lines(as.numeric(mm.ss[i,]),col=col.snp,lwd=lwds)
	}
	lines(snp_mid,col=col.mid,lwd=lwds*2)
	lines(sv_snp,col=col.sv,lwd=lwds*2)
	# axis
	axis(1,at=c(1,nrow(mm.ss)/2,nrow(mm.ss)),label=c(1,"r2 rank",nrow(mm.ss)),cex.axis=1.5)
	axis(2,at=c(0,0.5,1),label=c(0,"r2 value",1),cex.axis=1.5)
	dev.off()
    write(paste0("lines Fig: ",line_pdf,"\n"), stderr())
}

