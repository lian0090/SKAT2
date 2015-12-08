
`[.gds.class` <- 
    function(x, samples,snps)
{
       if(missing(samples)){
       	sample.id=NULL
       }else{
       	if(is.numeric(samples)){
       		sample.id=read.gdsn(index.gdsn(x,"sample.id"))[samples]
       	}else if (is.character(samples)) sample.id=samples
       	}       	
       if(missing(snps)){
       	snp.id=NULL
       }else{
        if(is.numeric(snps)){	
       	snp.id=read.gdsn(index.gdsn(x,"snp.id"))[snps]
       	}else if (is.character(j)) snp.id=snps
       }
    return(snpgdsGetGeno(x,sample.id=sample.id,snp.id=snp.id))
}

