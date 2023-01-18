import os,time,ConfigParser
def main():
		
	#	dataFileName ='~/SimingLab/library/diffdriver_anno/hmmOGpar_ASHmean.rds'
		for i1 in ['0','1','1.5']:
			for i2 in ['0','0.5','1']:
				for i3 in ['200','400','600','800','1000','1200']:
					identifier='pppower_betaf0=' + i1 + '_betagc=' + i2 + '_sample' + i3
					outfiles =  identifier + '.Rd'
					RFileName = 'phynotype'+identifier  +'.R' 
					RFile = open(RFileName, 'w')
					RFile.write('library("devtools")\n')
					#RFile.write('library("foreach")\n')
					#RFile.write('library("doParallel")\n')
					#RFile.write('registerDoParallel(cores=4)\n')
					RFile.write('load_all("../../")\n')
					RFile.write('i1='+i1+'\n')
					RFile.write('i2='+i2+'\n')
					RFile.write('i3='+i3+'\n')					
               				#RFile.write('hotspot=readRDS(file="'+dataFileName+ '")\n')
					#RFile.write('hmm[9]=0\n')
					RFile.write('Nsim=200\n') 
					RFile.write('set.seed(10)\n')               			
					RFile.write('simuresdiff=power_comparediff(binary=TRUE,Nite=Nsim,sgdata=sgdata,bmrpars=log(BMR),Nsample=i3,betaf0=i1,beta_gc=c(i2,2.2),para=c(0.8,0.2),hotseq=hotseq2,hmm=hmm)\n')
					RFile.write('save(simuresdiff,file="'+outfiles +'")\n')
					RFile.write('q("no")\n')
					RFile.close()
					shFileName = 'phynotype'+identifier  +'.sh'
					shFile = open(shFileName, 'w')
					shFile.write('#!/bin/bash\n')
					shFile.write('#SBATCH --job-name=dd'+i1+'_'+i2+'_'+i3+'\n')
					shFile.write('#SBATCH --nodes=1\n')
					shFile.write('#SBATCH --ntasks-per-node=1\n')
					shFile.write('#SBATCH --time=30:00:00\n')
					shFile.write('#SBATCH --mail-type=BEGIN,END,FAIL\n')
					shFile.write('#SBATCH --output=real.out\n')
					shFile.write('#SBATCH --error=real.err\n')
					shFile.write('time R CMD BATCH '  + RFileName + ' a' +  identifier  +'.out\n')
					shFile.close()
					os.system ('sbatch ' + shFileName)
					#time.sleep(1) #delay 1/10 of a second between job submissions

if __name__=="__main__":
	main()
	
