#!/usr/bin/python
import sys
import re
def main():
	try :
		Flag=sys.argv[1]
		VCF=sys.argv[2]
	except: 
		help()
	ReadVCF(VCF,Flag)

def ReadVCF(vcf,Flag):
	if Flag =='SNP':
		LIST=open("snp.list","w")
		INFO=open("snp.info","w")
	else :
		LIST=open("indel.list","w")
		INFO=open("indel.info","w")
	for line in open(vcf):
		if "##" in line:
			continue
		elif "#" in line:
			sample=str.split(line)[9:]
			if Flag=='SNP':	
				LIST.write("#Chr\tPos\t"+'\t'.join(sample)+"\n")
				INFO.write("#Chr\tPos\tRef\tAlt\tDeepMean\tCallingRate\tRefRate\tAltRate\tRR-Rate\tRA-Rate\tAA-Rate\n")
			else :
				LIST.write("#Chr\tPos\tRef\tAlt\t"+'\t'.join(sample)+"\n")
				INFO.write("#Chr\tPos\tLength\tRef\tAlt\tDeepMean\tCallingRate\tRefRate\tAltRate\tRR-Rate\tRA-Rate\tAA-Rate\n")
		else:
			arry=str.split(line) 
			chrom,pos,ref,alt=arry[0],arry[1],arry[3],arry[4]
			length=str(len(alt)-len(ref))
			genos=arry[9:]
			if Flag =='SNP':
				LIST.write(chrom+"\t"+pos)
				INFO.write(chrom+"\t"+pos+"\t"+ref+"\t"+alt)
			else :
				LIST.write(chrom+"\t"+pos+"\t"+ref+"\t"+alt)
				INFO.write(chrom+"\t"+pos+"\t"+length+"\t"+ref+"\t"+alt)
			Allen={'0':ref,'1':alt}	
			Total=Calling=Deep=Type=AltCount=0
			Genotype={0:0,1:0,2:0}
			for geno in genos:
				Total += 1
				info=re.split(':|/|,',geno)	
				if info[0] == '.' or info[4] == '.':
					if Flag =='SNP':	
						LIST.write("\tNN")
					else :
						LIST.write("\t./.")
					continue
				Calling += 1
				Deep += int(info[4])
				Type=int(info[0])+int(info[1])
				AltCount += int(Type)
				Genotype[int(Type)] += 1
				if Flag =='SNP':	
					LIST.write("\t"+Allen[info[0]]+Allen[info[1]])
				else:
					LIST.write("\t"+info[0]+"/"+info[1])
			LIST.write("\n")		
			CallingRate= '%.4f'%(float(Calling)/Total)
			DeepMean='%.2f'% (Deep/float(Calling))
			AltRate='%.4f'% (0.5*AltCount/Calling)
			RefRate=1-float(AltRate)
			RRrate='%.4f'% (Genotype[0]/float(Calling))
			RArate='%.4f'% (Genotype[1]/float(Calling))
			AArate='%.4f'% (Genotype[2]/float(Calling))
			txt="\t"+str(DeepMean)+"\t"+str(CallingRate)+"\t"+str(RefRate)+"\t"+str(AltRate)+"\t"+str(RRrate)+"\t"+str(RArate)+"\t"+str(AArate)+"\n"
			INFO.write(txt)

def help():
	print "\tUsage: python ./ExtractSNPvcf.py SNP/INDEL [your SNPvcf file]"
	sys.exit(0)

if __name__ == "__main__":
	main()
