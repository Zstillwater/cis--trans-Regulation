#!/usr/bin/env python 
import sys

def main():

    fasta = readFasta(sys.argv[1])

    variation = readVcf(sys.argv[2],sys.argv[3]+".Heterozygous.mark.xls")

    of1=open(sys.argv[3]+".change.fasta","w")
    of2=open(sys.argv[3]+".change.record","w")

    for contig in sorted(fasta.keys()):

        contig_list = list("".join(fasta[contig]))

        if contig not in variation:
            of1.write(">%s_Presude\n%s\n" % (contig,"".join(contig_list)))
            continue

        for site in variation[contig]:
            if len(site[1]) <= len(site[2]):    # snp or insert  
                of2.write("Reference %s:%s is %s, change %s -> %s\n" %(contig,site[0],contig_list[int(site[0])-1],site[1],site[2]))
                contig_list[int(site[0])-1] = site[2]

            if len(site[1]) > len(site[2]):     # del
                ref = contig_list[int(site[0])-1]
                for x in range(1,len(site[1])):
                    ref += contig_list[int(site[0])+x-1]
                    contig_list[int(site[0])+x-1] = ""
                of2.write( "Reference %s:%s is %s, change %s -> %s\n" %(contig,site[0],ref,site[1],site[2]))
        of1.write(">%s_Presude\n%s\n" % (contig,"".join(contig_list)))

    of1.close()
    of2.close()



def readVcf(file,outfile):

    var = {}
    of  = open(outfile,"w")
    of.write("#Chr\tPos\tRef\tAlt\n")

    for line in open(file):
        if line.startswith("#"):
            continue
        else:          
            info = line.strip().split()
            alt = info[3]
            if info[0] not in var:
                var[info[0]]=[]
            var[info[0]].append([info[1],info[2],info[3]])
        of.write("%s\t%s\t%s\t%s\n" %(info[0],info[1],info[2],info[3]))

    return var 



def readFasta(file):

    print "# Reading fasta dict from %s" %(file)
    fasta    = {}
    fasta_id = ''
    for line in open(file):
        if line.startswith(">"):
            fasta_id = line.strip().replace(">","")
            fasta[fasta_id] = []
        else:
            fasta[fasta_id].append(line.strip())

    return fasta 



if __name__ == '__main__':
    main()
