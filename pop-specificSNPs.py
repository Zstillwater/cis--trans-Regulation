#!/usr/bin/env python 
import sys
import argparse

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass


def main():
    
    args=get_args()
    snplist=args.INPUT

    (samplesInGroup,groupOfSample)=readGroup(args.g)
    out=open("result.groupAllele.%s_%s.xls" % (args.Q1,args.Q2),"w")
    sample=[]
    out.write( "#Chr\tPos\t%s\t%s\n" % (args.Q1,args.Q2) )

    (t_number,q1_number,q2_number,ud_number)=(0,0,0,0)
    homAllen=["AA","GG","CC","TT"]

    for line in open(snplist,"r"):
        info=line.strip().split()
        if line.startswith("#"):
            sample = info[2:]
            #print sample
            checkGroup(sample,groupOfSample)
        else: 
            cid=info.pop(0)
            pos=info.pop(0)
            StatQ1 = statQ(info,sample,samplesInGroup[args.Q1])
            StatQ2 = statQ(info,sample,samplesInGroup[args.Q2])
            #print StatQ1,StatQ2
            if args.hom:
                if StatQ1[0] not in homAllen or StatQ2[0] not in homAllen:
                    continue
            if "N" in StatQ1[0]  or "N" in StatQ2[0]:
                continue

            if StatQ1[0] == StatQ2[0] : 
                continue

            if  StatQ1[1] < 3 or float(StatQ1[1])/StatQ1[2] <= 0.66  or  StatQ2[1] < 3 or float(StatQ2[1])/StatQ2[2] <= 0.66:
            #print StatQ1,StatQ2,homAllen
            #if float(StatQ1[1])/len(samplesInGroup[args.Q1]) < args.qr or float(StatQ2[1])/len(samplesInGroup[args.Q2]) < args.qr:
                continue
            t_number += 1
            out.write(  "%s\t%s\t%s(%s/%s)\t%s(%s/%s)\n" % (cid,pos,StatQ1[0],StatQ1[1],StatQ1[2],StatQ2[0],StatQ2[1],StatQ2[2]))
    out.close()


def statQ(info,sample,Q):
    # return the max allen
    # stat=["none","none"]
    count={}
    Q=(getGene(info,sample,Q))[0]
    stat=judgeMin(Q)  
    return stat


def getGene(info,sample,Q):
    gQ=[]
    dgQ={}
    for i in range(len(info)):
        if sample[i] in Q:
            gQ.append(info[i])
            dgQ[sample[i]] = info[i]
    return gQ,dgQ


def  checkGroup(sample,groupOfSample):
    for sampleName in groupOfSample.keys():
        if sampleName not in sample:
            print "#Error! Sample ID %s in group file not in snplist file !" % (sampleName)
            exit(1) 


def readGroup(file):
    group={}
    sample={}
    for line in open(file,"r"):
        if line.startswith("#"):
            continue
        info=line.strip().split()
        if info[1] not in group:
            group[info[1]]=[]
        group[info[1]].append(info[0])
        sample[info[0]] = info[1]
    return (group,sample)



def accumulateDict(myDict,val,keyA,keyB="no"):
    '''
        accumulate a 1~2d Dictionary  
    '''
    if keyB == "no":
        if keyA in myDict:
            myDict[keyA] += val
        else:
            myDict[keyA] = val
    else:
        if keyA in myDict:
            if keyB in myDict[keyA]:
                myDict[keyA][keyB] += val
            else:
                myDict[keyA].update({keyB: val})
        else:
            myDict.update({keyA:{keyB: val}})

    return myDict


def judgeMin(mylist):
    count={}
    allnum=0
    for ss in mylist:
        if ss=="NN":
            continue
        allnum+=1
        if ss not in count:
            count[ss]=1
        else:
            count[ss]+=1
    if allnum == 0:
        return(["NN",0,0])
    s = list(sorted(count.items(), key=lambda x: x[1],reverse=True)[0])
    s.append(allnum)
    return s


def get_args():

    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''

'''
    )
    parser.add_argument('-g',metavar='group',help='input of group file',type=str)
    parser.add_argument('-Q1',metavar='Q1',help='name of group1',type=str)
    parser.add_argument('-Q2',metavar='Q2',help='name of group2',type=str)
    parser.add_argument('-qr',metavar='qr',help='min-rate to define the SNP type of Q1 and Q2 ',default="0.8",type=float)
    parser.add_argument('-hom',help='only keep hom allen for parents',action='store_true',default=False)
    parser.add_argument('INPUT',metavar='input',help='cmd input vcf file must be given')


    args = parser.parse_args()
    return args



if __name__ == '__main__':
    main()
