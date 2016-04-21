#!/usr/bin/python
#Oct 12, 2015
#Make FASTA file from gtf and genome sequence file
#

#Imports
import re
import sys, time
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
########## Support Functions #########

  
######### Functions #########

def open_genome(gen):
  handle=open(gen, "rU")
  record_dict=SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
  #for record in SeqIO.parse(handle, "fasta"):
  #  print record
  handle.close()
  return record_dict

def get_gene_name(attr,ftype): ##Needs work to make it usable in lots of circumstances
  ret={}
  if ftype=='gtf': #use gene_id
    for att in attr.split(";"):
      if att:
        try:
          att=att.strip()
          key, value = att.split(' ') #
          key,value=key.replace('"',''),value.replace('"','')
          ret[key]=(value)
        except ValueError:
          sys.exit('\033[31mNot a standard GTF file format!\033[0m...')
    try:
      return ret['gene_id']
    except IndexError:
      sys.exit('\033[31mNot a standard GTF file format!\033[0m...')
  elif ftype=='gff': #use parent or name (protein_id for ensembl gff3)
    for att in attr.split(";"):
      try:
        key, value = att.split("=") #Works for GFF3 but not GTF (use ' ' instead of '=')
        key,value=key.replace('"',''),value.replace('"','')
        ret[key]=(value)
      except ValueError:
        sys.exit('\033[31mNot a standard GFF file format!\033[0m...')
    try:
      if 'protein_id' in ret: #ensembl gff3 format
        return ret['protein_id']
      elif 'Parent' in ret:
        return ret['Parent']
      elif 'Name' in ret:
        return ret['Name']
    except IndexError:
      sys.exit('\033[31mNot a standard GFF format or non-standard attributes!\033[0m...')
  #Alternative:name=attr[attr.lower().find(x)+len(x):].strip('"').strip('=').strip() #slice attribute by word

def gtfingene(gtffile,filetype):
  GTF_Genes = {}
  with open(gtffile) as f:
    for line in f:
      if line[0]=="#":
        pass
      else:
        Chrm,Meth,Type,Start,End,score,strand,frame,attribute = [n for n in line.strip().split('	')]
        if Type=='CDS' or Type=='exon':
          #attr=[n for n in attribute.split(';')] #creates list of all attribute pieces split by ;
          #Key=attr[0].split('=')[1]#Key = attribute.split(';')[5].split('"')[1].strip('"')+'.1'
          Key=get_gene_name(attribute,filetype)
          if Key in GTF_Genes.get(Chrm,{}):
            GTF_Genes.setdefault(Chrm, {})[Key] = ((GTF_Genes[Chrm][Key])+(int(Start),int(End),strand))
          else:
            GTF_Genes.setdefault(Chrm, {})[Key] = (int(Start),int(End),strand)
        else:
          pass
  print '\033[92mTotal Chromosomes Added\033[0m.....', len(GTF_Genes)
  print '\033[92mTotal Gene Lines Added\033[0m....\033[0m',sum([len(x) for x in GTF_Genes.values()])
  #print GTF_Genes
  return GTF_Genes
  
def Extractor(Genome,GeneDict):
  '''Takes the Gene Dict and the genome to output a new dictionary of SeqRecords for genes'''
  Flipper,curchrm = {},''
  for k, v in Genome.items():
    newchrm=k
    #if newchrm!=curchrm:
      #print newchrm
    try:
      Sequence,Genes,curgene = v.seq.upper(),GeneDict[k],[] #get gene info for current chromosome
      for gname, locus in Genes.items(): #for each gene (which has multiple exons)
        for i in range(0,len(locus),3): #Go through exons
        #print gname,locus[i],locus[i+1],locus[i+2]
          coords1,coords2,strand = locus[i],locus[i+1],locus[i+2]
          if i == 0: #first time through loop
            GeneSeq = SeqRecord(Sequence[(coords1-1):coords2],name=gname)
          #print gname,len(Sequence[(coords1-1):coords2]),len(GeneSeq.seq),GeneSeq.seq[0:3],GeneSeq.seq[-3:]
          else: #not first iteration
            GeneSeq = GeneSeq+Sequence[(coords1-1):coords2]
        if strand == "+": #deal with strand info after looping through exons
          Flipper[gname] = GeneSeq
        elif strand == '-':
          Reverse = GeneSeq.reverse_complement()#Get RvComp of SeqRecord
	  Flipper[gname] = Reverse
        else:
          print "Incorrect strand info!!!"
    except KeyError:
      pass
    curchrm=k
  print '\033[92mGenes Indexed\033[0m...',len(Flipper)
  #print Flipper,[v.seq for k,v in Flipper.items()]
  return Flipper

def get_protein(record): #not currently being used
  '''Extract amino acid sequence from SeqRecord'''
  minlen=1
  for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
    prolist=[]
    for frame in range(3):
      length=3*((len(record)-frame)//3)
      for pro in nuc[frame:frame+length].translate().split("*"):
        if len(pro) >= minlen: #and pro[0]=='M':
          info=(pro,len(pro),strand,frame)
          prolist.append(info)
          #print pro[:30],"...",pro[-3:]," - ",len(pro),strand,frame
    for n in prolist:
      a=re.findall('M.*',str(n[0]))
      if a:
        a=Seq(a)
        print a
      #for a,b,c,d in n:
        #print a,b,c,d
  #print prolist

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer

def Filter_proteins(proteins):
  error=open("Error.protein.list.txt","w")
  prot={}
  for k,v in proteins.items():
    orf_list=find_orfs_with_trans(v.seq,1,30)
    l,store=0,[]
    for start,end,strand,pro in orf_list:
      a=re.findall('(?=(M.*))',pro) #Get proteins from start codon
      if a:
        store=store+a
      #print("%s...%s - length %i, strand %i, %i:%i" % (pro[:30], pro[-3:], len(pro), strand, start, end))
    if len(store)==1:
      prot[k]=(store[0],len(store[0]))
    elif len(store)>1:
      p=max(store,key=len)
      prot[k]=(p,len(p))
    else: #empty protein list
      towrite=k+str(store)
      error.write(towrite)
      #print "Check length of protein list!!!",k,v,len(store)
      pass
  error.close()
  return prot

def Load_splices(splfile):
  '''Load splice file data, to be incorporated into annotations'''
  splices = []
  with open(splfile) as f:
    for line in f:
      #tsv 5col +2/dataset
      if "gene_name" in str(line): #header line
        pass
      else:
        fline = [n for n in line.strip().split('	')]
        #does not currently deal with unknown number of data fields
        #gid,score,jnum,jtype,a1,n1,a2,n2,a3,n3,a4,n4,a5,n5,postype=fline
        splices.append(fline) #append tuple of line to list
  return splices

def calc_min(alist,anum): ##in progress
  minl=[]
  for x in alist:
    n=x[0]
    minl.append(n-anum)
  return min(minl)

def check_junc(n1,n2,jtype): ###Haven't implemented this function yet (Jan 13 '16)
  '''Check to make sure junc numbers make sense'''
  diff=abs(n1-n2)
  if jtype=='IR' or jtype=='ALT':
    if diff==1:
      pass
    else:
      print "Error in junction numbering, should be == 1"
  elif jtype=='SK':
    if diff>1:
      pass
    else:
      print "Error in junction numbering, should be >1 for SKIP"    

def deal_with_splice(alist,exons):
  '''Takes splice event list and list of exons (tuples) and modifies exons involved in AS'''
  unk,newex={},[]
  print alist[0]
  if 'IR' in alist[3]: #intron retention
    junc=alist[2].split('-') #split junction number ie."2-3"
    #print exons
    if exons[0][2]=='-': #neg strand gene
      junc1,junc2=int(junc[1])-1,int(junc[0])-1
    else: #
      junc1,junc2=int(junc[0])-1,int(junc[1])-1
    j1,j2=exons[int(junc1)][0],exons[int(junc2)][1] 
    new=(j1,j2,exons[0][2]) #new exon tuple
    #print alist,'\n',exons,'\n',new #de-bugging#
    pos1,pos2=int(junc[0])-1,int(junc[1])
    del exons[pos1:pos2]
    exons.insert(pos1,new) #insert new exon tuple where two exons used to exist
    #print exons #de-bug purposes#
    #sys.exit(0) #de-bug purposes#
    return exons #returns whole list of exons with modification
  elif 'ALT' in alist[3]: #ALTA/ALTD/ALTP
    #ATLA left coordinate is new, ALTD right coord is new
    junc=alist[2].split('-')
    if 'ALTA' in alist[3]:
      print "ALTA",junc
    #  j1,j2=exons[int(junc[0])-1][0],exons[int(junc[1])-1][1]
    #  new=(j1,j2,exons[0][2]) #new exon tuple
    #  pos1,pos2=int(junc[0])-1,int(junc[1])
    #  del exons[pos1:pos2]
    #  exons.insert(pos1,new) #insert new exon tuple where two exons used to exist
    #  return exons #returns whole list of exons with modification
    elif alist[3].startswith("ALTD"):
      print "ALTD",junc
      #
    elif alist[3].startswith("ALTP"):
      print "ALTP",junc
      #
  elif alist[3].startswith("SK"): #SKIP SKAT SKAD etc
    junc=alist[2].split('-')
    
  else:
    unk[alist[0]]=(alist[1:])

###How to deal with complex events? SKAT SKAD SKAP CREX etc

def gff_plus_splice(spl,annt):
  '''Compbine splice data and annotation to modify annotations'''
  #splice info doesn't have chrm, only gene names
  newgff,glist,genes={},[],[]
  for k,v in annt.iteritems(): #for each chromosome
    newgff[k]={} #initialize nested dict, is this correct?
    for l,m in v.iteritems(): #for each gene on that chrm, m=list of all exons (start,end,strand)
      for i in spl: #go through list of splice data
        if l==i[0]: #if gene names match, modify annotation
          glist=[m[t:t+3] for t in range(0,len(m),3)] #turns flat list into list of tuples for exons
          ##Not all splice descriptions have position info
          #switch exon order for '-' strand genes
          if glist[0][2]=='-': #neg strand
            glist=glist[::-1] #remember: newex will be in reverse order
          newex=deal_with_splice(i,glist) #input splice line and exon tuples, gives new modified exon
      #    print i,'\n',newex
          if newex: #if new exons returned from "deal_with_splices"
            if l in newgff[k]:
              #create dict values with lists of tuples for each transcript isoform
              newgff[k][l]=newgff[k][l]+[newex] #same order as splice event info
              #print "Already in newgff",'\n',newgff[k][l]
            else:
              newgff[k][l]=[newex]
              #print "First time in newgff",'\n',newgff[k][l]

###(Jan15) Need to finish (below) to add new exons to list and then into nested dict
 #         if [item for item in genes if item[0] == l]: #if gene in list, name should be first item of tuple
 #           newname=str(nm[0])+'_'+str(int(nm[1])+1)
 #           mytup=([newname]+newex) #add lists inside tuple
 #           genes.append(mytup)
 #         else: #add to list for first time
 #           newname=str(l)+'_'+str(1) #start naming gene with '_1' '_2' etc
 #           try:
 #             mytup=([newname]+newex) #add name and list of exons into one tuple
 #           except TypeError:
 #             print newname,'\n',newex,'\n',i
 #             sys.exit(0)
 #           genes.append(mytup)
 #         newgff[k][l]=(newname,newex) #
          #

          #print i[0],i[-1],type(i[-1].split(":")[0])
         # mynum1,mynum2=i[-1].split(":")[0],i[-1].split(":")[1] #location of AS event beg./end
         # exonid=min([x[0] for x in glist], key=lambda x:abs(x-mynum1)) #find nearest match to v[0] in list of tuples
         # print exonid
          ##print l,m,'\n',i
          #have one line slice info (i) and list (m) of all gene annotation info

          #exonid=min([x[0] for x in m], key=lambda x:abs(x-mynum))
          #m[exonid]
          #newgff[k][l]=
          #tmplist.append(i) #append spl info to temp list
    ###at end of each genes' loop, modify annotation
  print "gff_plus_splice dict"
  for k,v in newgff.iteritems():
    if v:
      print k,":",v

def Output(adict,output):
  outfile=open(output,'w')
  for k,v in adict.items():
    try:
      myline=">"+k+"_length="+str(v[1])+"\n"+v[0]+"\n"
      outfile.write(myline)
    except TypeError:
      print k,v

if __name__ == '__main__':
  starttime=time.time()
  #Intialized variables...
  filetype=sys.argv[1]
  thegenome = sys.argv[2]
  gtffile = sys.argv[3]
  splicefile = sys.argv[4]
  outputfile=sys.argv[5]
  #Calls
  #print '\033[93mReading in GTF file\033[0m...'
  GENES = gtfingene(gtffile,filetype)
  #print '\033[93mReading In Genome\033[0m...'
  Genome=open_genome(thegenome)
  #print '\033[93mExtracting Genes from Genome\033[0m...'
  seqs=Extractor(Genome,GENES)
  pro=Filter_proteins(seqs)
  print '\033[93mRun analysis for splice-modified annotations\033[0m...'
  splices=Load_splices(splicefile)
  gff_plus_splice(splices,GENES) #run Extractor and Filter_proteins on output gff of spliced modified txn
  #Output protein sequences to file and logs errors to another file
  Output(pro,outputfile)
  ##Need another output for sequences after splice modifications
  endtime=(time.time()-starttime)
  print "Script run-time (sec): ",endtime

'''
for k, v in seqs.items():
  print k#,v.seq
  orf_list=find_orfs_with_trans(v.seq,1,1)
  #report=max(len(x) for s,e,d,x in orf_list)
  l,store=0,[]
  for start,end,strand,pro in orf_list:
    a=re.findall('(?=(M.*))',pro)
    if a:
      #print a
      store=store+a
      #print("%s...%s - length %i, strand %i, %i:%i" % (pro[:30], pro[-3:], len(pro), strand, start, end))
  #for x in store:
  if len(store)<2:
    print store
  else:
    p=max(store,key=len)
    print p,len(p)
'''
