import sys,getopt,os,commands,SVCNV_set
#parameter setting
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"i:r:d:c:p:t:o:f:g:")
inFile = ""
cutoff = 0
percent = 0.8
distance = 100
tabixPath = ".../bin/tabix-0.2.4/tabix"
for op, value in opts:
	if op == "-i":
	    inFile = value
	if op == "-d":
	    distance = value
	if op == "-r":
	    allele_freq = value
	if op == "-c":
	    cutoff = float(value)
	if op == "-p":
	    percent = float(value)
	if op == "-t":
	    template_database = value
	if op == "-o":
	    filter_pass = value
	if op == "-f":
	    filter_fail = value
	if op == "-g":
	    filter_gap = value         
if inFile == "":
	print("-i invalid")
	sys.exit()
      
fp=open(filter_pass,'w')
ff=open(filter_fail,'w')
fg=open(filter_gap,'w')
# read sv_list from analysis file
with open(inFile,'r') as f:
	for line in f:
	    if line.startswith('#') or line[:3]!="chr":
#		print(line.strip())
		continue
            item = line.strip().split("\t")
            chr = item[0]    #chrom
            start = min(int(item[1]),int(item[2]))     #start_pos        
            end = max(int(item[1]),int(item[2]))     #end_pos
            svtype=item[3]   
            start1 = int(start-(end-start)*percent)
            if start1<0:
                start1=0
            end1 = int(end+(end-start)*percent)
            start2 = start-1000
            start3 = start+1000
            if start2<0:
                start2=0            
            end2 = end+500
            sub_len=int(start-end)
            svcnv1=SVCNV_set.SVCNV(line)
            sml=[]
            svcnv_list=[]
# read sv_list from database
            if svtype=="DEL" or svtype=="DUP": 
                svcnv_list = commands.getoutput(tabixPath + " -f " + template_database +" "+ chr + ":" + str(start) + "-" + str(end) +"|grep "+svtype).split("\n")
            else:
                svcnv_list = commands.getoutput(tabixPath + " -f " + template_database +" "+ chr + ":" + str(start2) + "-" + str(start3)+"|grep "+svtype).split("\n")
            if len(svcnv_list)==1 and svcnv_list[0]=="":
                fp.write(line.strip()+"\tNot_in_database\n")
                continue
            else:
                sv_dict=[]
                for svcnvs in svcnv_list:
                    result = svcnvs.strip().split("\t")
                    if result[4]>=allele_freq or template_database==".../trio_cnv/GRCh37_hg19_variants_2016-05-15.txt":
                        s = SVCNV_set.SVCNV(svcnvs)
                        continue
                    sv_dict.append(s)
# simplify sv_list from database
                if svtype=="DEL" or svtype=="DUP": 
                    sm1=SVCNV_set.simplify_by_overlap(sv_dict)
                else:
                    sm1=SVCNV_set.simplify_by_breakpoint(sv_dict)
# calculate svcnv exclude ratio
                sm2=[]
                sv_dict_ori=[]
                sv_dict_sim=[]
                sv_dict_gap=[]
                if svtype=="DEL" or svtype=="DUP":
                    sm2= SVCNV_set.subtract_by_overlap(sm1,svcnv1,percent)
                else:
                    sm2 = SVCNV_set.subtract_by_breakpoint(sm1,svcnv1,distance)
                subtract_list_len=0
                svcnv1=SVCNV_set.SVCNV(line)
                for smo in sv_dict:
                    sv_dict_ori.append(smo.chr + ":" + str(smo.start_pos) + ":" + str(smo.end_pos) + ":" + str(smo.length) + ":" + smo.svcnv_type  + ":" + smo.info)                      
                for sms in sm1:               
                    sv_dict_sim.append(sms.chr + ":" + str(sms.start_pos) + ":" + str(sms.end_pos) + ":" + str(sms.length) + ":" + sms.svcnv_type  + ":" + sms.info)  
                for sm in sm2:
                    subtract_list_len+=sm.length 
                    sv_dict_gap.append(sm.chr + ":" + str(sm.start_pos) + ":" + str(sm.end_pos) + ":" + str(sm.length) + ":" + sm.svcnv_type  + ":" + sm.info)  
                    fg.write(sm.chr + "\t" + str(sm.start_pos) + "\t" + str(sm.end_pos) + "\t" + str(sm.length) + "\t" + sm.svcnv_type + "\t" + sm.info + "\n" )
                if svcnv1.length==0:    
                    svcnv_ratio=len(filter(None,sm2))
                else:                            
                    svcnv_ratio=float(1-subtract_list_len/svcnv1.length)                
                if svcnv_ratio<percent:
                    fp.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + str(svcnv1.length) + "\t" + svcnv1.svcnv_type + "\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio)+ "\n" )
                    ff.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + str(svcnv1.length) + "\t" + svcnv1.svcnv_type + "\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio)+ "\t" + str(subtract_list_len)+"\t" +";".join(sv_dict_sim)+ "\n" )
                if svcnv_ratio>=percent:
                    ff.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + str(svcnv1.length) + "\t" + svcnv1.svcnv_type + "\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio)+ "\t" + str(subtract_list_len)+"\t" +";".join(sv_dict_sim)+ "\n" )
