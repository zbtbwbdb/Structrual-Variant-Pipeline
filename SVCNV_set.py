#chr: chromsome
#start_pos: start position
#end_pos: end position
#svcnv_type: svcnv type(DUP,DEL,INV,INS...)

class SVCNV:
    def __init__(self,cnv_line):
        cols = cnv_line.strip().split('\t')
        self.info=";".join(cols[4:])
        self.chr = cols[0]
        self.start_pos = int(min(int(cols[1]),int(cols[2])))
        self.end_pos = int(max(int(cols[1]),int(cols[2])))
        self.svcnv_type = cols[3][:3]
        if self.svcnv_type =="DEL" or self.svcnv_type =="DUP": 
            self.length=abs(self.end_pos-self.start_pos)
        else:
            self.length=0    
                
    def __lt__(self, other):
        if self.chr < other.chr:
            return True
        elif self.chr == other.chr:
            if self.start_pos < other.start_pos:
                return True
            else:
                return False
        else:
            return False    

            
class SVCNV_set:
    def __init__(self,cnv_line):
        cols = cnv_line.split('\t')
        self.info=";".join(cols[4:])
        self.chr = cols[0]
        self.start_pos = int(min(int(cols[1]),int(cols[2])))
        self.end_pos =int(max(int(cols[1]),int(cols[2])))
        self.svcnv_type = cols[3][:3]
        self.length=abs(self.end_pos-self.start_pos)
                
#class SVCNV_sort:

#Function DELDUP length

# SVCNV shared

# SVCNV different
# SVCNV different
class SVCNV_transform:
    def __init__(self,svcnv):
        self.chr = svcnv.chr
        self.info=svcnv.info
        self.start_pos = svcnv.start_pos
        self.end_pos = svcnv.end_pos
        self.svcnv_type = svcnv.svcnv_type
        self.transform_list = [svcnv]
        if self.svcnv_type =="DEL" or self.svcnv_type =="DUP": 
            self.length=abs(self.end_pos-self.start_pos)
        else:
            self.length=0   

def subtract_by_overlap(svcnv_list,svcnv,percent):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sub_list = []
    SVCNV_transform_current = None
    sim_len=len(svcnv_list)
    for idx,sim_current in enumerate(svcnv_list_sorted):
        if sim_current.chr==svcnv.chr and (sim_current.start_pos>=svcnv.end_pos or sim_current.end_pos<=svcnv.start_pos):
            continue
        elif sim_current.chr==svcnv.chr and sim_current.start_pos<svcnv.start_pos and sim_current.end_pos>svcnv.start_pos and sim_current.end_pos<svcnv.end_pos:
            svcnv.start_pos=sim_current.end_pos
        elif sim_current.chr==svcnv.chr and sim_current.start_pos<svcnv.start_pos and sim_current.end_pos>=svcnv.end_pos:
            sim_current.start_pos=svcnv.end_pos
            sim_current.end_pos=svcnv.end_pos
            sim_current.length=abs(sim_current.end_pos-sim_current.start_pos)
            sub_list.append(sim_current)
            break         
        elif sim_current.chr==svcnv.chr and sim_current.start_pos>svcnv.start_pos and sim_current.end_pos<svcnv.end_pos:
            SVCNV_transform_current=SVCNV_transform(svcnv)
            SVCNV_transform_current.end_pos =sim_current.start_pos
            SVCNV_transform_current.length=abs(SVCNV_transform_current.end_pos-SVCNV_transform_current.start_pos)
            SVCNV_transform_current.info=sim_current.info
            sub_list.append(SVCNV_transform_current)
            svcnv.start_pos=sim_current.end_pos
        elif sim_current.chr==svcnv.chr and sim_current.start_pos>svcnv.start_pos and sim_current.start_pos<=svcnv.end_pos and sim_current.end_pos>=svcnv.end_pos:
            SVCNV_transform_current=SVCNV_transform(svcnv)
            SVCNV_transform_current.end_pos =sim_current.start_pos
            SVCNV_transform_current.length=abs(SVCNV_transform_current.end_pos-SVCNV_transform_current.start_pos)   
            SVCNV_transform_current.info=sim_current.info         
            sub_list.append(SVCNV_transform_current)                
            svcnv.start_pos=svcnv.end_pos
            break                
        if svcnv.chr != sim_current.chr:
            continue
    return sub_list

def subtract_by_breakpoint(svcnv_list,svcnv,distance):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sub_list = []
    SVCNV_transform_current = None
    for idx,sim_current in enumerate(svcnv_list_sorted):
        if  sim_current.chr==svcnv.chr and abs(sim_current.start_pos-svcnv.start_pos)<distance and abs(sim_current.end_pos-svcnv.end_pos)<distance:
            SVCNV_transform_current = SVCNV_transform(sim_current)
            sub_list.append(SVCNV_transform_current)
        if  sim_current.chr==svcnv.chr and abs(sim_current.start_pos-svcnv.start_pos)>=distance and abs(sim_current.end_pos-svcnv.end_pos)>=distance:
            continue
        if svcnv.chr != sim_current.chr:
            continue
    return sub_list

# SVCNV union
#Function SVCNV simplify: 
def merge_by_overlap(svcnv_list,percent):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sm_list = []
    SVCNV_transform_current = None
    for idx,svcnv_current in enumerate(svcnv_list_sorted):
        if idx == 0:
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
            continue
        if SVCNV_transform_current.chr == svcnv_current.chr and svcnv_current.start_pos <= SVCNV_transform_current.end_pos:
            overlap_ratio = float(min(svcnv_current.end_pos,SVCNV_transform_current.end_pos) - svcnv_current.start_pos) / float(max(svcnv_current.end_pos,SVCNV_transform_current.end_pos) - SVCNV_transform_current.start_pos)
            if overlap_ratio >= percent:
                SVCNV_transform_current.end_pos = max(svcnv_current.end_pos,SVCNV_transform_current.end_pos)
                SVCNV_transform_current.transform_list.append(svcnv_current)
            else:
                sm_list.append(SVCNV_transform_current)
                SVCNV_transform_current = SVCNV_transform(svcnv_current)
        else:
            sm_list.append(SVCNV_transform_current)
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
    sm_list.append(SVCNV_transform_current)
    return sm_list

def merge_by_breakpoint(svcnv_list,distance):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sm_list = []
    SVCNV_transform_current = None
    for idx,svcnv_current in enumerate(svcnv_list_sorted):
        if idx == 0:
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
            continue
        if SVCNV_transform_current.chr == svcnv_current.chr and (abs(svcnv_current.start_pos - SVCNV_transform_current.start_pos) <= distance or abs(svcnv_current.end_pos - SVCNV_transform_current.end_pos) <= distance):
            SVCNV_transform_current.end_pos = max(svcnv_current.end_pos,SVCNV_transform_current.end_pos)
            SVCNV_transform_current.transform_list.append(svcnv_current)
        else:
            sm_list.append(SVCNV_transform_current)
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
    sm_list.append(SVCNV_transform_current)
    return sm_list

#Function SVCNV simplify:                  
def simplify_by_overlap(svcnv_list):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sm_list = []
    SVCNV_transform_current = None
    for idx,svcnv_current in enumerate(svcnv_list_sorted):
        if idx == 0:
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
            continue
        if svcnv_current.start_pos >= SVCNV_transform_current.start_pos and svcnv_current.start_pos <= SVCNV_transform_current.end_pos :
                SVCNV_transform_current.end_pos = max(svcnv_current.end_pos,SVCNV_transform_current.end_pos)
                SVCNV_transform_current.length=abs(SVCNV_transform_current.end_pos-SVCNV_transform_current.start_pos)
        else:
            sm_list.append(SVCNV_transform_current)
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
    sm_list.append(SVCNV_transform_current)
    return sm_list
    
def simplify_by_breakpoint(svcnv_list):
    if len(svcnv_list) == 0:
        return []
    svcnv_list_sorted = sorted(svcnv_list)
    sm_list = []
    SVCNV_transform_current = None
    for idx,svcnv_current in enumerate(svcnv_list_sorted):
        if idx == 0:
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
            continue
        if abs(svcnv_current.start_pos-SVCNV_transform_current.start_pos)<100 and abs(svcnv_current.end_pos-SVCNV_transform_current.end_pos)<100 :
                SVCNV_transform_current.start_pos = int((svcnv_current.start_pos+SVCNV_transform_current.start_pos)/2)
                SVCNV_transform_current.end_pos = int((svcnv_current.end_pos+SVCNV_transform_current.end_pos)/2)
        else:
            sm_list.append(SVCNV_transform_current)
            SVCNV_transform_current = SVCNV_transform(svcnv_current)
    sm_list.append(SVCNV_transform_current)
    return sm_list
        
