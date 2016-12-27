#-------------------------------------------------------------------------------

# For scraping data from CB-API

# Copyright (C) 2014 by Yun Ling

#-------------------------------------------------------------------------------

# libraries

import os, codecs, re, json, time, unirest, math

#-------------------------------------------------------------------------------

# functions

def getdistinct(lista):
    listb = []
    for li in lista:
        if li.__class__==list:
            for l in li:
                if l not in listb:
                    listb.append(l)
        else:
            if li not in listb:
                listb.append(li)
    return listb
    
def fetchitem(data,key):
    if key in data.keys():
        if data[key]==None:
            return ""
        else:
            if data[key].__class__==int:
                ret=str(data[key])
            if data[key].__class__==long:
                ret=str(data[key])
            if data[key].__class__==str:
                ret=data[key]
            if data[key].__class__==unicode:
                ret=data[key]
            if data[key].__class__==str:
                ret=str(data[key])
            if data[key].__class__==bool:
                ret=str(data[key])
            if data[key].__class__==float:
                ret=str(data[key])
            if data[key].__class__==list:
                ret=str(data[key])
            return re.sub('\n+',' ',ret)
    else:
        return ""
        
def getnumber(data,key):
    if key in data.keys():
        if data[key]==None:
            return 0
        else:
            return int(data[key]['paging']['total_items'])
    else:
        return 0
        
def getfirsturl(data,key):
    if key in data.keys():
        return data[key]['paging']['first_page_url'] 
        
def getnexturl(firstpagerawtext):
    try:
        js = json.loads(firstpagerawtext)
        #total_items = js['data']['paging']['total_items']
        #items_per_page = js['data']['paging']['items_per_page']
        #pages = int(math.ceil(total_items*1.0/items_per_page))
        next_page_url = js['data']['paging']['next_page_url']
        return next_page_url
    except:
        return None
    
def unfolddict(dicta):
    if dicta.__class__!=dict:
        return dicta
    else:
        dictb = {}
        for k in dicta.keys():
            if unfolddict(dicta[k]).__class__!=dict:
                dictb[k] = unfolddict(dicta[k])
            else:
                for s in unfolddict(dicta[k]).keys():
                    dictb[k+'/'+s] = unfolddict(dicta[k])[s]
        return dictb  
                
#-------------------------------------------------------------------------------

# api keys

###### privacy.....

keys = [key1,key2,key3,key4,key5,key6,key7,key8,key9,key10,key11,key12,key13, \
    key14,key15,key16,key17,key18,key19,key20,key21,key22,key23,key24]

#-------------------------------------------------------------------------------

# macro variables

base = 'http://api.crunchbase.com/v/2/'
raw_data = 'D:/Dropbox/Dissertation/Databases/crunchbase/raw data'
themes = ['acquisition','category','funding-round','fund-raise','image','ipo','organization','person','product']
theme = themes[0] 
wants = ['acquisition','category','funding-round','fund-raise','ipo','person','product','organization']
print(theme)

################################################################################

# read theme data

os.chdir(raw_data+'/'+theme+'/json')
files = os.listdir('.')
matrix = []
for i in range(len(files)):
    if i%5000==0:
        print(i)
    with codecs.open(files[i],'rb','utf-8') as f:
        txt = f.read()
    try:
        js = json.loads(txt)
        if ('data' in js.keys()) and ('metadata' in js.keys()):
            matrix.append([re.sub('.json','',files[i]),js['data']])
    except:
        pass
keys_properties = map(lambda x: x[1]['properties'].keys(),matrix)
keys_properties = getdistinct(keys_properties)
keys_relationships = map(lambda x: x[1]['relationships'].keys(),matrix)
keys_relationships = getdistinct(keys_relationships)

#-------------------------------------------------------------------------------

# organize properties

os.chdir(raw_data+'/'+theme)
with codecs.open('['+theme+'] '+theme+'_properties.csv','wb','utf-8') as f:
    f.write('uuid,type')
    for p in keys_properties:
        f.write(','+p)
    for p in keys_relationships:
        f.write(','+p)
    f.write('\n')
    for i in range(len(matrix)):
        f.write(matrix[i][1]['uuid']+',')
        f.write(matrix[i][1]['type'])
        for p in keys_properties:
            f.write(',"'+re.sub('"',"'",fetchitem(matrix[i][1]['properties'],p))+'"')
        for p in keys_relationships:
            f.write(',"'+str(getnumber(matrix[i][1]['relationships'],p))+'"')
        f.write('\n')            

#-------------------------------------------------------------------------------

# get long relationships

os.chdir(raw_data+'/'+theme)
if os.path.isdir('./relationships')==False:
    os.mkdir('./relationships')
os.chdir('./relationships')
for p in keys_relationships:
    if os.path.isdir(p)==False:
        os.mkdir(p)
longlist_relationships = []
for p in keys_relationships:
    longlist_relationships.extend(filter(lambda x: getnumber(x[1]['relationships'],p)>8,matrix))
with codecs.open('relationships_long.csv','wb','utf-8') as f:
    for p in keys_relationships:
        for x in longlist_relationships:
            if getnumber(x[1]['relationships'],p)>8:
                f.write(getfirsturl(x[1]['relationships'],p)+'\n')

# download long relationships
del longlist_relationships
with codecs.open('relationships_long.csv','rb','utf-8') as f:
    ll = f.readlines()
ll  = map(lambda x: x.strip(),ll)
print('Download Relationships...'+str(len(ll))+'\n')
for i in range(len(ll)):
    if i%1000==0:
        print(i)
    if os.path.exists(ll[i].split('/')[-1]+'/'+ll[i].split('/')[-2]+'.json')==False:
        time.sleep(1.5)
        unirest.timeout(60)
        try:
            res = unirest.get(ll[i]+'?&user_key='+keys[i%len(keys)])
            with codecs.open(ll[i].split('/')[-1]+'/'+ll[i].split('/')[-2]+'.json','wb','utf-8') as f:
                f.write(res.raw_body)
            nexturl = getnexturl(res.raw_body)
            k = 1
            while nexturl!=None:
                time.sleep(1.5)
                unirest.timeout(60)
                try:
                    res = unirest.get(nexturl+'?&user_key='+keys[(i+k)%len(keys)])
                    with codecs.open(ll[i].split('/')[-1]+'/'+ll[i].split('/')[-2]+'_'+str(k+1)+'.json','wb','utf-8') as f:
                        f.write(res.raw_body)
                    nexturl = getnexturl(res.raw_body)
                    k = k+1
                except:
                    pass
        except:
            pass

# organize relationships
relationships = {}
subkeys_relationships = {}
for p in keys_relationships:
    relationships[p] = []
for p in keys_relationships:
    print(p)
    i = 0
    for x in matrix:
        if i%20000==0:
            print(i)
        if getnumber(x[1]['relationships'],p)>8:
            pages = int(math.ceil(getnumber(x[1]['relationships'],p)/1000.0))
            try:
                with codecs.open(p+'\\'+x[0]+'.json','rb','utf-8') as f:
                    js = json.loads(f.read())
                relationships[p].append([x[1]['uuid'],js['data']['items']])
                for k in range(2,pages+1):
                    with codecs.open(p+'\\'+x[0]+'_'+str(k)+'.json','rb','utf-8') as f:
                        js = json.loads(f.read())
                    relationships[p].append([x[1]['uuid'],js['data']['items']])
            except:
                relationships[p].append([x[1]['uuid'],x[1]['relationships'][p]['items']])
        elif getnumber(x[1]['relationships'],p)>0:
            relationships[p].append([x[1]['uuid'],x[1]['relationships'][p]['items']])
        i = i+1
for p in keys_relationships:
    relationships[p] = map(lambda x: map(lambda y: [x[0],unfolddict(y)],x[1]),relationships[p])
    relationships[p] = reduce(lambda x,y: x+y,relationships[p])
    subkeys_relationships[p] = getdistinct(map(lambda x: x[1].keys(),relationships[p]))
    
os.chdir(raw_data+'/'+theme)
for p in keys_relationships:
    print(p)
    with codecs.open('['+theme+'] '+theme+'_relationships ('+p+').csv','wb','utf-8') as f:
        f.write('uuid')
        for pp in subkeys_relationships[p]:
            f.write(','+pp)
        f.write('\n')
        for x in relationships[p]:
            f.write(x[0])
            for pp in subkeys_relationships[p]:
                f.write(',"'+re.sub('"',"'",fetchitem(x[1],pp))+'"')
            f.write('\n')
            
#-------------------------------------------------------------------------------

# recursively get theme(s) paths
paths = []
for p in keys_relationships:
    print(p)
    for pp in subkeys_relationships[p]:
        for s in wants:
            if ((s+'/') in fetchitem(relationships[p][0][1],pp)):
                paths.extend(map(lambda x: fetchitem(x[1],pp),relationships[p]))
paths_a = []
paths_b = []
i = 0
for x in paths:
    if i%10000==0:
        print(i)
    if len(x)>0:
        if len(x.split('/'))==2:
            paths_a.append(x)
        if (len(x.split('/'))==3) and (x.split('/')[0]=='category'):
            paths_b.append(x)
    i = i+1

# a is non-category
os.chdir(raw_data)
paths_a = filter(lambda x: os.path.exists(x.split('/')[-2]+'/json/'+x.split('/')[-1]+'.json')==False,paths_a)
print(len(paths_a))
print(len(paths_b))
paths_a = getdistinct(paths_a)
paths_a = map(lambda x: x+'\n',paths_a)
os.chdir(raw_data)
with codecs.open("uuid_a.csv","wb","utf-8") as f:
    f.writelines(paths_a)
    
# b is category
sub_categories = getdistinct(map(lambda x: x.split('/')[1],paths_b))
os.chdir(raw_data+'/category/json')
for s in sub_categories:
    if os.path.isdir(s)==False:
        os.mkdir(s)
paths_b = filter(lambda x: os.path.exists(x.split('/')[-2]+'/'+x.split('/')[-1]+'.json')==False,paths_b)
paths_b = getdistinct(paths_b)
paths_b = map(lambda x: x+'\n',paths_b)
os.chdir(raw_data)
with codecs.open("uuid_b.csv","wb","utf-8") as f:
    f.writelines(paths_b)

# house cleaning   
del matrix
del relationships
print('paths_a: '+len(paths_a))
print('paths_b: '+len(paths_b))

# download additional non-category data  
'''  
os.chdir(raw_data)
with codecs.open("uuid_a.csv","rb","utf-8") as f:
    paths = f.readlines()
paths = map(lambda x: x.strip(),paths)
print('Download Uuid...\n')
for i in range(len(paths)):
    if i%1000==0:
        print(i)
    if os.path.exists(paths[i].split('/')[-2]+'/json/'+paths[i].split('/')[-1]+'.json')==False:
        time.sleep(1.5)
        unirest.timeout(60)
        try:
            res = unirest.get(base+paths[i]+'?&user_key='+keys[i%len(keys)])
            with codecs.open(paths[i].split('/')[-2]+'/json/'+paths[i].split('/')[-1]+'.json','wb','utf-8') as f:
                f.write(res.raw_body)
        except:
            pass
'''

################################################################################

# summary
os.chdir(raw_data)
files = {}
for theme in themes:
    files[theme] = os.listdir(theme)
    files[theme] = filter(lambda x: x.endswith('.csv'),files[theme])
lengths = {}
for theme in themes:
    lengths[theme] = []
    for fi in files[theme]:
        with codecs.open(theme+'/'+fi,'rb','utf-8') as f:
            x = f.readlines()
        lengths[theme].append(len(x))
os.chdir(raw_data)
with codecs.open('master.csv','wb','utf-8') as f:
    for theme in themes:
        for i in range(len(files[theme])):
            f.write('"'+files[theme][i]+'",'+str(lengths[theme][i])+'\n')
              