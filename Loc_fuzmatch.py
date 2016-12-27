#-------------------------------------------------------------------------------

# For fuzzy matching locations of VC and startup

# Copyright (C) 2015 by Yun Ling

#-------------------------------------------------------------------------------

# fuzzy match two documents
# left: needs full match (internal)
# right: dictionary (external)
# first var: the name for fuzzy matching ("name")
# second var: the time for precise matching (yyyymm)

# example: python "$py\\fuzmatch.py" "ma_left.csv" "ma_right.csv"
# - for ma

# example: pathon "$py\\fuzmatch.py" "vc_left.csv" "vc_right.csv"
# - for vc (startup)

import os, sys, re
# split " ", ",", "!", "."
def splitw(w):
    vec = []
    last = 0
    if (" " not in w) and ("," not in w) and ("!" not in w) and ("." not in w) and ("(" not in w) and (")" not in w):
        return w
    else:
        for i in range(len(w)):
            if (w[i]==" ") or (w[i]==",") or (w[i]=="!") or (w[i]==".") or (w[i]=="(") or (w[i]==")"):
                vec.append(w[last:i])
                last = i+1    
        return vec
os.chdir("D:/Dropbox/Dissertation/Data/CrunchBase/external")
#left = "ma_left.csv"
#right = "ma_right.csv"
left = sys.argv[1]
right = sys.argv[2]
output = sys.argv[3]
print(left)
print(right)
print(output)
with open(left) as f:
    left_list = f.readlines()
with open(right) as f:
    right_list = f.readlines()
left_vars = left_list.pop(0).strip().split(",")
right_vars = right_list.pop(0).strip().split(",")
left_n = map(lambda x: splitw(re.sub('"','',x.split(",")[0]).lower()),left_list)
right_n = map(lambda x: splitw(re.sub('"','',x.split(",")[0]).lower()),right_list)
# wordlist build for the 
# entire corpus 
# translate name list to word list
# for each name
word_all = []
count_all = []
left_w = []
for i in range(len(left_n)):
    # a word
    if isinstance(left_n[i],list)==False:
        if left_n[i] in word_all:
            index = word_all.index(left_n[i])
            count_all[index] = count_all[index]+1
            left_w.append([index])
        else:
            word_all.append(left_n[i])
            count_all.append(1)
            left_w.append([len(word_all)-1])
    else:        
        left_w.append([])
        for l in range(len(left_n[i])):
            if left_n[i][l] in word_all:
                index = word_all.index(left_n[i][l])
                count_all[index] = count_all[index]+1
                left_w[i].append(index)
            else:
                word_all.append(left_n[i][l])
                count_all.append(1)
                left_w[i].append(len(word_all)-1)
right_w = []
for i in range(len(right_n)):
    # a word
    if isinstance(right_n[i],list)==False:
        if right_n[i] in word_all:
            index = word_all.index(right_n[i])
            count_all[index] = count_all[index]+1
            right_w.append([index])
        else:
            word_all.append(right_n[i])
            count_all.append(1)
            right_w.append([len(word_all)-1])
    else:        
        right_w.append([])
        for l in range(len(right_n[i])):
            if right_n[i][l] in word_all:
                index = word_all.index(right_n[i][l])
                count_all[index] = count_all[index]+1
                right_w[i].append(index)
            else:
                word_all.append(right_n[i][l])
                count_all.append(1)
                right_w[i].append(len(word_all)-1)
# for each left, calculate distance 
# and return the most close ones 
match1 = []
match1n = []
for i in range(len(left_w)):
    maxs12 = [0,-1]
    maxs12n = [0,-1]
    s1 = pow(sum(map(lambda x: pow(1.0/count_all[x],2), left_w[i])),0.5)
    for j in range(len(right_w)):
        s2 = pow(sum(map(lambda x: pow(1.0/count_all[x],2), right_w[j])),0.5)
        s12 = sum(map(lambda x: pow(1.0/count_all[x],2)*int(x in right_w[j]), left_w[i]))
        s12n = s12/(s1*s2) ############## ADJUST WEIGHT HERE ############
        if s12 > maxs12[0]:
            maxs12[0] = s12
            maxs12[1] = j
        if s12n > maxs12n[0]:
            maxs12n[0] = s12n
            maxs12n[1] = j
    match1.append(maxs12)
    match1n.append(maxs12n)
# write output file as the "merge"
# command in stata
right_vars1 = map(lambda x: x+"_m1",right_vars)
right_vars2 = map(lambda x: x+"_m2",right_vars)
blank = ",".join([""]*len(right_vars))
with open(output,"w") as f:
    f.write(",".join(left_vars)+","+",".join(right_vars1)+",_merge1,_score1,"+",".join(right_vars2)+",_merge2,_score2\n")
    for i in range(len(left_list)):
        f.write(left_list[i].strip()+",")
        if match1[i][1] == -1:
            f.write(blank+",1,0,")
        else:
            f.write(right_list[match1[i][1]].strip()+",3,"+str(match1[i][0])+",")
        if match1n[i][1] == -1:
            f.write(blank+",1,0\n")
        else:
            f.write(right_list[match1n[i][1]].strip()+",3,"+str(match1n[i][0])+"\n")
            
        