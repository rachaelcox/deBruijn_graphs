file=open("eEF2_full_text_file.txt","r")
import Bio
from Bio.SubsMat import MatrixInfo as matrixFile
matrix=matrixFile.blosum50

P1P2comp_list=[]             #matrix portion
P2P3comp_list=[]
P1P3comp_list=[]
count=0.0                



dic={}                  #dictionary portion
lists=[]                

for lines in file:
    l=lines.rstrip("\r\n")
    line=l.split('\t')
    if dic.has_key(line[0]):
        lists.append(line[1])
        #print lists
        dic[line[0]]=lists
    else:
        lists=[]
        lists.append(line[1])
        dic[line[0]]=lists

new_dict = {}
#for 1 entry in dictionary list
for i in dic:
    val_len=len(dic[i])
    if val_len==1:                                      #1 value for y coordinate key
        new_dict[i] = dic[i]
        #print("All converge {0} {1}".format(i,dic[i]))


#for 2 entries in dictionary list
for i in dic:
    val_len=len(dic[i])
    if val_len==2:
        P1=dic[i][0]
        P2=dic[i][1]
        P1len=len(P1)

        #print(i,P1,P2)
        for x in range(P1len):
            if P1[x]==P2[x]:
                P1P2comp_list.append("same")
            else:
                if P1[x]!=P2[x]:
                    for entry in matrix.keys():
                        if (P1[x],P2[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P1P2comp_list.append("similar")
                            elif matrix[entry]<0:
                                P1P2comp_list.append("different")
                        elif (P2[x],P1[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P1P2comp_list.append("similar")
                            elif matrix[entry]<0:
                                P1P2comp_list.append("different")

        # decide if P2 and P3 are similar enough to add together
        #print(P1P2comp_list)
        count=0.0
        percent_different=0.0
        P1P2=""
        for l in P1P2comp_list:
            if l=="different":
                count=count+1
        #print("diff count=",count)
        #print("kmere length=",P1len)
        percent_different=(count/P1len)*100
                
        #print(percent_different)    
        if percent_different < 16:
            P1P2 = P1+P2
            #print(P1)
            #print(P2)
                    
            new_dict[i]=[P1P2]
        else:
            new_dict[i]=[P1,P2]
                                              
        P1P2comp_list=[]



#for 3 entries in dictionary list
for i in dic:
    length=len(dic[i])
    if length==3:
        P1=dic[i][0]
        P2=dic[i][1]
        P3=dic[i][2]

        # compare P2 and P3, do matrix eval
        P3len=len(P3)
        #print(i,P2,P3)
        for x in range(P3len):                      
            if P2[x]==P3[x]:
                P2P3comp_list.append("same")
            else:
                if P2[x]!=P3[x]:
                    for entry in matrix.keys():
                        if (P2[x],P3[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P2P3comp_list.append("similar")
                            elif matrix[entry]<0:
                                P2P3comp_list.append("different")
                        elif (P3[x],P2[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P2P3comp_list.append("similar")
                            elif matrix[entry]<0:
                                P2P3comp_list.append("different")
            

        # decide if P2 and P3 are similar enough to add together
        #print(P2P3comp_list)
        count=0
        percent_different=0.05
        P2P3=""
        for l in P2P3comp_list:
            if l=="different":
                count=count+1

        #print("diff count=",count)
        #print("kmer len=",P3len)

        percent_different=(count/P3len)*100
        #print(percent_different)       

        if percent_different < 16:
            P2P3=P2+P3
            #print(i,P2P3)
                                                 
        P2P3comp_list=[]
     
        # compare P1 and P2, do matrix eval
        #print(i,P1,P2)
        for x in range(P3len):
            if P1[x]==P2[x]:
                P1P2comp_list.append("same")
            else:
                if P1[x]!=P2[x]:
                    for entry in matrix.keys():
                        if (P1[x],P2[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P1P2comp_list.append("similar")
                            elif matrix[entry]<0:
                                P1P2comp_list.append("different")
                        elif (P2[x],P1[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P1P2comp_list.append("similar")
                            elif matrix[entry]<0:
                                P1P2comp_list.append("different")

        # decide if P2 and P3 are similar enough to add together
        #print(P1P2comp_list)
        count=0.0
        percent_different=0.0
        P1P2=""
        for l in P1P2comp_list:
            if l=="different":
                count=count+1
        #print("diff count=",count)
        #print("kmere length=",P3len)
        percent_different=(count/P3len)*100
                
        #print(percent_different)    
        if percent_different < 16:
            P1P2 = P1+P2    
            #print(i,P1P2)
                                              
        P1P2comp_list=[]



        #compare P1 and P3
        for x in range(P3len):
            if P1[x]==P3[x]:
                P1P2comp_list.append("same")
            else:
                if P1[x]!=P3[x]:
                    for entry in matrix.keys():
                        if (P1[x],P3[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P1P3comp_list.append("similar")
                            elif matrix[entry]<0:
                                P1P3comp_list.append("different")
                        elif (P3[x],P1[x])==(entry[0],entry[1]):
                            #print(matrix[i])
                            if matrix[entry]>=0:
                                P1P3comp_list.append("similar")
                            elif matrix[entry]<0:
                                P1P3comp_list.append("different")

        # decide if P1 and P3 are similar enough to add together
        #print(P1P3comp_list)
        count=0.0
        percent_different=0.0
        P1P3=""
        for l in P1P3comp_list:
            if l=="different":
                count=count+1
        #print("diff count=",count)
        #print("kmere length=",P3len)
        percent_different=(count/P3len)*100
                
        #print(percent_different)    
        if percent_different < 16:
            P1P3 = P1+P3    
            #print(i,P1P3)
                                              
        P1P3comp_list=[]


        if len(P1P2)>0 and len(P2P3)>0 and len(P1P3)>0:
            new_dict[i]=[P1+P2+P3]
        elif len(P1P2)>0:
            new_dict[i]=[P1P2,P3]
        elif len(P2P3)>0:
            new_dict[i]=[P1,P2P3]
        elif len(P1P3)>0:
            new_dict[i]=[P1P3,P2]
        elif len(P1P2)==0 and len(P2P3)==0:
            new_dict[i]=[P1,P2,P3]
        
for i in new_dict:
    #print(i,new_dict[i])
    print(len(new_dict[i]))
    #if len(new_dict[i])==12:
        #print(i,new_dict[i])

#print(new_dict)
