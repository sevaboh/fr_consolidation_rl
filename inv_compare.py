#!/usr/bin/python
import sys
import math
if len(sys.argv)!=3:
    quit()
n1=sys.argv[1]
n2=sys.argv[2]
print(n1)
print(n2)
f1=open(n1,"rt")
f2=open(n2,"rt")
ls1=f1.readlines()
ls2=f2.readlines()
f1.close()
f2.close()
for i in range(len(ls1)):
    ls1[i]=ls1[i].split(" ")
for i in range(len(ls2)):
    ls2[i]=ls2[i].split(" ")
err=0
nv=0
last_t=0
for i in range(len(ls1)):
    t1=float(ls1[i][1])
    for j in range(len(ls2)):
        t2=float(ls2[j][1])
        if abs(t2-t1)<0.01:
            for k in range(len(ls1[i])):
                if k>=14 and k<len(ls2[j]) and ls1[i][k]!="\n" and ls2[j][k]!="\n":
                    err=err+(float(ls2[j][k])-float(ls1[i][k]))*(float(ls2[j][k])-float(ls1[i][k]))
                    nv=nv+1
    if last_t<1 and t1>=1:
        print(str(t1)+" "+str(err)+" "+str(nv)+" "+str(math.sqrt(err/nv)))
    last_t=t1
print(str(t1)+" "+str(err)+" "+str(nv)+" "+str(math.sqrt(err/nv)))
