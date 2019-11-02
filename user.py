import LLL as LL
import numpy as np                                                           
import random
import math
def generate_prib(n):                                                            
    a=np.identity(n,dtype='int64')                                               
    r=np.random.randint(-4,4,size=(n,n),dtype='int64')                           
    #print(r)
    basis=4*(math.ceil(math.sqrt(n)))*a+r                                        
    return basis
def generate_pub(basis):
    a=basis
    n=len(a)
    for i in range (2*n):                                                            
        a[(i+1)%n]+=a[i%n]*(random.randint(-2,2))                                
    return a                                                                 
def ggh_encrypt(m,pu,r):
    c=m.dot(pu)                                                                  
    #print("c1",c)                                                                
    c=c+r
    return c
def distant(a,b):
    c=a-b
    print("a,b",a,b,c)
    length=np.linalg.norm(c)
    #print(length)
    return length
n=int(input("输入n"))
pri_basis=generate_prib(n)                                                 
print(pri_basis,LL.Hdemcoeff(pri_basis))                                     
pub_basis=generate_pub(pri_basis)                                            
print(pub_basis,LL.Hdemcoeff(pub_basis))
lll_reduce,coeff=LL.LLL(pub_basis,0.75)
print("lll\n",lll_reduce,LL.Hdemcoeff(lll_reduce))                           
inv_pubasis=np.linalg.inv(pub_basis)
t=np.dot(inv_pubasis,pri_basis)
#print("t\n",t,np.linalg.det(t))
inv_basis=np.linalg.inv(pri_basis)                                           
inv_norm=[math.sqrt(np.dot(i,i)) for i in inv_basis]                         
inv_pubasis=np.linalg.inv(pub_basis)                                         
maxno=max(inv_norm)
#print("inv\n",inv_basis,inv_norm,max(inv_norm))                              
#print("inv_pub",inv_pubasis)                                                 
m=np.random.randint(-n,n,size=(n),dtype='int64')
r=np.random.randint(round(-1/(2*max(inv_norm))),round(1/(2*max(inv_norm))),size=(n))
r2=np.random.randint(-3,3,size=(n))
co=(maxno*math.sqrt(8*math.log(2*120/0.00001)))**-1
print("co",co)
r1=round(1/(2*max(inv_norm)))
#print(r,round(1/(2*max(inv_norm))))
c=ggh_encrypt(m,pub_basis,r2)
print("c",c)
x=LL.babai(np.transpose(lll_reduce),c)
y=LL.babai(np.transpose(pub_basis),c)
x1=LL.babai(np.transpose(pub_basis),x.dot(lll_reduce))
print("m,x1\n",m,x1)
print("y",y)
lengthc_x=distant(c,x1.dot(pub_basis))
lengthc_y=distant(c,y.dot(pub_basis))
print("lengthx,y is \n ",lengthc_x,lengthc_y)
#print(LL.Hdemcoeff(pri_basis))
#print("pri",pri_basis)
#print(LL.Hdemcoeff(inv_basis))
#print("pub\n",pub_basis,LL.Hdemcoeff(pub_basis))
#print("1,2\n",r*(np.linalg.inv(pub_basis)),r*(np.linalg.inv(lll_reduce)))
if m.all()==x1.all():
    print("good")                                                            
else:
    if m.all()==y.all():
         print("bad")
    else:
        print("awful")   
