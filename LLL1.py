%reset
import  numpy as np
import numpy.linalg as LA
import math
import random
def dotProduct(a,b):
    n=len(a)
    val=0
    for i in range(n):
        val+=a[i]*b[i]
    #print("val",val)
    return val
def scale(a,k):
    m=[x*k for x in  a]
    #print("m",m)
    return m
def sub(a,b):
    n=len(a)
    #print("sub a",a,b)
    c=list(a)
    for i in range(n):
        c[i]=c[i]-b[i]
    #print("c",c)
    return c
def GramShmid(B):
    ''''coef=[0 for i in range(len(B))]
    basis=B
    n=len(B)
    print(n)
    for i in range(n):
        te=[0]*n
        #print("te",te)
        te[i]=1
        bi=list(B[i])
        #print("i",i)
        #print(bi)'''
    coeff = []
    basis = []
    #print("gbasis",basis,B)
    n = len(B)
    for i in range(n):
        mu = [0] * n
        mu[i] = 1
        bi = list(B[i])
        for j in range(i):
            mu[j] = dotProduct( B[i], basis[j] ) / float(dotProduct(basis[j],basis[j]))
            bi = sub( bi, scale( basis[j], mu[j] ) )
            #print("t",te,"bi",bi,"dot",np.dot(basis[j],basis[j]))
            #print("j",j)
        basis.append(bi)
        #print("basis",bi)
        coeff.append(mu)
        #print(basis[1],coef)
    return basis,coeff
def Lovaz(norm1,norm2,coef,r):
    return (norm2<(r-coef**2)*norm1)
def size_reduce(basis,coeff,idx1,idx2):
    m=round(coeff[idx1][idx2])
    #if(m>1/2):
    basis[idx1]=sub(basis[idx1],scale(basis[idx2],m))
    for j in range(idx2+1):
        coeff[idx1][j]-=m*coeff[idx2][j]
def multiply_list(a):
    i=1
    for x in a:
        i*=x
    return i
def LLL(inbasis,r):
    n=len(inbasis)
    basis=[]
    for i in range(n):
        basis.append(list(inbasis[i]))
    outset_norm=[math.sqrt(dotProduct(i,i)) for i in inbasis]
    vector_all=multiply_list(outset_norm)
    det=abs(LA.det(inbasis))
    HMcoeff=pow(det/vector_all,1/n)
    #print("1",HMcoeff)
    gsbasis,coeff=GramShmid(basis)
    #print(gsbasis)
    #print("\n",coeff)
    gsnorm=[np.dot(x,x) for x in gsbasis]
    #print("gsnorm",gsnorm)
    k=1
    while k<n:
        size_reduce(basis,coeff,k,k-1)
        #print("k",k)
        #print("latter basis",basis,coeff)
        if (Lovaz(gsnorm[k-1],gsnorm[k],coeff[k][k-1],r)):
            mu=coeff[k][k-1]
            no=gsnorm[k]+mu*mu*gsnorm[k-1]
            coeff[k][k-1]=(mu*gsnorm[k-1])/float(no)
            gsnorm[k]*=gsnorm[k-1]/float(no)
            gsnorm[k-1]=no
            basis[k-1],basis[k]=basis[k],basis[k-1]
            for j in range(k-1):
                coeff[k-1][j],coeff[k][j]=coeff[k][j],coeff[k-1][j]
            for j in range(k+1,n):
                coeff[j][k-1],coeff[j][k]=coeff[k][k-1]*coeff[j][k-1]+(1-mu*coeff[k][k-1])*coeff[j][k],coeff[j][k-1]-mu*coeff[j][k]
            k=max(k-1,1)
        else:
            for j in range(k-2,-1,-1):
                size_reduce(basis,coeff,k,j)
            k+=1
    return basis,coeff
def Hdemcoeff(a):
    det=abs(LA.det(a))
    #print("det",det)
    terminal_set_norm=[math.sqrt(dotProduct(i,i)) for i in a]
    #print("ttt",terminal_set_norm)
    vector_all=multiply_list(terminal_set_norm)
    HMcoeff=pow(det/vector_all,1/len(a))
    return HMcoeff
def babai(a,v):
    x=LA.solve(a,v)
    for i in range(len(x)):
        x[i]=round(x[i])
    return x
def Key_creation(a):
    e=np.identity(3,dtype='int64')
    for i in range (20):
        e[(i+1)%3]+=e[i%3]*(random.randint(-3,3))
    w=np.matmul(e,a)
    return w
def GGH_Encry(m,w):
                                                                                                                                                                                                                                                                                                      
    
inbasis=np.array([[19,2,32,46,3,33],[15,42,11,0,3,24],[43,15,0,24,4,16],[20,44,44,0,18,15],[0,48,35,16,31,31],[48,33,32,9,1,29]],dtype='int64')
print("aaa",Hdemcoeff(inbasis))
#mat=np.linspace(0.25,0.99,100)
goodbasis1=np.array([[-97,19,19],[-36,30,86],[-184,-64,78]],dtype='int64')
goodbasis2=np.array([[58,53,-68],[-110,-112,35],[-10,-119,123]],dtype='int64')
print("hdemcoeff good is",Hdemcoeff(goodbasis1))
bad_basis=np.array([[-4179163,-1882253,583183],[-3184353,-1434201,444361],[-5277320,-2376852,736426]],dtype='int64')
bad_basis1=np.array([[324850,-1625176,2734951],[165782,-829409,1395775],[485054,-2426708,4083804]],dtype='int64')
print("hedcoe bad is",Hdemcoeff(bad_basis))
outgoodbasis,outcoeff=LLL(bad_basis1,0.75)
print(outgoodbasis)
print(Hdemcoeff(outgoodbasis))
e=np.array([-79081427,-35617462,11035473],dtype='int64')
e1=np.array([8930810,-44681748,75192665],dtype='int64')
x=babai(np.transpose(outgoodbasis),e1)
x1=babai(np.transpose(bad_basis1),x.dot(outgoodbasis))
print("solve",x,x1)
'''for i in range(100):
    outbasis,outcoeff=LLL(inbasis,mat[i])
    det=abs(LA.det(outbasis))
    terminal_set_norm=[math.sqrt(np.dot(i,i)) for i in outbasis]
    vector_all=multiply_list(terminal_set_norm)
    HMcoeff=pow(det/vector_all,1/len(outbasis))
    print("2",mat[i],HMcoeff)
    '''
