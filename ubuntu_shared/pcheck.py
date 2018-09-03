import numpy as np

X = np.array([[35, 35.5, 0],
     [0, 0, 0],
     [0,35.5,0],
     [-7,-9,2],
     [-19, 43.5,-10],
     [50.5,37.75,-1.75]])
x = np.array([[124, 182],
     [315, 348],
     [315, 174],
     [362,411],
     [433,100],
     [37,172]])
'''
print("Test values")
X = np.array([[124, 182,1],
     [315*2, 348*2,2],
     [315*3, 174*3,3],
     [362*4,411*4,4],
     [433*5,100*5,5],
     [37*6,172*6,6]])
x = np.array([[124, 182],
     [315, 348],
     [315, 174],
     [362,411],
     [433,100],
     [37,172]])
'''
m,n=X.shape
A = np.zeros([2*m,12])
for i in range(m):
	A[2*i,0:3]=X[i,0:3]	
        A[2*i,3]=1
	A[2*i,8:12]=[-X[i,0]*x[i,0],-X[i,1]*x[i,0],-X[i,2]*x[i,0],-x[i,0]]
	A[2*i+1,4:7]=X[i,0:3]	
        A[2*i+1,7]=1
	A[2*i+1,8:12]=[-X[i,0]*x[i,1],-X[i,1]*x[i,1],-X[i,2]*x[i,1],-x[i,1]]
print("A=",A)

u, s, vh = np.linalg.svd(A,full_matrices=True)
print("s=",s)

print("vh=",vh)
#print(vh[:,11])
P=np.reshape(vh[11,:],(3,4))
print("P=",P)

PN=P/P[0,0]
print("PN=",PN)

print("test")
v1 = vh[11,:]
v2 = np.matmul(A,v1)
#print("v1.shape=",v1.shape)
print("v2=",v2)

a1 = np.array([1])
for i in range(m):
	a2=X[i,:]
	a3=np.concatenate([a2,a1])
	print("a3=",a3)
	xth=np.matmul(P,a3)
	xth=xth[0:2]/xth[2]
	print("xth=",xth)

XT = np.array([[-39.75, 40.5, 0],
     [14, 18, 6.5]])

mt,nt=XT.shape

xt = np.zeros([mt,2])
a1 = np.array([1])
for i in range(mt):
	a2=XT[i,:]
        a3 = np.concatenate([a2,a1])
	xth = np.matmul(P,a3)
	#print("xth=",xth)
	xt[i,:]=xth[0:2]/xth[2]
print("xt=",xt)
