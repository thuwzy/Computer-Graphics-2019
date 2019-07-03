
cntv=0
cntvn=0
cntf=0
v=[]
vt=[]
ff=[]
x=[]
y=[]
z=[]
with open("my_bezier.obj",mode="r") as f:
	lines=f.readlines()
	for line in lines:
		line=line[:-1]
		L=line.split(" ")
		if L[0]=="v":
			v.append((L[1],L[2],L[3]))
			x.append(float(L[1]))
			y.append(float(L[2]))
			z.append(float(L[3]))
		if L[0]=="vt":
			vt.append((L[1],L[2]))
		if L[0]=="f":
			ff.append((L[1].split("//")[0],L[2].split("//")[0],L[3].split("//")[0]))
print(max(x),min(x))
print(max(y),min(y))
print(max(z),min(z))		

with open("bezier2.txt",mode="w+") as f:
	f.write(str(v.__len__())+" "+str(len(vt))+" "+str(len(ff))+"\n")
	for i in v:
		f.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
	for i in vt:
		f.write(str(i[0])+" "+str(i[1])+"\n")
	for i in ff:
		f.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")


