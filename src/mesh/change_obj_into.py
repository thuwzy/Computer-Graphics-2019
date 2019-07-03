
cntv=0
cntvn=0
cntf=0
v=[]
vn=[]
vt=[]
ff=[]
x=[]
y=[]
z=[]
group=-1
with open("dragon/dragon2.obj",mode="r") as f:
	lines=f.readlines()
	for line in lines:
		line=line[:-1]
		L=[x for x in line.split(' ') if x]
		if len(L)==0:
			continue
		if L[0]=="usemtl":
			group+=1
		if L[0]=="v":
			v.append((L[1],L[2],L[3]))
			x.append(float(L[1]))
			y.append(float(L[2]))
			z.append(float(L[3]))
		if L[0]=="vn":
			vn.append((L[1],L[2],L[3]))
		if L[0]=="vt":
			vt.append((L[1],L[2]))
		if L[0]=="f":
			if len(L)==5:
				'''
				ff.append((group,L[1].split("/")[0],L[2].split("/")[0],L[3].split("/")[0]))
				ff.append((group,L[1].split("/")[-1],L[2].split("/")[-1],L[3].split("/")[-1]))
				ff.append((group,L[1].split("/")[-2],L[2].split("/")[-2],L[3].split("/")[-2]))

				ff.append((group,L[3].split("/")[0],L[4].split("/")[0],L[1].split("/")[0]))				
				ff.append((group,L[3].split("/")[-1],L[4].split("/")[-1],L[1].split("/")[-1]))				
				ff.append((group,L[3].split("/")[-2],L[4].split("/")[-2],L[1].split("/")[-2]))
				'''
				pass
			else:
				
				ff.append((group,L[1].split("/")[0],L[2].split("/")[0],L[3].split("/")[0]))
				ff.append((group,L[1].split("/")[-1],L[2].split("/")[-1],L[3].split("/")[-1]))
				#ff.append((group,L[1].split("/")[-2],L[2].split("/")[-2],L[3].split("/")[-2]))
				
				#ff.append((group,L[1],L[2],L[3]))
print(max(x),min(x))
print(max(y),min(y))
print(max(z),min(z))		
#print(max(vt),min(vt))
with open("dragon/dragon2.txt",mode="w+") as f:
	f.write(str(v.__len__())+" "+str(len(vn))+" "+str(len(vt))+" "+str(len(ff)//2)+"\n")
	for i in v:
		f.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
	for i in vn:
		f.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+"\n")
	for i in vt:
		f.write(str(i[0])+" "+str(i[1])+"\n")
	for i in ff:
		f.write(str(i[0])+" "+str(i[1])+" "+str(i[2])+" "+str(i[3])+"\n")


