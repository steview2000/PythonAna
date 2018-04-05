import Image
import numpy as np

# More convenient parameter input
def askVal(s):
	while 1:
		try:
			x = input(s)
			break
		except NameError:
			print "Only numbers alowed"					
	return x

# Function definition for exponential fit
def ExpFit(p,x,y):
	A1 = p[0]
	b1 = p[1]
	x01= p[2]
	A2 = p[3]
	b2 = p[4]
	x02= p[5]
	fitError = (y-A1*np.exp(-b1*(x-x01))-A2*np.exp(-b2*(x-x02)))**2
	return fitError

# Loading all images into one big array
def LoadImages(files):
	stack=0
	if np.size(files)==1:
		#print "All Images in a single file."
		M = Image.open(files[0])
		count = 0
		while True:
			try:
				M.seek(count)
			except EOFError:
				break
			count=count+1
		N=count
		stack=1
	else:
		# Sort the list of file names 
		files.sort()
		N = np.size(files)
	 
#	T = zeros((N,5)) 
	M = Image.open(files[0]); 
	arr = np.asarray(M.getdata()).reshape(M.size[1], M.size[0]) 
	
	del M
	
	NY = np.shape(arr)[0]
	NX = np.shape(arr)[1]
	
	MAll = np.zeros((N,NY,NX),dtype='int16')
	
#	print "Loading all files ..."
	if stack==0:
		for l in range(N):
			M = Image.open(files[l]);
			arr = np.asarray(M.getdata()).reshape(M.size[1], M.size[0])
			MAll[l,:,:] = arr[:,:]; 
		
	else:
		M = Image.open(files[0]);
		for l in range(N):
			M.seek(l)
			arr = np.asarray(M.getdata()).reshape(M.size[1], M.size[0])
			MAll[l,:,:] = arr[:,:]; 
	
	return MAll	

def onclick(event): 
	global count,run 
	global x,y 
	if event.key=='right': 
		count=count+1 
	if event.key=='left': 
		count=count-1 
	if event.key=='up': 
		count=count+5 
	if event.key=='down': 
		count=count-5 
	if count>(shape(MAll1)[0]-1): 
		count=(shape(MAll1)[0]-1) 
	if count<0: 
		count=0 
	figure(1);clf() 
	subplots_adjust(top=0.95,bottom=0,left=0,right=1) 
	imshow(MAll1[count,:,:],cmap='jet') 
	gca() 
	if y<max(shape(MAll1[0,:,:])): 
		gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red')) 
	for k in range(size(xOld)): 
		gca().add_patch(Rectangle((xOld[k]-RW,yOld[k]-RW),2*RW,2*RW,facecolor='NONE',edgecolor='yellow')) 
	title("Ceru - %d"%count) 
	draw(); 
	figure(2);clf();imshow(MAll2[count,:,:],cmap='jet'); subplots_adjust(top=0.95,bottom=0,left=0,right=1) 
	title("P - %d"%count);draw();show(block=False) 
	#figure(3);clf();imshow(MAll3[count,:,:],cmap='gray'); 
	#title("PS - %d"%count);draw();show() 
#	if event.key=='q': 
#		print "disconnect" 
#		figure(1);disconnect(cid1);disconnect(cid2) 
#		run=0 

def mouseClick(event): 
	global x,y 
	x,y =  event.xdata,event.ydata 
	YN,XN = shape(MAll1[0,:,:]) 
	if (x<XN)*(y<YN)*(x>=0)*(y>=0): 
		figure(1);clf();imshow(MAll1[count,:,:],cmap='jet') 
		subplots_adjust(top=0.95,bottom=0.05,left=0.05,right=1) 
		gca() 
		gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red')) 
		title("frame: %d  x: %d y: %d"%(count,round(x),round(y))) 
		for k in range(size(xOld)): 
			gca().add_patch(Rectangle((xOld[k]-RW,yOld[k]-RW),2*RW,2*RW,facecolor='NONE',edgecolor='yellow')) 
		draw(); 
	else: 
		y=10000 

	figure(2);clf();imshow(MAll2[count,:,:],cmap='jet') ;subplots_adjust(top=0.95,bottom=0.05,left=0.05,right=1) 
	gca() 
	gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red')) 
	title("x: %d\t y: %d"%(round(x),round(y))) 
	draw(); 
	#figure(3);clf();imshow(MAll3[count,:,:],cmap='gray') 
	#gca() 
	#gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red')) 
	#title("x: %d\t y: %d"%(round(x),round(y))) 
	#draw(); 

def CPSfromFiles(files):
	MA = LoadImages(files)

	Cer  = MA[0::3,:,:] # Cerul
	P = np.array(MA[1::3,:,:],dtype=float) # S-polarization
	S  = np.array(MA[2::3,:,:],dtype=float) # P-polarisation

	return Cer,P,S
