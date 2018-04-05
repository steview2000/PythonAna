#!/u/home/stwe/anaconda/bin/python
#This is only to analyse existing files 

import glob
import os.path
import Image
import string as str
from pylab import *

close("all")
curr_pos= 0

# Name of the DataSet
#DataSet = raw_input("DataSet: ")
#DataSet='BDNF_Phl1024-4'
DataSet='NPY_PHL_DI1024-1'

# The Path of the DataSet - should be changed if you use your external hard drive
#PATH = 'F:\\HOLZLABDATA\ANALYZEDPTIRFData\SecretionBDNF_PHL\BDNF_Phl10232014\\'+DataSet+'/'
PATH = '/u/home/stwe/Desktop/NPY_PHL_DI10242014/'+DataSet+'/'
#PATH='/media/My Passport/HOLZLABDATA/ANALYZEDPTIRFData/SecretionBDNF_PHL/BDNF_Phl10232014/'+DataSet+"/"

#Region of Interest 
ROI = 5

ROI = 2*ROI+1
RW = round((ROI-1.)/2.)


# Read find and read all files, whose filename consists of the following pattern
files1 = glob.glob(PATH+'*HeC*.tif')
files2 = glob.glob(PATH+'*P2S*.tif')
files3 = glob.glob(PATH+'*PS*.tif')

print "size:"

HecAll = -1000*ones((3000,200))
P2SAll = -1000*ones((3000,200))
PSAll = -1000*ones((3000,200))

if size(files1)==0:
  print "\n\t *** Error *** !"
  print '\n \tNo .tif files found under: \n '+PATH+' \n\n '
  quit()

# Sort the list of file names
files1.sort()
files2.sort()
files3.sort()

N=size(files1)

T = zeros((N,5))

M1 = Image.open(files1[0]);
arr = np.asarray(M1.getdata()).reshape(M1.size[1], M1.size[0])

NY=shape(arr)[0]
NX=shape(arr)[1]

MAll1 = zeros((N,NY,NX))
MAll2 = zeros((N,NY,NX))
MAll3 = zeros((N,NY,NX))

print "Loading all files ..."
for l in range(N):
	print "Files loaded: %d/%d"%(l,N)
	M1 = Image.open(files1[l]);
	arr1 = np.asarray(M1.getdata()).reshape(M1.size[1], M1.size[0])
	MAll1[l,:,:] = arr1[:,:]; 

	M2 = Image.open(files2[l]);
	arr2 = np.asarray(M2.getdata()).reshape(M2.size[1], M2.size[0])
	MAll2[l,:,:] = arr2[:,:]; 
	
	M3 = Image.open(files3[l]);
	arr3 = np.asarray(M3.getdata()).reshape(M3.size[1], M3.size[0])
	MAll3[l,:,:] = arr3[:,:]; 

run=1
	
# finding x and y
found = 1;
count=0
xOld=[]
yOld=[]


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
	imshow(MAll1[count,:,:],cmap='gray',interpolation='none')
	gca()
	if y<max(shape(MAll1[0,:,:])):
		gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red'))
	for k in range(size(xOld)):
		gca().add_patch(Rectangle((xOld[k]-RW,yOld[k]-RW),2*RW,2*RW,facecolor='NONE',edgecolor='yellow'))
	title("Hec - %d"%count)
	draw();

def mouseClick(event):
	global x,y
	x,y =  event.xdata,event.ydata
	YN,XN = shape(MAll1[0,:,:])
	if (x<XN)*(y<YN)*(x>=0)*(y>=0):
		figure(1);clf();imshow(MAll1[count,:,:],cmap='gray',interpolation='none')
		subplots_adjust(top=0.95,bottom=0.05,left=0.05,right=1)
		gca()
		gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red'))
		title("frame: %d  x: %d y: %d"%(count,round(x),round(y)))
		for k in range(size(xOld)):
			gca().add_patch(Rectangle((xOld[k]-RW,yOld[k]-RW),2*RW,2*RW,facecolor='NONE',edgecolor='yellow'))
		draw();
	else:
		y=10000


SaveCount=0
HecN_max=0
P2SN_max=0
PSN_max=0

DataDir = PATH+'GreyData/'

files = glob.glob(DataDir+"*.txt")
files.sort()

FileNumber = 1
HecFileName = PATH+"/"+DataSet+"Hec.txt"
P2SFileName = PATH+"/"+DataSet+"P2S.txt"
PSFileName = PATH+"/"+DataSet+"PS.txt"

HalfTimeFileName = PATH+"/"+DataSet+"_HalfTimes.txt"
StartFrameSaveName = PATH+"/"+DataSet+"StartFrame.txt"


while os.path.isfile(HecFileName):
	HecFileName = PATH+"/"+DataSet+"Hec-%d.txt"%FileNumber
	P2SFileName = PATH+"/"+DataSet+"P2S-%d.txt"%FileNumber
	PSFileName = PATH+"/"+DataSet+"PS=%d.txt"%FileNumber
	StartFrameSaveName = PATH+"/"+DataSet+"StartFrame-%d.txt"%FileNumber
	HalfTimeFileName = PATH+"/"+DataSet+"_HalfTimes-%d.txt"%FileNumber
	FileNumber += 1

StartFrameSave = open(StartFrameSaveName,"w")				
StartFrameSave.write("# x\t|y\t|startHec\t|startP2S_PS\n")
StartFrameSave.close()

HalfOut = open(HalfTimeFileName,"w")
HalfOut.write("# x:\t| y\t| d_inc\t| d_dec\t| d_half\n")
HalfOut.close()

for run in range(len(files)):
	SaveFileName = files[run]
	# Extract location
	tempS1 = SaveFileName[len(DataDir+"/"+DataSet+"-x"):] 
	tempS2 = tempS1[:-4]
	x = int(str.split(tempS2,sep="-y=")[0])
	y = int(str.split(tempS2,sep="-y=")[1])

	figure(1);clf();
	subplots_adjust(top=0.95,bottom=0,left=0,right=1)
	imshow(MAll1[count,:,:],cmap='gray',interpolation='none')
	gca()
	if y<max(shape(MAll1[0,:,:])):
		gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red'))

	for k in range(size(xOld)):
		gca().add_patch(Rectangle((xOld[k]-RW,yOld[k]-RW),2*RW,2*RW,facecolor='NONE',edgecolor='yellow'))
	title("Hec - %d"%count)
	cid1 = connect('key_press_event',onclick)	
	cid2 = connect('button_press_event',mouseClick)	
	show(block=False)
	print "\n*****************"	
	print "Choose the point of interest by clicking on Figure 1!\n\n"
	
#	figure(10);clf();
#	subplots_adjust(top=0.95,bottom=0.05,left=0.05,right=1)
#	imshow(MAll1[count,:,:],cmap='gray',interpolation='none')
#	gca()
#	gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red'))
#	title("Click on a reference area!")
#	draw();show(block=False)
#	print "\n*************"
#	print "Choose reference area on figure 10! "
#	p=ginput(n=1,timeout=0)
#	x0,y0=p[0]
#	print "x0: %.1f y0: %.1f"%(x0,y0)	
#	gca().add_patch(Rectangle((x0-RW,y0-RW),2*RW,2*RW,facecolor='NONE',edgecolor='blue'))
#	draw();
#	pause(2)
	#input("Enter to continue")
	
#	SAll = zeros((N,4)) # the actual data
	SAll = genfromtxt(SaveFileName)
		

	figure(5);clf()
	plot(SAll[:,1],'r-')
	plot(SAll[:,2],'b-')
	plot(SAll[:,3],'g-')
	plot(SAll[:,1],'r.')
	plot(SAll[:,2],'b.')
	plot(SAll[:,3],'g.')
	legend(['Hec','P2S(+5000)','PS(+3000)'],loc='right', bbox_to_anchor=(1.1,1),prop={'size':10})
	xlabel("image no")
	ylabel("gray value")
	title("Data Location: x=%d  y=%d"%(x,y))
	draw();show(block=False)

	# Define the directories, where the plots and the data should be stored
	
	# Check whether these directories exist, and create them if they don't.
	if not os.path.exists(DataDir):
	    os.makedirs(DataDir)
	

	a='s'
	if (a=='s') or (a== 'q') or (a=='Q') or (a=='S'):

		xOld.append(x)
		yOld.append(y)
		
		print "\nSaved!\n"
		print "\n*************"
		ana = raw_input("Analyze (y/n)? ")
				
		if (ana=='y') or (ana=='Y'):
			figure(5);#clf()					
			#plot(SAll[:,0],SAll[:,1],'r-')
			#plot(SAll[:,0],SAll[:,1],'r.')
			#title("HeC - grey value")
			#xlabel("frame no")
			#ylabel("grey value")
			#draw();show(block=False);
			print "\n ****************"
			print "Go to Figure 5!"
			while 1:
				try:
					x1 = int(round(input("Choose the starting frame (x-value): ")))	
					break;
				except NameError:
					print "Only numbers alowed"					
		
			y1 = SAll[x1,1]
			print "x1: %d y1: %d"%(x1,y1)
			print ""
			while 1:
				try:
					x_max1 = int(round(input("Choose 1st maximum (x-value): ")))	
					break;
				except NameError:
					print "Only numbers alowed"					

			y_max1 = SAll[x_max1,1]
			print "x_max1: %d y_max1: %d"%(x_max1,y_max1)
			print ""
			while 1:
				try:
					x_max2=int(round(input("Choose 2nd maximum (x-value): ")))	
					break
				except NameError:
					print "Only numbers alowed"					

			y_max2 = SAll[x_max2,1]
			print ""
			print "x_max2: %d y_max2: %d"%(x_max2,y_max2)

			x_half = SAll[:,0]

			half1 = (y1+y_max1)/2.+0*x_half
			half2 = (y1+y_max2)/2.+0*x_half

			figure(5)
			plot(x_half,half2,'-k')
			draw()
			
			while 1:
				try:
					x_half2 = input("x-value for the 2nd intersect (black): ")
					break
				except NameError:
					print "Only numbers alowed"					
			print ""
			
			d_inc = x_max1-x1
			d_dec = x_half2-x_max2
			d_half = x_half2-(x1+x_max1)/2.

			print "d_inc: %.1f\t d_dec: %.1f\t d_half: %.1f"%(d_inc,d_dec,d_half)
			
			cut = max(x1-8,0)
			Hec = (SAll[cut:,1]-SAll[x1,1])/(SAll[x_max1,1]-SAll[x1,1])

			figure(5);clf()
			plot(SAll[:,2],'b-')
			plot(SAll[:,3]+2000,'g-')
			plot(SAll[:,2],'b.')
			plot(SAll[:,3]+2000,'g.')
			legend(['P2S','PS'],loc='right', bbox_to_anchor=(1.1, 1),prop={'size':10})
			draw();

			x4 = int(round(input("Choose the starting frame for PS (x-value): ")))	
			print "x4: %d"%x4
			print "Good choice!!"
			
			cut = max(x4-8,0)
			
			SAll[:,2] = 5000+SAll[:,2]
			SAll[:,3] = 3000+SAll[:,3]
			
			P2S = (SAll[cut:,2]-SAll[x4,2])/SAll[x4,2]
			PS = (SAll[cut:,3]-SAll[x4,3])/SAll[x4,3]
			

			figure(5);clf()
			plot(Hec,'r-')
			plot(P2S,'b-')
			plot(PS,'g-')
			plot(Hec,'r.')
			plot(P2S,'b.')
			plot(PS,'g.')
			legend(['Hec','P2S','PS'],loc='right', bbox_to_anchor=(1.1, 1),prop={'size':10})
			draw();
			
			# Save the dx:
			HalfOut = open(HalfTimeFileName,"a")
			HalfOut.write("%d\t| %d\t| %.2f\t| %.2f\t| %.2f\n"%(x,y,d_inc,d_dec,d_half))
			HalfOut.close()

			HecAll[:(size(Hec)),SaveCount] = Hec
			P2SAll[:(size(P2S)),SaveCount]= P2S
			PSAll[:(size(PS)),SaveCount] = PS 

			HecN_max = max(HecN_max,size(Hec))
			P2SN_max= max(P2SN_max,size(P2S))
			PSN_max= max(PSN_max,size(PS))

			SaveCount +=1
			
			# Save the Hec

			HecFile = open(HecFileName,"w")
			P2SFile = open(P2SFileName,"w")
			PSFile = open(PSFileName,"w")

			StartFrameSave = open(StartFrameSaveName,"a")				
			StartFrameSave.write("%d  |  %d  | %d  |  %d \n"%(x,y,x1,x4))
			StartFrameSave.close()

			HecFile.write("# ")
			P2SFile.write("# ")
			PSFile.write("# ")

			for l1 in range(SaveCount):
				HecFile.write("x=%d y=%d|"%(xOld[l1],yOld[l1]))
				P2SFile.write("x=%d y=%d|"%(xOld[l1],yOld[l1]))
				PSFile.write("x=%d y=%d|"%(xOld[l1],yOld[l1]))

			HecFile.write("\n")
			P2SFile.write("\n")
			PSFile.write("\n")
		
			for l1 in range(HecN_max):
				for l2 in range(SaveCount):
					if HecAll[l1,l2] !=  -1000:
						HecFile.write("%.4f| "%HecAll[l1,l2])
					else:
						HecFile.write("| ")

				HecFile.write("\n")
			for l1 in range(P2SN_max):
				for l2 in range(SaveCount):
					if P2SAll[l1,l2] != -1000:
						P2SFile.write("%.4f| "%P2SAll[l1,l2])
					else:
						P2SFile.write("| ")
				P2SFile.write("\n")
			for l1 in range(PSN_max):
				for l2 in range(SaveCount):
					if PSAll[l1,l2] != -1000:
						PSFile.write("%.4f| "%PSAll[l1,l2])
					else:	
						PSFile.write("| ")
				PSFile.write("\n")
		
			HecFile.close()
			P2SFile.close()
			PSFile.close()

			print "\n*******************"
			a = raw_input("\'Enter\' to continue, \'q\' to quit! ")
			
	elif (a=='q') or (a=='x') or (a=='Q') or (a=='X'):
		print "Stop!"
		break;
	elif a=='':
		print "Continue! \n\n\n\n"


	disconnect(cid1);
	disconnect(cid2)
