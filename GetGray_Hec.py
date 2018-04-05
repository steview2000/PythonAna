#!/u/home/stwe/anaconda/bin/python
# by Stephan Weiss
# changes addes 

#TODO: - Change input from single tiff file to one big tiff file
#	   - Determine maximum, beginning and end, calculated position where it has 25 and 75% of the maximum
#		 and fit a line. From this get the slope. Also get the integral of the curve from the beginning to the
#		 end.
# Ask Annita for background subtraction
# Todo: Change any StartFram ... to Content

import glob
import os.path
import Image
import os
import sys
from scipy.optimize import leastsq
from pylab import *
from GetGrayLib import *

close("all")
curr_pos= 0


# User parameters :
####################################################################################

# Name of the DataSet
#DataSet = raw_input("DataSet: ")
#DataSet='BDNF_Phl1024-4'
DataSet='Secretion_TPASynWT_PHL0112'

# The Path of the DataSet - should be changed if you use your external hard drive

#PATH = 'F:\\HOLZLABDATA\ANALYZEDPTIRFData\SecretionBDNF_PHL\BDNF_Phl10232014\\'+DataSet+'/'

PATH = '/home/stevie/Desktop/'+DataSet+'/'
#PATH='/media/My Passport/HOLZLABDATA/ANALYZEDPTIRFData/SecretionBDNF_PHL/BDNF_Phl10232014/'+DataSet+"/"


# FrameRate (in hertz)
FR = 27.0

# Define the directories, where the plots and the data should be stored 
PlotDir = PATH+'Plots/'
DataDir = PATH+'GreyData/'

# Check whether these directories exist, and create them if they don't. 
if not os.path.exists(PlotDir):
    os.makedirs(PlotDir)
if not os.path.exists(DataDir):
    os.makedirs(DataDir)

#Region of Interest  
ROI = 5 
ROI = 2*ROI+1 
RW = round((ROI-1.)/2.) 



# Loading all images 
##################################################################################

# Read find and read all files, whose filename consists of the following pattern
files1 = glob.glob(PATH+'/*Cer*.tif')

HecAll = -1000*ones((3000,200))
 
stack=0
if size(files1)==0:
  print "\n\t *** Error *** !" 
  print '\n \tNo .tif files found under: \n '+PATH+' \n\n ' 
  raise SystemExit



run=1

print "Load Cer files. Please be patient ......"
MAll1 = LoadImages(files1)

N = shape(MAll1)[0]
times = 1.0*arange(N)/FR

# Finding x and y
##################################################################################
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
	#figure(2);clf();imshow(MAll2[count,:,:],cmap='gray'); 
	#title("P2S - %d"%count);draw();show() 
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

	#figure(2);clf();imshow(MAll2[count,:,:],cmap='gray') 
	#gca() 
	#gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red')) 
	#title("x: %d\t y: %d"%(round(x),round(y))) 
	#draw(); 
	#figure(3);clf();imshow(MAll3[count,:,:],cmap='gray') 
	#gca() 
	#gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red')) 
	#title("x: %d\t y: %d"%(round(x),round(y))) 
	#draw(); 


SaveCount=0
HecN_max=0
P2SN_max=0
PSN_max=0

# Save the Hec
FileNumber = 1

StartFrameSaveName = PATH+"/"+DataSet+"GraphContent.txt"

try:
	C = genfromtxt(StartFrameSaveName)
	for i in range(size(C[:,0])):
		FileNumber=FileNumber+1
		xOld.append(C[i,1])
		yOld.append(C[i,2])
except IOError:
	C=[]


if os.path.isfile(StartFrameSaveName):
	StartFrame = open(StartFrameSaveName,"a")				
else:
	StartFrame = open(StartFrameSaveName,"w")				
	StartFrame.write("#No :\tx :\ty:\t max:\n")

# Main loop
###########################################################

while run==1:
	#figure(2);clf();imshow(MAll2[count,:,:],cmap='gray');
	#figure(3);clf();imshow(MAll3[count,:,:],cmap='gray');
	y=100000
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
	while y>max(shape(MAll1[0,:,:])): 
		pause(1) 
	print "got it!" 
	
	figure(10);clf(); 
	subplots_adjust(top=0.95,bottom=0.05,left=0.05,right=1) 
	imshow(MAll1[count,:,:],cmap='gray',interpolation='none') 
	gca() 
	gca().add_patch(Rectangle((x-RW,y-RW),2*RW,2*RW,facecolor='NONE',edgecolor='red')) 
	title("Click on a reference area!") 
	draw();show(block=False) 
	print "\n*************" 
	print "Choose reference area on figure 10! " 
	p=ginput(n=1,timeout=0) 
	x0,y0=p[0] 
	print "x0: %.1f y0: %.1f"%(x0,y0)	
	gca().add_patch(Rectangle((x0-RW,y0-RW),2*RW,2*RW,facecolor='NONE',edgecolor='blue'))
	draw();
	pause(2)
	#input("Enter to continue")

	# The actual analysis start
	##############################################################################
	SAll = zeros((N,2)) # the actual data
	BAll = zeros((N,2)) # data for the background (to be subtracted)
	
	print "Calculating ..."
	for l in range(N):
		SAll[l,0] = l
		SAll[l,1] = mean(double(MAll1[l,y-RW:y+RW,x-RW:x+RW]))
	
		BAll[l,0] = l
		BAll[l,1] = mean(double(MAll1[l,y0-RW:y0+RW,x0-RW:x0+RW]))
	
	# Plot the Cer data	
	figure(2);clf()
	plot(times,SAll[:,1],'r.-')
	plot(times,BAll[:,1],'b.-')
	plot(times,SAll[:,1]-BAll[:,1],'g.-') # Background substraction
	legend(['Raw ','Ref','Raw-Ref'],loc='right', bbox_to_anchor=(1.1, 1.),prop={'size':10})
	xlabel("time (s)")
	ylabel("gray value")
	title("Hec- Data Location: x=%d  y=%d"%(round(x),round(y))) 
	draw();show(block=False) 
	

	# Background subtraction only for Cer
	SAll[:,1]=SAll[:,1]-BAll[:,1]
	
	# Figure for the Analysis
	figure(5);clf() 
	plot(SAll[:,1],'r.-') 
#	legend(['Hec','P2S(+5000)','PS(+3000)'],loc='right', bbox_to_anchor=(1.1,1),prop={'size':10}) 
	xlabel("image no") 
	ylabel("gray value") 
	title("Analyse in this figure 5! (x=%d  y=%d)"%(x,y)) 
	draw();show(block=False) 

	
	# These strings contain the filename where the plot and the data are stored
	SaveFileName = DataDir+"/"+DataSet+"-%02d-x=%d-y=%d.txt"%(FileNumber,round(x),round(y))
	SaveFileNameFig = PlotDir+"/"+DataSet+"-%02d-x=%d-y=%d.pdf"%(FileNumber,round(x),round(y))
	# x,y should not be changed anymore by an accidental click
	xs = x
	ys = y
	wrongInput = 1
	while(wrongInput):
		print "\n******************"	
		print "What do you want to do next?" 
		print "\t \'Enter\' - Continue without saving"	 
		print "\t \'r\'\t - Change Size of region of interest and Continue" 
		print "\t \'s\'\t - Save and Continue" 
		print "\t \'q\'\t - Save and Quite" 
		print "\t \'x\'\t - Quite without saving" 
		a=raw_input("Choose: ")
		if (a=='') or (a=='r') or (a=='R') or (a=='s') or (a=='q') or (a=='x') or (a=='S') or (a=='Q') or (a=='X'): 
			wrongInput=0

	if (a=='r') or (a=='R'): 
		print "" 
		print "Currently it is R=%d. (The ROI is a square with side length 2R+1)"%RW  
		Rnew = raw_input("New R (default is R=%d): "%RW) 
		if Rnew=='':
			RW=RW
		else: 
			RW=int(Rnew) 
		print "R=%d"%RW 

	elif (a=='s') or (a== 'q') or (a=='Q') or (a=='S'): 
		TxtFile = open(SaveFileName,'w') 
		TxtFile.write("# Location: x=%d\ty=%d Frame: %d\n"%(round(xs),round(ys),count)) 
		TxtFile.write("#ImNumber\tHeC\n")
		#savetxt(SaveFileName,SAll,fmt="%d")
		savetxt(TxtFile,SAll,fmt="%.2lf")
		TxtFile.close() 

		figure(5); 
		savefig(SaveFileNameFig,dpi=200,format='pdf')
		xOld.append(xs)
		yOld.append(ys)
		
		print "\nSaved!\n" 
		print "\n*************" 

		start = max(0,count-20)
		end = min(N,count+20)
		x_max = start+argmax(SAll[start:end,1])

		StartFrame.write("%d\t %d\t %d\t %d\n"%(FileNumber,round(xs),round(ys),count))
		FileNumber=FileNumber+1;
		if (a=='q') or (a=='Q'):		 
			print "Stop!"
			break;
	elif (a=='q') or (a=='x') or (a=='Q') or (a=='X'): 
		print "Stop!" 
		break; 
	elif a=='': 
		print "Continue! \n\n\n\n" 
 
	disconnect(cid1); 
	disconnect(cid2) 
StartFrame.close() 

# Ok, let's sort everything in a big file
# load the content file
C = genfromtxt(StartFrameSaveName)

# create and open the final file:
FOut = open(PATH+DataSet+"_All.txt","w")


# Dimensions of the super file:
SX = 1+shape(C)[0]
SY = N

# Build Super Matrix and Header H
M = zeros((SY,SX))
H = zeros((4,SX-1))

for l in range(N):
	M[l,0] = l

for k in range(1,SX):
	x = C[k-1,1]
	y = C[k-1,2]
	print "x: %d"%x
	SaveFileName = DataDir+"/"+DataSet+"-%02d-x=%d-y=%d.txt"%(k,x,y)
	T = genfromtxt(SaveFileName)
	H[0,k-1] = k 
	H[1,k-1] = x
	H[2,k-1] = y
	H[3,k-1] = C[k-1,3]
	for l in range(SY):
		M[l,k] = T[l,1]
	
# Write the file:
# First the header:

FOut.write("ROI:\t")
for k in range (SX-1):
	FOut.write("%.4f \t"%(1+k))
FOut.write("\n")

FOut.write("Frame: \t")
for k in range (SX-1):
	FOut.write("%.4f \t"%H[2,k])
FOut.write("\n")

FOut.write("xpos:\t")
for k in range (SX-1):
	FOut.write("%.4f \t"%H[1,k])
FOut.write("\n")

FOut.write("ypos:\t")
for k in range (SX-1):
	FOut.write("%.4f \t"%H[2,k])
FOut.write("\n")

savetxt(FOut,M,fmt="%.3f",delimiter=" \t")
FOut.close()
