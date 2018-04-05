#!/u/home/stwe/anaconda/bin/python
# by Stephan Weiss
# changes addes 

#TODO: - Change input from single tiff file to one big tiff file
#	   - Determine maximum, beginning and end, calculated position where it has 25 and 75% of the maximum
#		 and fit a line. From this get the slope. Also get the integral of the curve from the beginning to the
#		 end.
# Ask Annita for background subtraction
# Ask Annita for the devision by 10000


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
DataSet='BDNF_PH_synWT_1106'
#DataSet='Secretion_TPASynWT_PHL0112'

# The Path of the DataSet - should be changed if you use your external hard drive

#PATH = 'F:\\HOLZLABDATA\ANALYZEDPTIRFData\SecretionBDNF_PHL\BDNF_Phl10232014\\'+DataSet+'/'

PATH = '/home/stevie/ForAnnita/'+DataSet+'/'
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
files2 = glob.glob(PATH+'/*P2S*.tif')
files3 = glob.glob(PATH+'/*PSratioTHR*.tif')

HecAll = -1000*ones((3000,200))
P2SAll = -1000*ones((3000,200))
PSAll = -1000*ones((3000,200))
 
stack=0
if size(files1)==0:
  print "\n\t *** Error *** !" 
  print '\n \tNo .tif files found under: \n '+PATH+' \n\n ' 
  raise SystemExit



run=1

print "Load Cer files. Please be patient ......"
MAll1 = LoadImages(files1)
print "Load P2S files. Please be patient ......"
MAll2 = LoadImages(files2)
print "Load PS files. Please be patient ........"
MAll3 = LoadImages(files3)

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

HecFileName = PATH+"/"+DataSet+"Hec.txt"
P2SFileName = PATH+"/"+DataSet+"P2S.txt"
PSFileName = PATH+"/"+DataSet+"PS.txt"
#SlopeFileName = PATH+"/"+DataSet+"_Slopes.txt"
StartFrameSaveName = PATH+"/"+DataSet+"StartFrame.txt"

ContentName = PATH+"/"+DataSet+"GraphContent.txt"

try:
	C = genfromtxt(ContentName)
	if max(shape(C))>0:
		for i in range(size(C[:,0])):
			FileNumber=FileNumber+1
			xOld.append(C[i,1])
			yOld.append(C[i,2])
except IOError:
	C=[]

FileN = 1
while os.path.isfile(HecFileName):
	HecFileName = PATH+"/"+DataSet+"Hec-%d.txt"%FileN
	P2SFileName = PATH+"/"+DataSet+"P2S-%d.txt"%FileN
	PSFileName = PATH+"/"+DataSet+"PS=%d.txt"%FileN
	StartFrameSaveName = PATH+"/"+DataSet+"StartFrame-%d.txt"%FileN
	FileN +=1

StartFrameSave = open(StartFrameSaveName,"w")				
StartFrameSave.write("#No :\t x\t y\t startHec\t startP2S_PS\n")
StartFrameSave.close()


if os.path.isfile(ContentName):
	Content = open(ContentName,"a")				
else:
	Content = open(ContentName,"w")				
	Content.write("#No :\tx :\ty:\t max:\n")

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
	SAll_O = zeros((N,4)) # the actual data
	BAll_O = zeros((N,4)) # data for the background (to be subtracted)
	
	SAll = zeros((N,4)) # the actual data after calculating <P> and <S> and then calculate <P>/<S> and <P>+2<S>
	BAll = zeros((N,4)) # data for the background (to be subtracted)
	
	print "Calculating ..."
	for l in range(N):
		SAll[l,0] = l
		SAll[l,1] = mean(double(MAll1[l,y-RW:y+RW,x-RW:x+RW]))
		SAll[l,2] = mean(double(MAll2[l,y-RW:y+RW,x-RW:x+RW]))
		SAll[l,3] = mean(double(MAll3[l,y-RW:y+RW,x-RW:x+RW]))
	
		BAll[l,0] = l
		BAll[l,1] = mean(double(MAll1[l,y0-RW:y0+RW,x0-RW:x0+RW]))
		BAll[l,2] = mean(double(MAll2[l,y0-RW:y0+RW,x0-RW:x0+RW]))
		BAll[l,3] = mean(double(MAll3[l,y0-RW:y0+RW,x0-RW:x0+RW]))

	# Plot the Cer data	
	figure(2);clf()
	plot(SAll[:,1],'r.-')
	plot(BAll[:,1],'b.-')
	plot(SAll[:,1]-BAll[:,1],'g.-') # Background substraction
	legend(['Raw ','Ref','Raw-Ref'],loc='right', bbox_to_anchor=(1.1, 1.),prop={'size':10})
	xlabel("frame")
	ylabel("gray value")
	title("Hec- Data Location: x=%d  y=%d"%(x,y)) 
	draw();show(block=False) 
	
	# Plot P2S and PS
	figure(3);clf()
	plot(SAll[:,2],'r.-')
	legend(['P+2S '],loc='right', bbox_to_anchor=(1.1, 1.),prop={'size':10}) 
	xlabel("frame")
	ylabel("gray value")
	title("P2S at: x=%d  y=%d"%(x,y))
	draw();show(block=False)
	
	figure(4);clf()
	plot(SAll[:,3],'r.-')
	#plot(BAll[:,3],'b.-')
	#plot(SAll[:,3]-BAll[:,3],'g.-')
	#plot(SAll_O[:,3],'k.-')
	legend(['P/S'],loc='right', bbox_to_anchor=(1.1, 1.),prop={'size':10}) 
	#xlabel("frame")
	#ylabel("gray value")
	#title("PS new and old at: x=%d  y=%d"%(x,y))
	#draw();show(block=False)

	
	# Background subtraction only for Cer
	SAll[:,1] = SAll[:,1]-BAll[:,1] # background substraction for Cer
#	SAll[:,2] = SAll[:,2]-BAll[:,2] # background subtraction for P+2S
#	SAll[:,3] = SAll[:,3]-BAll[:,3] # background subtraction for P/S

	SAll_O[:,1] = SAll_O[:,1]-BAll_O[:,1]
#	SAll_O[:,2] = SAll_O[:,2]-BAll_O[:,2]
#	SAll_O[:,3] = SAll_O[:,3]-BAll_O[:,3]

	
	# Figure for the Analysis
	figure(5);clf() 
	plot(SAll[:,1],'r.-') 
	plot(SAll[:,2]+5000,'b.-') 
	plot(SAll[:,3],'g.-') 
	legend(['Hec','P2S+5000','PS*3000'],loc='right', bbox_to_anchor=(1.1,1),prop={'size':10}) 
	xlabel("image no") 
	ylabel("gray value") 
	title("Analyse in this figure 5! (x=%d  y=%d)"%(x,y)) 
	draw();show(block=False) 

	
	# These strings contain the filename where the plot and the data are stored
	SaveFileName = DataDir+"/"+DataSet+"-%02d-x=%d-y=%d.txt"%(FileNumber,round(x),round(y))
	SaveFileNameFig = PlotDir+"/"+DataSet+"-x=%02d-y=%d.pdf"%(round(x),round(y))
	SaveFileNameFig2 = PlotDir+"/"+DataSet+"Norm_Cut-x=%02d-y=%d.pdf"%(round(x),round(y))
	
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

	if (a=='s') or (a== 'q') or (a=='Q') or (a=='S'): 
		TxtFile = open(SaveFileName,'w') 
		TxtFile.write("# Location: x=%d\ty=%d Frame: %d\n"%(round(xs),round(ys),count)) 
		TxtFile.write("#ImNumber\tHeC\tP2S\tPS\n")
		#savetxt(SaveFileName,SAll,fmt="%d")
		savetxt(TxtFile,SAll,fmt="%.4lf")
		TxtFile.close() 

		figure(5); 
		savefig(SaveFileNameFig,dpi=200,format='pdf')
		xOld.append(xs)
		yOld.append(ys)
		
		print "\nSaved!\n" 
		print "\n*************" 
		ana = raw_input("Analyze (y/n)? ")
		#ana = 'n' # Remove this for analysis 
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

			x1 = int(round(askVal("Choose the starting frame (x-value): ")))	
			y1 = SAll[x1,1]
			
			print "x1: %d y1: %d"%(x1,y1)
			print ""


			y_max = max(SAll[:,1])
			x_max = int(round(SAll[argmax(SAll[:,1]),0] ))

			# Normalization of Cer/Hec
			SAll[:,1] = (SAll[:,1]-y1)/(y_max-y1)
			print ""


			
			# Ask Annita for the 2nd maximum and whether it is necessary

			cut = max(x1-8,0)
			Hec = (SAll[cut:,1]-SAll[x1,1])/(SAll[x_max,1]-SAll[x1,1])

			figure(5);clf()
			plot(SAll[:,2],'b-')
			plot(SAll[:,3]*2000,'g-')
			plot(SAll[:,2],'b.') 
			plot(SAll[:,3]*2000,'g.') 
			legend(['P2S','PS'],loc='right', bbox_to_anchor=(1.1, 1),prop={'size':10}) 
			draw(); 
 
			x4 = int(round(input("Choose the starting frame for PS (x-value): ")))	 
			print "x4: %d"%x4
			print "Good choice!!"
			
			cut = max(x4-8,0)
			
			SAll[:,2] = SAll[:,2]
			SAll[:,3] = SAll[:,3]
			
			P2S = (SAll[cut:,2]-SAll[x4,2])/abs(SAll[x4,2])
			PS = (SAll[cut:,3]-SAll[x4,3])/abs(SAll[x4,3])
			

			figure(5);clf()
			plot(Hec,'r.-')
			plot(P2S,'b.-')
			plot(PS,'g.-')
			legend(['Hec','P2S','PS'],loc='right', bbox_to_anchor=(1.1, 1),prop={'size':10}) 
			draw();
			savefig(SaveFileNameFig2,dpi=200,format='pdf')
			# Save the dx:
			# Ask Annita what to save

			
			HecAll[:(size(Hec)),SaveCount] = Hec
			P2SAll[:(size(P2S)),SaveCount]= P2S
			PSAll[:(size(PS)),SaveCount] = PS 

			HecN_max = max(HecN_max,size(Hec))
			P2SN_max= max(P2SN_max,size(P2S))
			PSN_max= max(PSN_max,size(PS))

			SaveCount +=1
			

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
			 
		Content.write("%d\t %d\t %d\t %d\n"%(FileNumber,round(xs),round(ys),count))
		if (a=='q') or (a=='Q'):
			FileNumber -=1
		FileNumber +=1
	if (a=='q') or (a=='x') or (a=='Q') or (a=='X'): 
		print "Stop!" 
		break; 
	if a=='': 
		print "Continue! \n\n\n\n" 
 

	disconnect(cid1); 
	disconnect(cid2)

Content.close()

# Ok, let's sort everything in a big file
# load the content file
C = genfromtxt(ContentName)

# create and open the final file:
FOut_Cer = open(PATH+DataSet+"_All_Cer.txt","w")
FOut_P2S = open(PATH+DataSet+"_All_P2S.txt","w")
FOut_PS = open(PATH+DataSet+"_All_PS.txt","w")

# Dimensions of the super file:
SX = 1+shape(C)[0]
SY = N

# Build Super Matrix and Header H
M_Cer = zeros((SY,SX))
M_P2S = zeros((SY,SX))
M_PS = zeros((SY,SX))

H = zeros((4,SX-1))

for l in range(N):
	M_Cer[l,0] = l
	M_P2S[l,0] = l
	M_PS[l,0] = l

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
		M_Cer[l,k] = T[l,1]
		M_P2S[l,k] = T[l,2]
		M_PS[l,k] = T[l,3]

# Write the file:
# First the header:

FOut_Cer.write("ROI:\t")
FOut_P2S.write("ROI:\t")
FOut_PS.write("ROI:\t")
for k in range (SX-1):
	FOut_Cer.write("%d \t"%(1+k))
	FOut_P2S.write("%d \t"%(1+k))
	FOut_PS.write("%d \t"%(1+k))
FOut_Cer.write("\n")
FOut_P2S.write("\n")
FOut_PS.write("\n")

FOut_Cer.write("Frame: \t")
FOut_P2S.write("Frame: \t")
FOut_PS.write("Frame: \t")
for k in range (SX-1):
	FOut_Cer.write("%d \t"%H[2,k])
	FOut_P2S.write("%d \t"%H[2,k])
	FOut_PS.write("%d \t"%H[2,k])
FOut_Cer.write("\n")
FOut_P2S.write("\n")
FOut_PS.write("\n")

FOut_Cer.write("xpos:\t")
FOut_P2S.write("xpos:\t")
FOut_PS.write("xpos:\t")
for k in range (SX-1):
	FOut_Cer.write("%d \t"%H[1,k])
	FOut_P2S.write("%d \t"%H[1,k])
	FOut_PS.write("%d \t"%H[1,k])
FOut_Cer.write("\n")
FOut_P2S.write("\n")
FOut_PS.write("\n")

FOut_Cer.write("ypos:\t")
FOut_P2S.write("ypos:\t")
FOut_PS.write("ypos:\t")
for k in range (SX-1):
	FOut_Cer.write("%d \t"%H[2,k])
	FOut_P2S.write("%d \t"%H[2,k])
	FOut_PS.write("%d \t"%H[2,k])
FOut_Cer.write("\n")
FOut_P2S.write("\n")
FOut_PS.write("\n")

savetxt(FOut_Cer,M_Cer,fmt="%.3f",delimiter=" \t")
savetxt(FOut_P2S,M_P2S,fmt="%.3f",delimiter=" \t")
savetxt(FOut_PS,M_PS,fmt="%.4f",delimiter=" \t")

FOut_Cer.close()
FOut_P2S.close()
FOut_PS.close()
