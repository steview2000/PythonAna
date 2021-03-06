#!/u/home/stwe/anaconda/bin/python
# by Stephan Weiss
# changes addes 

#TODO: - Change input from single tiff file to one big tiff file
#	   - Determine maximum, beginning and end, calculated position where it has 25 and 75% of the maximum
#		 and fit a line. From this get the slope. Also get the integral of the curve from the beginning to the
#		 end.
# Ask Annita for background subtraction

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

def ExpFit(p,x,y):
	A1 = p[0]
	b1 = p[1]
	x01= p[2]
	A2 = p[0]
	b2 = p[1]
	x02= p[2]
	fitError = (y-A1*np.exp(-b1*(x-x01))-A2*np.exp(-b2*(x-x02)))**2
	return fitError

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
#MAll2 = LoadImages(files2)
#print "Load PS files. Please be patient ........"
#MAll3 = LoadImages(files3)

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
HalfTimeFileName = PATH+"/"+DataSet+"_HalfTimes.txt"
SlopeFileName = PATH+"/"+DataSet+"_Slopes.txt"
IntFileName = PATH+"/"+DataSet+"_Integral.txt"
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
	SAll = zeros((N,4)) # the actual data
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
	plot(times,SAll[:,1],'r.-')
	plot(times,BAll[:,1],'b.-')
	plot(times,SAll[:,1]-BAll[:,1],'g.-') # Background substraction
	legend(['Raw ','Ref','Raw-Ref'],loc='right', bbox_to_anchor=(1.1, 1.),prop={'size':10})
	xlabel("time (s)")
	ylabel("gray value")
	title("Hec- Data Location: x=%d  y=%d"%(x,y)) 
	draw();show(block=False) 
	
	# Plot P2S and PS
	figure(3);clf()
	plot(times,SAll[:,2],'r.-')
	plot(times,SAll[:,3],'b.-')
	legend(['P2S','PS'],loc='right', bbox_to_anchor=(1.1, 1.),prop={'size':10}) 
	xlabel("time (s)")
	ylabel("gray value")
	title("P2S and PS- at: x=%d  y=%d"%(x,y))
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
	SaveFileName = DataDir+"/"+DataSet+"-x=%d-y=%d.txt"%(round(x),round(y))
	SaveFileNameFig = PlotDir+"/"+DataSet+"-x=%d-y=%d.pdf"%(round(x),round(y))

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
		TxtFile.write("# Location: x=%d\ty=%d Frame: %d\n"%(round(x),round(y),count)) 
		TxtFile.write("#ImNumber\tHeC\tP2S\tPS\n")
		#savetxt(SaveFileName,SAll,fmt="%d")
		savetxt(TxtFile,SAll,fmt="%d")
		TxtFile.close() 

		figure(5); 
		savefig(SaveFileNameFig,dpi=200,format='pdf')
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

			x1 = int(round(askVal("Choose the starting frame (x-value): ")))	
			y1 = SAll[x1,1]
			
			print "x1: %d y1: %d"%(x1,y1)
			print ""

			x_end = int(round(askVal("Choose the last frame (x-value): ")))	
			

			y_max = max(SAll[:,1])
			x_max = int(round(SAll[argmax(SAll[:,1]),0] ))

			# Normalization of Cer/Hec
			SAll[:,1] = (SAll[:,1]-y1)/(y_max-y1)
			y_end = SAll[x_end,1]
			print "x_end: %d y_end: %d"%(x_end,y_end)
			print ""
			y_max=1

			y_25 = 0.25
			y_75 = 0.75
			for l in range(x1,x_max):
				if (SAll[l,1]<y_25) and	(SAll[l+1,1]>=y_25):
					x_25 = l + (y_25-SAll[l,1])/(SAll[l+1,1]-SAll[l,1])
				if (SAll[l,1]<y_75) and	(SAll[l+1,1]>=y_75):
					x_75 = l + (y_75-SAll[l,1])/(SAll[l+1,1]-SAll[l,1])

			slope = (y_75-y_25)/(FR*(x_75-x_25))
			print "slope: %lf"%slope	

			
			# Ask Annita for the 2nd maximum and whether it is necessary

			#x_max1 = int(round(askVal("Choose 1st maximum (x-value): ")))	
			#y_max1 = SAll[x_max1,1]
			#print "x_max1: %d y_max1: %d"%(x_max1,y_max1)
			#print ""
			
			x_max2 = int(round(askVal("Choose 2nd maximum (x-value): ")))
			y_max2 = SAll[x_max2,1]
			print ""
			print "x_max2: %d y_max2: %d"%(x_max2,y_max2)
		
			# Exponential fit
			p0=zeros(3)
			p0 = 0.015,2.5,times[x_max2]#,0.1,0.01,times[x_max2]+6
			pls = leastsq(ExpFit,p0,args=(times[x_max2:x_end],SAll[x_max2:x_end,1]))
			A1,b1,x01,A2,b2,x02 = pls[0]
			xfit = linspace(times[x_max2],times[x_end],100)
			yfit = A1*exp(-b1*(xfit-x01))+A2*exp(-b2*(xfit-x02))
			print pls[0]
			print ""

			x_half = SAll[:,0]

			half1 = (y_max)/2.+0*x_half
			half2 = (y_max2)/2.+0*x_half

			figure(5);clf() 
			plot(times,SAll[:,1],'r.-') 
			plot(times,half2,'-k')
			plot([times[x_25],times[x_75]],[y_25,y_75],'k-o')
			plot([times[x_max],times[x_max2]],[y_max,y_max2],'ko')
			plot(xfit,yfit,'b-')
			legend(['Cer','half of max','slope','maxima','exp fit'],loc='right', bbox_to_anchor=(1.1,1),prop={'size':10}) 
			xlabel("times in s") 
			ylabel("gray value") 
			title("Analyse in this figure 5! (x=%d  y=%d)"%(x,y)) 
			draw();show(block=False) 
			print("Please check whether both maxima, as well as the slope are marked correctly!")
			a = raw_input("\'Enter\' to continue! ") 
			x_half2 = FR*askVal("x-value for the 2nd intersect (black): ")

			print ""
			# TODO Check whether the frame rate is correct everywhere

			d_inc = (x_max-x1)/FR
			d_dec = (x_half2-x_max2)/FR
			d_half = (x_half2-(x1+x_max))/(FR*2.)

			print "d_inc: %.1f\t d_dec: %.1f\t d_half: %.1f slope: %.1f"%(d_inc,d_dec,d_half,slope)
			
			cut = max(x1-8,0)
			Hec = (SAll[cut:,1]-SAll[x1,1])/(SAll[x_max,1]-SAll[x1,1])

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
			
			SAll[:,2] = (SAll[:,2]-SAll[x4,2])/SAll[x4,2]
			SAll[:,3] = (SAll[:,3]-SAll[x4,3])/SAll[x4,3]
			
			P2S = (SAll[cut:,2]-SAll[x4,2])/SAll[x4,2]
			PS = (SAll[cut:,3]-SAll[x4,3])/SAll[x4,3]
			

			# integration
			integr1 = 1.0*sum(SAll[x1:x_end,1])/FR
			# Ask Annita for the other integration
			integr2 = 1.0*sum(SAll[x1:x_end,2])/FR
			integr3 = 1.0*sum(SAll[x1:x_end,3])/FR

			figure(5);clf()
			plot(Hec,'r.-')
			plot(P2S,'b.-')
			plot(PS,'g.-')
			legend(['Hec','P2S','PS'],loc='right', bbox_to_anchor=(1.1, 1),prop={'size':10}) 
			draw();
			
			# Save the dx:
			# Ask Annita what to save
			HalfOut = open(HalfTimeFileName,"a")
			HalfOut.write("%d\t| %d\t| %.2f\t| %.2f\t| %.2f\n"%(x,y,d_inc,d_dec,d_half))
			HalfOut.close()

			SlopOut = open(SlopeFileName,"a")
			SlopOut.write("%.2f\t| %.2f\t| %.2f\t| %.2f\t| %.2f \n"%(x_25,y_25,x_75,y_75,slope))
			SlopOut.close()
			
			IntOut = open(IntFileName,"a")
			IntOut.write("%.2f\t| %.2f\t| %.2f\t| %.2f\t| %.2f \n"%(integr1,integr2,integr3))
			IntOut.close()
			
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
			a = raw_input("\'Enter\' to continue, \'q\' to quit! ") 
			 
	elif (a=='q') or (a=='x') or (a=='Q') or (a=='X'): 
		print "Stop!" 
		break; 
	elif a=='': 
		print "Continue! \n\n\n\n" 
 
 
	disconnect(cid1); 
	disconnect(cid2) 
