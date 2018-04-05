#!/u/home/stwe/anaconda/bin/python
import glob
import os.path
import Image
from pylab import *

close("all")
curr_pos= 0

# Name of the DataSet
#DataSet = raw_input("DataSet: ")
DataSet='BDNF_Phl1024-3'

# The Path of the DataSet - should be changed if you use your external hard drive
PATH = '/media/My Passport/HOLZLABDATA/ANALYZEDPTIRFData/BDNF_Phl10242014/'+DataSet+"/"
DataDir = PATH+'GreyData/'

files = glob.glob(DataDir+'/*.txt')
files.sort()

if size(files)==0:
  print "\n\t *** Error *** !"
  print '\n \tNo .txt files found under: \n '+DataDir+' \n\n '
  quit()

N=size(files)
count=0
for file in files:
	count+=1
	print ("%d/%d Load: %s")%(count,N,file[-20:])
	T = genfromtxt(file)

	figure(1);clf();
	plot(T[:,0],T[:,1],'r-')
	plot(T[:,0],T[:,1],'r.')
	title("HeC - grey value")
	xlabel("frame no")
	ylabel("grey value")
	draw();show(block=False);
	print "\n ****************"
	x=input("Choose the starting frame (x-value): ")	
	print "x: %d"%x
	print "Good choice!!"
	
	TN = T[x-8:,:]
	TN[:,1]=TN[:,1]/T[x,1]-1
	TN[:,2]=TN[:,2]/T[x,2]-1
	TN[:,3]=TN[:,3]/T[x,3]-1

	figure(2);clf();
	plot(TN[:,0],TN[:,2],'b-',label='P2S')
	plot(TN[:,0],TN[:,3],'g-',label='PS')
	plot(TN[:,0],TN[:,2],'b.')
	plot(TN[:,0],TN[:,3],'g.')
	legend(['P2S','PS'])
	xlabel("original frame no")
	ylabel("grey value")
	draw();show(block=False)

	figure(1);clf();
	plot(TN[:,1],'r-')
	plot(TN[:,1],'r.')
	xlabel("new frame no")
	ylabel("grey value")
	title("HeC")
	draw()
	
	print "\n*************"
	xmax = input("Where is the maximum (x-value): ")
	ymax = TN[xmax,1]
	yHalf = 0.5*ymax+0*TN[:,1]
	plot(yHalf,'-k')
	draw()
	print "\n****************"
	char = raw_input("Read the x value and do whatever you want!\nEnter to continue! \'q\' to quit. ")
	if char == 'q':
		break
