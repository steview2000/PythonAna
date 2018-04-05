pro Ivstime3

;   6-24-08:  This program:
;       inputs a stack of time-sequence images of cells from the Andor in tif form,
;       allows the user to define a set of points (maximum # nroi) in any image in the stack (which apply to the whole stack),
;       calculates the average intensity in a circle around those points,
;       and plots average intensity vs. frame number.

;       The outputs are text files, screen graphs, and one image showing the locations of all the ROI's.  The screen graphs and image
;       can be saved by the keyboard "Print Screen"

;       This program is similar to parts of avgintcircle.pro

;   6-26-08  Add possibility of multiple composite movie read-in and display, up to three side-by-side, (e.g., Cer, P/S, and P2S),
;       allow ROI selection on the first (left-most) one, and calculate avg ROI intensity on all three stacks.
;   6-30-08  Some bugs removed.  File output now starts with frame # at which roi was selected, indicating (say) exocytosis frame.
;   9-20-08  Cursor boxes shown in middle and right panels during selection search.
;   10-13-08  Imageblock now defined as an unsigned long64 integer array to avoid P+2S from overflowing
;   10-13-08  The side-by-side display in animateselectmulti is now a composite of byte-scaled versions of each panel.
;                This avoids one panel saturating (usually the P+2S) while the others are too dark.
;   10-14-08  Negative textfile output at very bright P+2S fixed with ulonarr type for current1,current2,current3 arrays.
;   10-16-8   MPEG Movie option added.
;   10-21-08  Rows flipped (IDL command REVERSE) in mpeg file to agree with original tif files.
;   10-22-08  Add option (set by preavg=1) to calculate <P>/<S> instead of <P/S> where <> is avergae in each roi
;   10-25-08  add option for locking in contrast adjustments for mpg movie.
;   10-25-08  Fixed slight mistake in lockin procedure:  In the animatemultiselect subroutine,
;              lockin=1 locks in" the contrast adjustments to apply to the movie.  Set to 0 for no lock-in.
;   10-25-13  Option for IDL 8 (actually added at an earlier date) and option for multimage tiff files added
;   1-31-15   Add option to read in NON-thresholded P/S.  Change to always read in four stacks: ceru, P2Ssum,P/S w/thresh and P/S NONthresh
;				Major changes in all variable names and ordering.
;   2-16-15   Fixed timing issues with playback.  Recommend using Loupe while selecting ROI's.  Allow change of order (priority variable)
;               so ROI's can be selected on the P/S panel.  Allow reversal into a "negative" for easier viewing.  Make ROI choices erasable.
;   3-2-15    Error in read-in of PSratioNON now fixed.
;   3-10-15	  Removed option for IDL8
;				Image panel use for ROI selection expanded 2x
;				Only P/S and cer displayed during ROI selection but at double size
;               Jpg output of graphs, superimpose Cer, P/S, P+2S. Labeleed as (parent)-ALL-ROI (#).jpg
;               Input PSratNON as obtained from tif restored to actual values by dividing by 10,000 to counter multiplier in PSratiodual
;				Multicol output, one col for each roi, added.  Single col roi files still preserved.
;   3-11-15   Adjusted txt jpg and screen output formats
;   3-11-15   Additional format fix
;   3-13-15.  Calculate average over each ROI, not total as before.  Change to new method of isolating ROI's with masks, truncation,
;				and array (rather than pixel) math, increases calc speed by 25x !

path='C:\users\dan\My Documents\pstest\'
multimage=1  ;  =1 for input multi-image tif files;  =0 for stack of separate tif files for each image

aspect=1  ;  (=0 for screen and jpg graphical output square 1000x1000;  =1 for 1800x1000
priority=1 ;(=0 for roi selection on cerulean; =1 for roi selection on P/S)
BWreverse=0 ; (=0 for images to have black where countsare low; =1 for white where counts are low; i.e., a 'negative')

radius=4  ;  the radius of the area to be used for intensity calc around each roi, and 1/2 side length of the selected box display
nstacks=2  ;  # of side-by-side images to be shown in selecting ROI's
label1='Cerul'
label2='P2Ssum'
label3='PSratioTHR'
label4='PSratioNON'

device,decomposed=0
loadct,0

tryagain:

;  The stacks must be the same size in all three dimensions:
;  'file_' is the full file name including path and extension (.tif); parentname is the same but w/o the extension
print,'Select CERUL tif stack.
if multimage eq 0 then opentfiles,parentname1,file1,ceru,nframes1,xsize1,ysize1,1,path
if multimage eq 1 then opentmulti,parentname1,file1,ceru,nframes1,xsize1,ysize1,1,path

 ;  Cleave off the 5-character CERUL part of the file name for later use in identifying the set of files
length=STRLEN(parentname1)
parent=STRMID(parentname1,0,length-5)
print,'parent=',parent

print,'Select P2Ssum tif stack.'
if multimage eq 0 then opentfiles,parentname2,file2,P2Ssum,nframes2,xsize2,ysize2,2,path
if multimage eq 1 then opentmulti,parentname2,file2,P2Ssum,nframes2,xsize2,ysize2,2,path
if xsize2 ne xsize1 or ysize2 ne ysize1 or nframes2 ne nframes1 then begin
	print,'Size and depth of stacks are not equal'
    goto,tryagain
endif

print,'Select PSratioNON (no thresh) tif stack.'
if multimage eq 0 then opentfiles,parentname4,file4,PSratNON,nframes4,xsize4,ysize4,4,path
if multimage eq 1 then opentmulti,parentname4,file4,PSratNON,nframes4,xsize4,ysize4,4,path
if xsize4 ne xsize1 or ysize4 ne ysize1 or nframes4 ne nframes1 then begin
  	print,'Size and depth of stacks are not equal'
   	goto,tryagain
endif

;  Note that the input PSratNON.tif input file has probably had its values pre-multiplied by 10,000 before
;    recording by PSratiodual2.pro.  So restore the correct values here
PSratNON=PSratNON/10000.

print,'Select PSratioTHR (with/thresh) tif stack.'
if multimage eq 0 then opentfiles,parentname3,file3,PSratTHR,nframes3,xsize3,ysize3,3,path
if multimage eq 1 then opentmulti,parentname3,file3,PSratTHR,nframes3,xsize3,ysize3,3,path
if xsize3 ne xsize1 or ysize3 ne ysize1 or nframes3 ne nframes1 then begin
  	print,'Size and depth of stacks are not equal'
   	goto,tryagain
endif

;  Temp test of selection accuracy.  Make images with just two bright pixels: at (100,100) and (111,111) (one even, one odd)
;  ***************************

;xsize1=327
;ysize1=328
;xsize2=327
;ysize2=328
;xsize3=327
;ysize3=328
;xsize4=327
;ysize4=328
;nframes1=100
;nframes2=100
;nframes3=100
;nframes4=100

;ceru=UINTARR(xsize1,ysize1,nframes1)
;P2Ssum=UINTARR(xsize2,ysize2,nframes2)
;PSratTHR=UINTARR(xsize3,ysize3,nframes3)
;PSratNON=UINTARR(xsize4,ysize4,nframes4)

;ceru(100,100,*)=255
;P2Ssum(100,100,*)=255
;PSratTHR(100,100,*)=255
;PSratNON(100,100,*)=255

;ceru(111,111,*)=255
;P2Ssum(111,111,*)=255
;PSratTHR(111,111,*)=255
;PSratNON(111,111,*)=255

; ****************************

cerufrm=fltarr(xsize1,ysize1)
P2Ssumfrm=fltarr(xsize1,ysize1)
PSratNONfrm=fltarr(xsize1,ysize1)
Sfrm=fltarr(xsize1,ysize1)
Pfrm=fltarr(xsize1,ysize1)

proceed:

nframes=nframes1
print,'nframes=',nframes
ysize=ysize1
xsize2=xsize1*2  ;  width of array with side-by-side P/S and cerul

;  Make left-to-right composite viewblock for display during ROI selection;
;    viewblock will show the individually byte-scaled data in each of the stacks for better side-by-side viewing

viewblock=bytarr(xsize2,ysize,nframes)

if priority eq 0 then begin
	viewblock[0:xsize1-1,*,*]=bytscl(ceru[0:xsize1-1,*,*])
	viewblock[xsize1:2*xsize1-1,*,*]=bytscl(PSratTHR[0:xsize1-1,*,*])
endif
if priority eq 1 then begin
	viewblock[0:xsize1-1,*,*]=bytscl(PSratTHR[0:xsize1-1,*,*])
	viewblock[xsize1:2*xsize1-1,*,*]=bytscl(ceru[0:xsize1-1,*,*])
endif
if BWreverse eq 1 then viewblock=255-viewblock

;  SELECT ROI positions

;  Subroutine animateselectmulti is similar to movies, except it is a succession of IDL WINDOW 0 frames.
;  It is exited by hitting any keyboard key.

animateselectmulti,viewblock,nframes,xsize2,ysize,frame,xpos,ypos,radius,nroi,nstacks,modimage,priority,BWreverse

cer=fltarr(nroi,nframes)
sum=fltarr(nroi,nframes)  ; This is the final array for P+2S
rat=fltarr(nroi,nframes)  ; This is the final array for P/S
cerplot=fltarr(nframes)   ; for plotting
sumplot=fltarr(nframes)
ratplot=fltarr(nframes)

print,'# or ROIs selected = ',nroi
bigw=nroi*16
widthout=80
if bigw gt 80 then widthout=bigw
print,'widthout=',widthout

xpos=xpos/2
ypos=ypos/2

loadct,2  ;  enable use of colors on plots
tvlct,red,green,blue,/get
red[216]=0
green[216]=0
blue[216]=255
modifyct,2,'blue213',red,green,blue  ;  replace color 255 in rainbow palette with white
loadct,2

;Set up circular mask in 3D (i.e., a cylinder) for isolating an ROI

circle3=intarr(xsize1,ysize,nframes)
ntr=2*radius+1  ; number of pixels in the side of a square that trunactes around the circle
cershifttr=intarr(ntr,ntr,nframes)
circle3tr=intarr(ntr,ntr,nframes)
Pshifttr=intarr(ntr,ntr,nframes)
Sshifttr=intarr(ntr,ntr,nframes)

xmid=xsize1/2
ymid=ysize/2
for ix=xmid-radius,xmid+radius do for iy=ymid-radius,ymid+radius do if (ix-xmid)^2+(iy-ymid)^2 lt radius^2 then circle3[ix,iy,*]=1.
circle3tr[0:ntr-1,0:ntr-1,*]=circle3[xmid-radius:xmid+radius,ymid-radius:ymid+radius,*]
inmask=where(circle3tr eq 1.,npix)  ;  Figure out how many pixels npix in the circle. Inmask is a dummy, not used again

;  Treat the whole stack to recalculating P and S, including all frames
pos=where(PSratNON gt 0.)
S=P2Ssum    ;  This sets up the P and S array dimensions just like their sources
S[*,*,*]=0.
P=PSratNON
P[*,*,*]=0.

S[pos]=P2Ssum[pos]/(PSratNON[pos]+2.);S and P are pixel-by-pixel arrays over all frames in the stack.  Pixels outside of pos remain at 0
P[pos]=PSratNON[pos]*S[pos]

;  Shift the image files to coincide with the mask.  This must be done separately for each roi

for iroi=0,nroi-1 do begin
	x=xpos[iroi]
    y=ypos[iroi]

	cershift=SHIFT(ceru,xmid-x,ymid-y,0)  ; cerroi is the 3D masked cer fluor in each and all frames at a particular roi
	Sshift=SHIFT(S,xmid-x,ymid-y,0)		; '0' in last argument means no shift in the frame# dimension
	Pshift=SHIFT(P,xmid-x,ymid-y,0)

;		Truncate these arrays to the region around xmid,ymid only;
;
	cershifttr[0:ntr-1,0:ntr-1,*]=cershift[xmid-radius:xmid+radius,ymid-radius:ymid+radius,*]
	Sshifttr[0:ntr-1,0:ntr-1,*]=Sshift[xmid-radius:xmid+radius,ymid-radius:ymid+radius,*]
	Pshifttr[0:ntr-1,0:ntr-1,*]=Pshift[xmid-radius:xmid+radius,ymid-radius:ymid+radius,*]

;  Caculate avgs in product of mask and these truncated images at this roi, through whole stack

	cerroi=circle3tr*cershifttr ;  because circle3tr was mad floating point, this calc should be FP
	Sroi=circle3tr*Sshifttr   ;  These are still 3D arrays
	Proi=circle3tr*Pshifttr
	for iframe=0,nframes-1 do begin
		cerfrm=TOTAL(cerroi[*,*,iframe])/npix  ;  These cerfrm,Sfrm,Pfrm are scalars
		Sfrm=TOTAL(Sroi[*,*,iframe])/npix
		Pfrm=TOTAL(Proi[*,*,iframe])/npix
    	cer[iroi,iframe]=cerfrm
		sum[iroi,iframe]=Pfrm+2*Sfrm  ;  nthese array are (nroi,nframes)
		rat[iroi,iframe]=Pfrm/Sfrm
	endfor
endfor

;  Plot graphs in windows, three plots (cer, P/S, P+2S in one graph, new graph for each ROI.  The scale for P/S is absolute.
;     The other two: cer and P+2S are normalized to their max.

 for iroi=0,nroi-1 do begin
    plotname=label1+' norm (RED)   '+label4+' (BLACK)   '+label2+' norm (BLUE)'
    graphfilename=parent+'-ALL- ROI '+STRTRIM(iroi,2)+'.jpg'
    if aspect eq 0 then xaspect=1000 else xaspect=1800
    window,title=plotname,xsize=xaspect,ysize=1000,xpos=0,ypos=20*iroi,/free
    cerplot[*]=cer[iroi,*]
	ratplot[*]=rat[iroi,*]
	sumplot[*]=sum[iroi,*]
	cerplotnorm=cerplot/max(cerplot)
	sumplotnorm=sumplot/max(sumplot)
	;ratplotnorm=ratplot/max(ratplot)
	;  Plot based on the max values.  The norm cer and P+2S will have a max of 1, so use ratplot to set the max
	if max(ratplot) gt 1. then mxplot=max(ratplot) else mxplot=1.0
	plot,ratplot,background=255,color=0,thick=2,yrange=[0,mxplot],font=-1  ;  black
	oplot,cerplotnorm,color=72,thick=2				 ;	red
	oplot,sumplotnorm,color=216,thick=2				 ;  blue
    plots,[frame[iroi],frame[iroi]],[0,mxplot],color=15  ;  dark green
    XYOUTS,160,60,'ROI# '+STRTRIM(iroi,2)+'    Frame '+STRTRIM(frame[iroi],2)+'    ImageJ position: ('+STRTRIM(FIX(xpos[iroi]),2)+$
    	','+STRTRIM(ysize1-FIX(ypos[iroi]),2)+')',color=0,/device,charsize=2
	XYOUTS,160,90,plotname,color=0,/device,charsize=2
	XYOUTS,160,120,parent,color=0,/device,charsize=2
    allgraph=tvrd(true=1)
    write_jpeg,graphfilename,allgraph,true=1
endfor

;  Create a file of avgintroi for later replotting for each stack.
;  Two types:   single ROI files, and multi-ROI multicols.  Make format compatible with secrate3.pro multicol input

;  Single ROI: one file for each ROI

for iroi=0,nroi-1 do begin
    cerplot[*]=cer[iroi,*]
    textout1=parentname1+' vs t roi'+STRTRIM(iroi,2)+'.txt'
    openw,lun,textout1, /get_lun
    printf,lun,frame[iroi],xpos[iroi],ysize-ypos[iroi]
    for iframe=0,nframes-1 do printf,lun,iframe,cerplot[iframe]
    close,lun
    free_lun,lun
endfor

for iroi=0,nroi-1 do begin
    sumplot[*]=sum[iroi,*]
    textout2=parentname2+' vs t roi'+STRTRIM(iroi,2)+'.txt'
    openw,lun,textout2, /get_lun
    printf,lun,frame[iroi],xpos[iroi],ysize-ypos[iroi]
    for iframe=0,nframes-1 do printf,lun,iframe,sumplot[iframe]
    close,lun
    free_lun,lun
endfor

for iroi=0,nroi-1 do begin
    ratplot[*]=rat[iroi,*]
    textout4=parentname4+' vs t roi'+STRTRIM(iroi,2)+'.txt'
    openw,lun,textout4, /get_lun
    printf,lun,frame[iroi],xpos[iroi],ysize-ypos[iroi]
    for iframe=0,nframes-1 do printf,lun,iframe,ratplot[iframe]
    close,lun
    free_lun,lun
endfor

;  Multicol:  a single file containing all the ROI's in successive columns for each type: cerul, P/S, P+2S

textcerul=parent+'-ALL-'+label1+'.txt'
textPSrat=parent+'-ALL-'+label4+'.txt'
textP2Ssum=parent+'-ALL-'+label2+'.txt'

 ;  The arrays for the data at each time point is a vector of results at each roi

cerroi=fltarr(nroi)
PSratroi=fltarr(nroi)
P2Ssumroi=fltarr(nroi)
ROIroi=indgen(nroi)
frameroi=intarr(nroi)
frameroi[0:nroi-1]=frame[0:nroi-1]
xposroi=intarr(nroi)
xposroi[0:nroi-1]=xpos[0:nroi-1]
yposroi=intarr(nroi)
yposroi[0:nroi-1]=ypos[0:nroi-1]

;  Cerulean multicol file:

openw,lun,textcerul, /get_lun,width=widthout
printf,lun,'ROI              ',STRTRIM(ROIroi,2)+  '           '
printf,lun,'Frame            ',STRTRIM(frameroi,2)+'          '
printf,lun,'xpos           ',STRTRIM(xposroi,2)+   '         '
printf,lun,'ypos           ',STRTRIM(yposroi,2)+   '         '
for iframe=0,nframes-1 do begin
	cerroi[*]=cer[*,iframe]  ; where the * means all of the roi position
	printf,lun,iframe,cerroi
endfor
close,lun
free_lun,lun

;  Cerulean multicol file:

openw,lun,textPSrat, /get_lun,width=widthout
printf,lun,'ROI              ',STRTRIM(ROIroi,2)+  '           '
printf,lun,'Frame            ',STRTRIM(frameroi,2)+'          '
printf,lun,'xpos           ',STRTRIM(xposroi,2)+   '         '
printf,lun,'ypos           ',STRTRIM(yposroi,2)+   '         '
for iframe=0,nframes-1 do begin
	PSratroi[*]=rat[*,iframe]  ; where the * means all of the roi position
	printf,lun,iframe,PSratroi
endfor
close,lun
free_lun,lun

;  Cerulean multicol file:

openw,lun,textP2Ssum, /get_lun,width=widthout
printf,lun,'ROI              ',STRTRIM(ROIroi,2)+  '           '
printf,lun,'Frame            ',STRTRIM(frameroi,2)+'          '
printf,lun,'xpos           ',STRTRIM(xposroi,2)+   '         '
printf,lun,'ypos           ',STRTRIM(yposroi,2)+   '         '
for iframe=0,nframes-1 do begin
	P2Ssumroi[*]=sum[*,iframe]  ; where the * means all of the roi position
	printf,lun,iframe,P2Ssumroi
endfor
close,lun
free_lun,lun

print,'Do you want to make an mpeg movie ?
print,  '0=NO'
print,  '1=YES: WITH selection boxes'
read,   '2=YES: WITHOUT selection boxes',movie
if movie eq 0 then goto,finish
read, 'FIRST frame # in movie = ',istart
read, 'LAST  frame # in movie = ',iend
if movie eq 1 then mpegdirect,modimage,parentname1,nframes,xsize,ysize,istart,iend
if movie eq 2 then mpegdirect,viewblock,parentname1,nframes,xsize,ysize,istart,iend

finish:print,'Program ALL DONE !'
for i=0,31 do wdelete,i

end

;*********************************************************************************************
pro opentfiles,parentname,file,imageblock,nframes,xsize,ysize,istack,path

;  This program is adapted from openbfiles.pro to read Andor tif files, which have the same numbering scheme as Sensicam b16 files

;  modified 8-29-06 to obtain parentname by truncating off last 8 characters in filename for stacks, and correctly handle single
;     files without truncation.
;  12-22-06 Fixed reporting of nframes
;  1-25-08  Modified to read in stacks that might not start at filename # 0000
;  1-25-08  Fixed bug: ifirst # must be explicitly converted to integer
;  6-26-08  Modified to change title in fileselection window according to istack, to keep track of multiple stack read-ins

if istack eq 1 then file = DIALOG_PICKFILE(TITLE='Select CERU stack',PATH=path,FILTER = '*.tif')
if istack eq 2 then file = DIALOG_PICKFILE(TITLE='Select P2Ssum stack',PATH=path,FILTER = '*.tif')
if istack eq 3 then file = DIALOG_PICKFILE(TITLE='Select PSrat (w/thresh)',PATH=path,FILTER = '*.tif')
if istack eq 4 then file = DIALOG_PICKFILE(TITLE='Select PSrat NONthresh',PATH=path,FILTER = '*.tif')

;  Parameter single tracks whether we have a single (=1) or multiple set (=0) of files.
;  Single files do not have a four-digit xxxx frame number. If the total length of a filename is le 8 places (...xxxx.tif)
;  then it must be a single with a filename structure header.tif.  Single files also exist with longer than 8 places.
;  These will be discovered by adding 0000 to the header and finding that such a name doe not exist.

single=0
length=STRLEN(file)
help,file
help,length
if length le 8 then begin
    single=1
    header=STRMID(file,0,length-4)
endif else header=STRMID(file,0,length-8)

blah=query_tiff(file,info)
;info is a structure that contains all sorts of info about the image
;you can get this info when you use the query function
xsize=info.dimensions[0]
ysize=info.dimensions[1]
print,xsize,ysize

;  Determine # of images

;read,'Frame number of FIRST frame = ', ifirst
ifirst=0
ifirst=FIX(ifirst)
for iframe=ifirst,9999 do begin
    if iframe lt 10 then insert='000'
    if iframe ge 10 and iframe lt 100 then insert='00'
    if iframe ge 100 and iframe lt 1000 then insert='0'
    if iframe ge 1000 and iframe lt 10000 then insert=''
    filenew = header + insert + strtrim(iframe,2) + '.tif'
    print,filenew
    exist=query_tiff(filenew)
    if exist eq 0 then goto,getout
endfor

getout:if iframe eq 0 then begin
    single=1
    nframes=1
    header=STRMID(file,0,length-4)
    print,'single file header = ',header
endif else nframes=iframe-ifirst
print,'nframes = ', nframes

;Set up stack of images
imageblock=UINDGEN(xsize,ysize,nframes)
singleimage=UINDGEN(xsize,ysize)

;Place images in the stack
for iframe=ifirst,nframes-1+ifirst do begin
    if iframe lt 10 then insert='000'
    if iframe ge 10 and iframe lt 100 then insert='00'
    if iframe ge 100 and iframe lt 1000 then insert='0'
    if iframe ge 1000 and iframe lt 10000 then insert=''
    filenew = header + insert + strtrim(iframe, 2) + '.tif'
    exist=query_tiff(filenew)
    if exist eq 1 then begin
        singleimage=read_tiff(filenew)
        imageBlock[*,*,iframe-ifirst]=REVERSE(singleimage,2) ; an inversion of the y-axis to agree with Andor and TIF output
    endif else begin
        nframes=iframe+1-ifirst
        goto,moveon
    endelse
endfor
moveon:help,nframes
if single eq 1 then begin
    imageblock[*,*,0]=read_tiff(header+'.tif')
    imageblock=REVERSE(imageblock,2)
endif

parentname=header

end

; ***************************************************************************************************
pro opentmulti,parentname,file,imageblock,nframes,xsize,ysize,istack,path

;  This program is adapted from opentfiles.pro to read Andor multi-image tif files

if istack eq 1 then file = DIALOG_PICKFILE(TITLE='Select CERU stack',PATH=path,FILTER = '*.tif')
if istack eq 2 then file = DIALOG_PICKFILE(TITLE='Select P2Ssum stack',PATH=path,FILTER = '*.tif')
if istack eq 3 then file = DIALOG_PICKFILE(TITLE='Select PSrat (w/thresh)',PATH=path,FILTER = '*.tif')
if istack eq 4 then file = DIALOG_PICKFILE(TITLE='Select PSrat NONthresh',PATH=path,FILTER = '*.tif')


blah=query_tiff(file,info)
;info is a structure that contains all sorts of info about the image
;you can get this info when you use the query function
xsize=info.dimensions[0]
ysize=info.dimensions[1]
print,xsize,ysize

nframes=info.NUM_IMAGES
length=STRLEN(file)
parentname=STRMID(file,0,length-4)
print,'TIF filename  = ',parentname
print,'nframes = ',nframes
imageblock=UINDGEN(xsize,ysize,nframes)

for iframe=0,nframes-1 do begin
	singleimage=read_tiff(file,image_index=iframe)
	imageBlock[*,*,iframe]=REVERSE(singleimage,2)
endfor
help,imageblock

end


; ************************************************************************************

pro animateselectmulti,imageblock,nframes,xsize,ysize,frame,xpos,ypos,box,jroi,nstacks,modimage,priority,BWreverse

;
;   This makes a movie, like movies.pro, but instead of using xinteranimate, it sequentially
;   displays images with tvscl. It stops by hitting anything on the keyboard (e.g., the spacebar)
;   This also allows use of that image for identifying the center of an roi
;
;   10-19-07  This is modified from animate.pro (as found in avgintroi.pro).  It does the following:
;     Outputs the first frame # ff to be used in tracking a particular granule.
;     Remove parentname argument, which was not used anyway.
;     Add ability to manually step forward or back
;     Adds ability to adjust sppeed of animation
;     Add ability to select a granule with the cursor
;     Adds ability to select a point near a granule during an adjacent bright "cloud" burst
;       for later background subtraction.
;   6-24-08  Modified from animateselect to allow multiple ROI selection within the subroutine.
;       'Cloud' identification deleted.  Expansion options deleted.  Selected box shown within this subroutine.
;   10-14-08.  Note that "imageblock" in this subroutine is "viewblock", the byte-scaled version of imageblock in the main program.
;   10-15-08.  The modified stack (with boxes) is now outputted to the main program, for possible mpeg making.
;   2-16-15  Various changes described in comments for Ivstime2 main program
;   3-9-15	 Various changes in wset, wait, and statement order to make smoother operation


;  Possible future redefinition of cursor symbol
;strArray = [ $
;       'oooooooooooooooo', $
;       'o              o', $
;		'o              o', $
;		'o              o', $
;		'o              o', $
;		'o              o', $
;		'o              o', $
;		'o       $      o', $
;		'o              o', $
;		'o              o', $
;		'o              o', $
;		'o              o', $
;		'o              o', $
;		'o              o', $
; 		'o              o', $
; 		'oooooooooooooooo']
;
;cursor_image = CREATE_CURSOR(strArray, HOTSPOT=hotspot);, MASK=mask)
;window,31
;tvscl,cursor_image
;REGISTER_CURSOR, 'translate', cursor_image, HOTSPOT=hotspot,/overwrite; MASK=mask,
;self->RegisterCursor, strArray, 'LUT', /DEFAULT
;obj->IDLgrWindow::SetCurrentCursor,translate
;DEVICE,CURSOR_IMAGE=cursor_image; CURSOR_MASK=value{WIN, X}]

nroimax=100 ;  The max # of roi positions that can be chosen
maximorig=MAX(imageblock)
minimorig=MIN(imageblock)
lockin=1 ;  lockin=1 locks in" the contrast adjustments to apply to the movie.  Set to 0 for no lock-in.

print,''
print,'f = FORWARD direction'
print,'r = REVERSE direction'
print,''
print,'+ = AUTO play, FASTER'
print,'- = AUTO play, SLOWER'
print,''
print,'m = MANUAL play'
print,'SPACEBAR = manual SINGLE STEP'
print,''
print,'b = BRIGHTER:      decrease upper limit '
print,'d = DARKER:        increase upper limit  '
print,'c = MORE CONTRAST: increase lower limit  '
print,'z = LESS CONTRAST: decrease lower limit  '
print,''
print,'GO to a frame of interest(i.e., the frame just before exocytosis),'
print,'   then:'
print,'LEFT  click to CREATE a ROI'
print,'x = DELETE last ROI choice'
print,''
print,'RIGHT click when DONE making ROI selections'

waitdelay=0.; sets the default speed (0 is maximum)
direct=1 ;    sets the default direction (direct =+1 for forward, -1 for reverse)
auto=0;       sets the default mode to auto
done=0
i=0  ;  this keeps track of the displayed frame #
jroi=0  ;  this is the number assigned to the roi
frame=intarr(nroimax)  ;  the frame # at which an roi is selected.
xpos=fltarr(nroimax)
ypos=fltarr(nroimax)

xsize2x=2*xsize  ;  These steps prepare to enlarge the image by a factor of 2
ysize2x=2*ysize
box2x=2*box

modimage=imageblock  ;  modimage will have cumulative white boxes drawn upon it.  imageblock is input set of 3 side-by-side images
modimage=REBIN(modimage,xsize2x,ysize2x,nframes);  This expands the images by 2x because of the redefinition of xsize and ysize

if priority eq 0 then title='Hit SPACEBAR, then LEFT click on left to create ROIs   CERULEAN                                               P/S w thresh'
if priority eq 1 then title='Hit SPACEBAR, then LEFT click on left to create ROIs   P/S w thresh                                             CERULEAN'

window,0,xsize=xsize2x,ysize=ysize2x,title=title,xpos=0,ypos=50

minim=minimorig
maxim=maximorig
xuse=0  ;  set up this variable to count # of times that x option is pressed repeatedly
j=0
repeat B=GET_KBRD(0) until B eq '' ;  this clears the keyboard
repeat begin

    !mouse.button=0

    wset,0
    j=j+1
    running='frame '+STRTRIM(FIX(i),2)
    scaledimage=BYTSCL(modimage[*,*,i],max=maxim,min=minim)  ;  this resets scaled image to show box clicked on below
    tv,scaledimage
    XYOUTS,10,10,running,color=255,/device,charsize=2

    if j eq 1 then A=GET_KBRD() else A=GET_KBRD(0)

    if A eq 'b' then begin
        range=maxim-minim
        maxim=maxim-0.1*range
        if maxim le minim then maxim=maxim+0.1*range
    endif
    if A eq 'd' then begin
        range=maxim-minim
        maxim=maxim+0.1*range
        if maxim gt 32768 then maxim=maxim-0.1*range
    endif
    if A eq 'c' then begin
        range=maxim-minim
        minim=minim+0.05*range
        if minim ge maxim then minim=minim-0.05*range
    endif
    if A eq 'z' then begin
        range=maxim-minim
        minim=minim-0.05*range
        if minim lt 0 then minim=minim+0.05*range
    endif

    if A eq 'f' then direct=1
	if A eq 'r' then direct=-1
	if A eq 'm' then auto=0

    i=i+auto*direct

    if A eq ' ' and j gt 1 then i=i+direct
	if i eq nframes then i=0
    if i eq -1 then i=nframes-1

    ;auto speed control, also returns mode to auto if it was not already there

    if A eq '+' then begin
        if waitdelay ge 0.1 then waitdelay=waitdelay-0.1
        auto=1
    endif
    if A eq '-' then begin
         if waitdelay le 1.0 then waitdelay=waitdelay+0.1
         auto=1
    endif
    if auto eq 1 then wait,waitdelay

    cursor,x,y,0,/device ;  this activates the mouse buttons and reports the current x,y position
    ;
    for istack=2,nstacks do begin     ;  This sets up and shows the 3-wide image with moveable box that does not stick to a position
        xs=(istack-1)*xsize2x/nstacks+x
        if xs-box2x ge 0 and xs+box2x lt xsize2x and y-box2x ge 0 and y+box2x lt ysize2x then begin
            scaledimage[xs-box2x:xs+box2x,y-box2x]=255*(1-BWreverse) ; ensures that box is white for positive and black for a negative image
            scaledimage[xs-box2x:xs+box2x,y+box2x]=255*(1-BWreverse)
            scaledimage[xs+box2x,y-box2x:y+box2x]=255*(1-BWreverse)
            scaledimage[xs-box2x,y-box2x:y+box2x]=255*(1-BWreverse)
        endif
    endfor

    tv,scaledimage

    if !mouse.button eq 1 then begin  ; left click
        newselect:cursor,x,y,4,/device  ;  records x,y on up-click
        if x gt xsize2x/nstacks then begin
            print,'x coordinate must be within the left-most image'
            goto,newselect
        endif
        frame[jroi]=i
        xpos[jroi]=x
        ypos[jroi]=y
        loc='ROI #'+STRTRIM(jroi,2)+' selected at frame #'+STRTRIM(i,2)+'  at ImageJ location ('+STRTRIM(x/2,2)+','+STRTRIM((ysize2x-y)/2,2)+')'
        print,loc
        xuse=0
        wset,0

; Create box (LEFT click) in all side-by-side image stacks
		oldmodimage=modimage  ; the image before recording the box, in case we want to cancel
        for istack=1,nstacks do begin
            xs=(istack-1)*xsize2x/nstacks+x
            if BWreverse eq 0 then boxcol=maxim else boxcol=minim
            modimage[xs-box2x:xs+box2x,y-box2x,*]=boxcol
            modimage[xs-box2x:xs+box2x,y+box2x,*]=boxcol
            modimage[xs+box2x,y-box2x:y+box2x,*]=boxcol
            modimage[xs-box2x,y-box2x:y+box2x,*]=boxcol
        endfor

        jroi=jroi+1
    endif

;  EXIT selection procedure

    if !mouse.button eq 4 then done=1  ; right click
    if A eq 'x' and xuse eq 0 then begin  ;xuse counts how many times x was pressed. This forbids any response to it if x>0
    	modimage=oldmodimage
    	xuse=1
    	jroi=jroi-1
    	print,'   canceled'
    endif

endrep until done eq 1
wait,1.0
print,'Computing, PLEASE WAIT (could be several minutes)...'


;  Display last frame with ALL of the ROI's in it

window,1,xsize=xsize2x,ysize=ysize2x,title='All ROI locations',xpos=0,ypos=600
imagedisp=BYTSCL(modimage[*,*,frame[jroi]],max=maxim,min=minim)
tv,imagedisp

for iroi=0,jroi-1 do xyouts,xpos[iroi]-box2x,ypos[iroi]-12-box2x,STRTRIM(iroi,2),color=255,/device,charsize=1
wait,0.5

if lockin eq 1 then modimage=BYTSCL(modimage,max=maxim,min=minim) ;  This "locks in" the contrast adjustmentments into the movie.
end

; **************************************************

pro mpegdirect,imagestack,parentname,nframes,xsize,ysize,istart,iend

;   This creates a mpeg movie from the input imagestack.  It is "direct" because it creates the movie from MPEG commands,
;    not from within XINTERANIMATE.  See subroutine moviesmpg or movies.pro for XINTERANIMATE-based mpeg making.

; Open an MPEG sequence:

dim=[xsize,ysize]
mpegID = MPEG_OPEN(dim)
singleimage=fltarr(xsize,ysize)

read,'SLOW playback by an integer factor of :  ',slow
slow=FIX(slow)
istart=FIX(istart)
iend=FIX(iend)

for iframe=istart,iend do begin
    print,'Working on mpeg frame # '+STRTRIM(iframe)
    singleimage[*,*]=REVERSE(imagestack[*,*,iframe],2)
    for islow=0,slow-1 do MPEG_PUT, mpegID, IMAGE=singleimage, FRAME=(iframe-istart)*slow +islow
endfor

MPEG_SAVE, mpegID, FILENAME=parentname+' movie.mpg'

MPEG_CLOSE, mpegID

end

