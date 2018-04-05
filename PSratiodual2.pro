 pro PSratiodual2

;   This program is modifies from psratio2color.pro to optionally read DualView frames, the left side being alternating 514p and 514s
;      and the right side always being 442.  The revision notes start with those in psratio2color.
;   6-26-06 Written by Dan for "classic" excitation-based polarized TIRF.
;   10-1-07 Renamed from PSratio.pro for P-pol plus S-pol plus a third color, with frames in triplet sets

;   This program:

;   (a)   imports full images of sequentially one p-excited 514p, one s-excited 514s, one different color 442
;           OR alternating DualView images of 514p/442 and 514s/442
;   (b)   In the DualView case, the images are separated into 442, and 514s and 514p, the latter two being padded to
;           have the same # of frames as the 442 stack
;   (c)   takes a ratio of the two images, pixel by pixel
;   (d)   displays the ratio image and records a 16-bit tif image file of the ratio.
;   (e)   takes ratios for each pair of images in a time-sequence stack and makes tif files for all, for later animation in ImageJ

;  6-07-06  modified to flip IDL display and TIF display so they agree with CamWare:  90 degree CW rotation from actual stage
;  8-29-06  Modified to handle stacks correctly.
;  9-1-06   Error corrected in offcell background subtraction
;  9-2-06   S-pol properly defined
;  9-10-06  Preview of frames 1 and 1 added to aid in determining which one is s or p
;  10-4-06  Modified to also read in two single frames, one S and one P, with different filenames.
;           Be sure the difference is NOT in the last four characters just to the left of the .b16 because those are truncated off
;           when determining whether more than one file has the same parentname.
;  11-2-06  Option to read tif files from Andor camera.
;  11-7-06  Option to normalize PS ratio to the PS ratio of Rhodamine 6G recorded on the same day.
;  11-9-06  Improve thresholding display section so it starts with a P/S range of 0 to 2.
;  11-9-06  Repaired bug in fixed background subtraction option
;  11-9-06  Made all orientations consistent:
;           With Andor software set to acquire with vertical and horizontal flipped, stage agrees with Andor.
;           Now this IDL program agrees with Andor, and its tif output agrees with the default Image J (i.e., NOT deliberately flipped)
;  11-14-06 Fixed bug in opentfiles: now reads in single files with correct orientation.
;  11-17-06 Fixed error in setting multiplicative factor for output ratio files

;  12-1-06  Fixed error which inadvertantly set normRP equal to normRS.
;  12-22-06 Fixed reporting of nframes in opentfiles
;  12-31-06 Changed SUM from P+S to P+2S
;  10-1-07  Adapt read-in to a sequence of three exposures (514s, 514p, 442).
;  5-21-08  Normalization modified:  if rhod images are differently sized, use average over image.  Otherwise, do pixel-by-pixel norm.
;  12-20-08  Fixed bug in non-normalization mode
;  On the date that the camera was switched,the blockage method for identifying early polarization was also switched.
;  When the Sensicam was used, early frames had s-pol blocked.
;  When the Andor was installed, early frames had p-pol blocked.
;  Therefore, the filetype variable can be used as a guide for both camera and polarization situations.
;  1-22-09  Introduce option for DualView of alternating cer/514p and cer/514s

;  Set filetype = 1 for Sensicam .b16.  Set filetype = 2 for Andor .tif.
;  The Andor file system concatenates all the frames.  We first must separate them into a stack with ImageJ

; 1-25-09   Revised for optional DualView use (dual=1).  Also should work for single frames, with minor bugs corrected (dual=0)
;           Normalization can be done in either Dual View or in full frame, regardless of other format choices.
;           This program completely supercedes psratio2color.pro
; 2-5-09    Bug fixed in separation of p and s frames for dual=1 and type=-1 (sequential files).  Intermediate frames are now repeats.
; 2-09-09   DualView mode now allows alignment with dualalign subroutine, for S and HeCd; same alignment then used for P and HeCd
;  9-27-09  Allow correction for spillover diI excitation at 442, based on amount seen at 514
;  9-28-09  Background subtraction for diI spillover correction added.
; 10-23-13  Allow multi-image tifs to be read in (by setting multimage=1)
; 10-24-13	Bugs in definition of normRP and normRS in certain modes corrected
; 1-25-15   Allow frame forward and reverse in P/S display for setting threshold, brightness, contrast.
; 1-31-15	Output P/S tif stack without threshold modification, in addition to all other oututs. Call it Rfloatnon
;           Eliminate Sensicam camera option.  Change "HeCd" to Ceru".  Previously, Pfloat was used for P/S ratio to save memory.
;           Here, call it Rfloat.  Rfloat has thresholding in it.
;  2-2-15	Option for writing multi-image files
;  2-7-15   Fixed bug in requesting rhodamine normalization

ifirst=0 ;  In previous protocols, the first few frames were blocked.  But here the first good frame IS frame 0
again=1 ;   Again=1 offers a repeat with the same normalization values
dual=0  ;   dual=1 for alternating DualView images 514p/442 and 514s/442 doublets; dual=0 for sequential full frame 442,514p,514s triplets.
correct=0;  1 turns on the diI correction procedure; 0 leaves it off
multread=1   ;  Set =1 to read a multi-image input tif file; =0 to read a stack of single image tif files
multwrite=1  ;  Set =1 to write a multi-image output tif file;  =0 to output a stack of single image tif files.
newalign=1 ;set newalign=1 for experiments (typically beads) to determine a new set of Dual View left/right shift parameters
           ;set newalign=0 to use old typed-in results from prior (bead) experiments on actual cell experiments
xshift=-5  ; Worry about these shift parameters ONLY if newalign=0
yshift=-22
irepeat=0  ;  irepeat keeps track of the repeat #, starting with 1

startover:
irepeat=irepeat+1

print,'SELECT tif file with cell images'
if multread eq 0 then opentfiles,parentname,file,nframestot,imageblock,xsize,ysize
if multread eq 1 then opentmulti,parentname,file,nframestot,imageblock,xsize,ysize
if correct eq 1 then begin
      print,'SELECT tif file with diI spillover CORRECTION images'
      opentfiles,parentnamecorr,filecorr,nframestotcorr,imageblockcorr,xsizecorr,ysizecorr
      read,'442  background for diI correction = ',backcorr442
      read,'514s background for diI correction = ',backcorr514
endif
help,imageblock

;  We assume that the 442 (usually it is Ceru) image is s-pol
;  so we will compare it with the first s-pol 514 image in the stack of correction images
;  The 442 is the first image (frame=0) and the s-pol 514 is frame=1

wrongsize=0
if correct eq 1 then begin
    if xsizecorr ne xsize or ysizecorr ne ysize then begin
       print,'Pure diI image dimensions MUST equal data image dimensions for pixel-by-pixel correction.'
       print,'This is not the case here, so the correction is a single scalar factor over the whole image.'
       wrongsize=1
    endif
    diI442=fltarr(xsizecorr,ysizecorr)
    diI514=fltarr(xsizecorr,ysizecorr)
    diI442[*,*]=imageblockcorr[*,*,0]-backcorr442
    diI514[*,*]=imageblockcorr[*,*,1]-backcorr514
    ratiocorr=diI442/diI514
    ratiocorrmean=mean(ratiocorr)
    print,'Mean diI spillover factor = ',ratiocorrmean
endif

type=-1 ;  set this variable to -1.  It will remain -1 for multi-file stacks, and reset to 0 or 1 for single files
if nframestot eq 1 then begin  ; nframestot will equal 1 only if the selected file has a unique parentname.
;   Imageblock will still have three dimensions, but only containing a frame 0
    nframes=1
    read,'Does that single file contain S-pol (=0) or P-pol (=1) ? ',type
    spol=type
    print, 'Thanks for the info! Now read in the OTHER polarization single file'
    opentfiles,parentnameB,file,nframestot,imageblockB,xsize,ysize
    print,'Good job! Now read in the single file for the corresponding 442nm-excited image'
    opentfiles,parentnameC,file,nframestot,imageblockC,xsize,ysize
    ifirst=0
endif else begin  ; for nframestot gt 1

;  Are s-pol frames the odd or even frames?

;  This newly inserted section replaces the old way of identifying which image is p or s
;         Once the experiments become routine and are always the same, we can replace the choices with the truth and skip the questions.

;  Experiments do 514s,514p, and 442 nm, in some arrangement. Ask what is that sequence by displaying the first two or three frames

    singleframe0=intarr(xsize,ysize)
    singleframe1=intarr(xsize,ysize)
    singleframe2=intarr(xsize,ysize)

    singleframe0[*,*]=imageblock[*,*,0]
    singleframe1[*,*]=imageblock[*,*,1]
    singleframe2[*,*]=imageblock[*,*,2]
    window,0,xsize=xsize,ysize=ysize, title='Frame 0'
    tvscl,singleframe0
    window,1,xsize=xsize,ysize=ysize, title='Frame 1'
    tvscl,singleframe1
    if dual eq 0 then begin  ;  In the full frame triplet case, 442 is a separate frame
        window,2,xsize=xsize,ysize=ysize, title='Frame 2'
        tvscl,singleframe2
    endif
	wshow,0
; ppol=0 or 1 depending on whether it is in frame 0 or 1.  Likewise for spol
    read, 'Which Frame # contains the 514 nm P-pol? ',ppol
    if dual eq 1 then spol=1-ppol
    if dual eq 0 then read, 'Which Frame # contains the 514 nm S-pol? ',spol

; Ceru is the first frame # of the 442 nm-excited Ceru images
    if dual eq 0 then for i=0,2 do if spol ne i and ppol ne i then Ceru=i

endelse

;  NORMALIZATION OPTIONS

if again ne 1 or irepeat eq 1 then begin
    read,'NORMALIZE images against Rhodamine 6G (2=YES  DualView half frame R6G;  1=YES full frame R6G;  0=NO) ? ',norm

;   Read in Rhodamine 6G images to be used for normalization.5-21-08: Allow rhod images to be different size (xsizeR,ysizeR)
;   Normalization is pixel-by-pixel for same sized images.  Otherwise, it is by average over rhod images.

    if norm eq 2 then begin
        print,'Select DualView half frame P image of Rhodamine 6G'
        if filetype eq 1 then openbfiles,parentnamenormRP,filenormRP,nframesnorm,imageblockRP,xsizeR,ysizeR
        if filetype eq 2 then opentfiles,parentnamenormRP,filenormRP,nframesnorm,imageblockRP,xsizeR,ysizeR

        print,'Select DualView half frame S image of Rhodamine 6G'
        if filetype eq 1 then openbfiles,parentnamenormRS,filenormRS,nframesnorm,imageblockRS,xsizeR,ysizeR
        if filetype eq 2 then opentfiles,parentnamenormRS,filenormRS,nframesnorm,imageblockRS,xsizeR,ysizeR
        xsplitR=xsizeR/2
        normRS=fltarr(xsplitR,ysizeR)
        normRP=fltarr(xsplitR,ysizeR)
;  "Left is long" so 0:xsplitR contains the relavant R6G normalization
        normRS[*,*]=FLOAT(imageblockRS[0:xsplitR-1,*,0])
        normRP[*,*]=FLOAT(imageblockRP[0:xsplitR-1,*,0])

;  If image sizes are unequal, then use average values.

        if xsize ne xsizeR or ysize ne ysizeR then begin
            normRS[*,*]=mean(normRS)
            normRP[*,*]=mean(normRP)
        endif else begin

;  Possible smoothing of Rhodamine 6G images
            read,'SMOOTH Rhodamine 6G normalization images by n x n pixels where n = ', nsmooth
            normRS=smooth(normRS,nsmooth)
            normRP=smooth(normRP,nsmooth)
            window,3,title='P-pol Rhodamine 6G for normalization',xsize=xsplitR,ysize=ysizeR,xpos=525,ypos=0
            tvscl,normRP
            window,4,title='S-pol Rhodamine 6G for normalization',xsize=xsplitR,ysize=ysizeR,xpos=525,ypos=550
            tvscl,normRS
        endelse
    endif

    if norm eq 1 then begin
        print,'Select full frame P image of Rhodamine 6G.  For DualView mode, choose a D'
        opentfiles,parentnamenormRP,filenormRP,nframesnorm,imageblockRP,xsizeR,ysizeR

        print,'Select full frame S image of Rhodamine 6G'
        opentfiles,parentnamenormRS,filenormRS,nframesnorm,imageblockRS,xsizeR,ysizeR

        normRS=fltarr(xsizeR,ysizeR)
        normRP=fltarr(xsizeR,ysizeR)
        normRS[*,*]=FLOAT(imageblockRS[*,*,0])
        normRP[*,*]=FLOAT(imageblockRP[*,*,0])
		meannormRS=mean(normRS)
		meannormRP=mean(normRP)
;  If image sizes are unequal, then use average values.

        if xsize ne xsizeR or ysize ne ysizeR then begin
         	normRS=fltarr(xsize,ysize)
        	normRP=fltarr(xsize,ysize)
            normRS[*,*]=meannormRS
            normRP[*,*]=meannormRP
        endif else begin

;  Possible smoothing of Rhodamine 6G images
            read,'SMOOTH Rhodamine 6G normalization images by n x n pixels where n = ', nsmooth
            normRS=smooth(normRS,nsmooth)
            normRP=smooth(normRP,nsmooth)
            window,3,title='P-pol Rhodamine 6G for normalization',xsize=xsizeR,ysize=ysizeR,xpos=525,ypos=0
            tvscl,normRP
            window,4,title='S-pol Rhodamine 6G for normalization',xsize=xsizeR,ysize=ysizeR,xpos=525,ypos=550
            tvscl,normRS
        endelse
    endif

    if norm eq 0 then begin
        normRS=fltarr(xsize,ysize)+1.
        normRP=fltarr(xsize,ysize)+1.
    endif
endif

;  SEPARATE 442,514p, and 514s images

if dual eq 0 and type eq -1 then nframes=nframestot/3
if dual eq 1 and type eq -1 then begin
    nframes=nframestot  ; Ceru (i.e., 442) appears in every frame if dual=1
    nframes=2*(nframes/2)  ;  Convert nframes to an even #, discarding highest frame if odd
endif
print, 'nframes = ', nframes

if dual eq 1 and type eq -1 then xsize=xsize/2

Pfloat=fltarr(xsize,ysize,nframes)
Sfloat=fltarr(xsize,ysize,nframes)
Cerufl=fltarr(xsize,ysize,nframes)
Rfloat=fltarr(xsize,ysize,nframes)  ;  To be used for P/S ratio, thresholded.
RfloatNON=fltarr(xsize,ysize,nframes)  ;  To be used for P/S ratio, NONthresholded.

;  FULL frame triplets

if dual eq 0 and type eq -1 then begin ;  this is for multi-frame sequences of full frame triplets
    for iframe=0,nframes-1 do begin
        Pfloat[*,*,iframe]=FLOAT(imageblock[*,*,3*iframe+ppol])
        Sfloat[*,*,iframe]=FLOAT(imageblock[*,*,3*iframe+spol])
        Cerufl[*,*,iframe]=FLOAT(imageblock[*,*,3*iframe+Ceru])
    endfor
endif else begin
    if type eq 0 then begin;  this is for single frames that are S-pol first
       Sfloat=FLOAT(imageblock)
       Pfloat=FLOAT(imageblockB)
       Cerufl=FLOAT(imageblockC)
    endif
    if type eq 1 then begin;  this is for single frames that are P-pol first
       Pfloat=FLOAT(imageblock)
       Sfloat=FLOAT(imageblockB)
       Cerufl=FLOAT(imageblockC)
    endif
endelse

;  DualView half frame doublets

if dual eq 1 and type eq -1 then begin ;  this is for multi-frame sequences of DualView half frame doublets
    for iframe=0,nframes-1 do begin    ;  Pfloat and Sfloat are assigned repeats to fill out nframes
        Pfloat[*,*,iframe]=FLOAT(imageblock[0:xsize-1,*,2*(iframe/2)+ppol])
        Sfloat[*,*,iframe]=FLOAT(imageblock[0:xsize-1,*,2*(iframe/2)+spol])
        Cerufl[*,*,iframe]=FLOAT(imageblock[xsize:2*xsize-1,*,iframe])
    endfor
    ;  Do correction for diI spillover
    if correct eq 1 then if wrongsize eq 0 then Cerufl=Cerufl-Sfloat*ratiocorr else Cerufl=Cerufl-Sfloat*ratiocorrmean
endif else begin
    if type eq 0 then begin;  this is for single frames that are S-pol first
        Sfloat[*,*,0]=FLOAT(imageblock[0:xsize-1,*,0])
        Pfloat[*,*,0]=FLOAT(imageblockB[0:xsize-1,*,0])
        Cerufl[*,*,0]=FLOAT(imageblockC[xsize:2*xsize-1,*,0])
    endif
    if type eq 1 then begin;  this is for single frames that are P-pol first
        Pfloat[*,*,0]=FLOAT(imageblock[0:xsize-1,*,0])
        Sfloat[*,*,0]=FLOAT(imageblockB[0:xsize-1,*,0])
        Cerufl[*,*,0]=FLOAT(imageblockC[xsize:2*xsize-1,*,0])
    endif
endelse

;  ALIGN Ceru (right side of DualView) to match left (S or P). Best done on beads, but possible on granules
;    First, recombine (Sfloat,Cerufl) to make composite integer arrays that will be re-separated
;    by the dualalign subroutine as it is already written.   Re-enter dualalign again to perform same treatment
;    on (Pfloat,Cerufl).  Only the right side  (Cerufl in this case) gets shifted, rotated, and warped, so doing with
;    two calls to dualalign is redundant.  But this ensures that Pfloat gets smoothed and expanded similarly to Sfloat
;    and Cerufl.

if dual eq 1 then begin
    Pimageblock=FIX([Pfloat,Cerufl])
    Simageblock=FIX([Sfloat,Cerufl])
    rpt=0
    dualalign,parentname,file,nframes,Simageblock,2*xsize,ysize,xsize,Sfloat,Cerufl,xpand,newalign,xshift,yshift,$
       nx,ny,width,degtot,magtot,ntie,xout,yout,xin,yin,rpt
    rpt=1
    dualalign,parentname,file,nframes,Pimageblock,2*xsize,ysize,xsize,Pfloat,Cerufl,xpand,newalign,xshift,yshift,$
       nx,ny,width,degtot,magtot,ntie,xout,yout,xin,yin,rpt
endif

;  Background

quot0=fltarr(xsize,ysize)

print, 'Type -1 for GRAPHICAL background ROI selection.'
read,  'Type a nonnegative value for FIXED background : ',back

wdelete,0
wdelete,1
wdelete,2

window,0,xsize=xsize,ysize=ysize,title='S-POL image.  Mark OFF-CELL REGION'
tvscl,REFORM(Sfloat[*,*,0])
if back eq -1 then begin
    print,'Define OFF-CELL background, on S-POL image'
    offcell=defroi(xsize,ysize)

;  Subtract the off-cell background mean in each frame from the original floating point stacks

    for iframe=ifirst,nframes-1 do begin
        Pfloatslice=REFORM(Pfloat[*,*,iframe])
        Sfloatslice=REFORM(Sfloat[*,*,iframe])
        Pfloat[*,*,iframe]=(Pfloat[*,*,iframe]-mean(Pfloatslice[offcell]))/normRP[*,*]
        Sfloat[*,*,iframe]=(Sfloat[*,*,iframe]-mean(Sfloatslice[offcell]))/normRS[*,*]
    endfor
endif else begin

    for iframe=ifirst,nframes-1 do begin
        Pfloat[*,*,iframe]=(Pfloat[*,*,iframe]-back)/normRP[*,*]
        Sfloat[*,*,iframe]=(Sfloat[*,*,iframe]-back)/normRS[*,*]
    endfor
endelse

; Calculate P/S and P+2S with threholding

maxS=max(Sfloat)

maximum=1.0  ;  This sets the starting maximum within range of the largest P/S ratio expected, for successful initial display
minimum=0.
nincmin=0
nincmax=0

alpha=0.2  ;  fraction of the way from off-cell to brightest pixel in S-pol image, for initial guess at display thresholding

;  Display quotient for first good frame, to set threshold values for display

repeat A=GET_KBRD(0) until A eq '' ;  this clears the keyboard

ifrme=ifirst
repeat begin

    Sthresh=alpha*maxS
    incmin=0.05*(maximum-minimum)*nincmin
    incmax=0.05*(maximum-minimum)*nincmax
    mindisplay=minimum+incmin
    maxdisplay=maximum+incmax
    if mindisplay lt 0. then mindisplay=0.
    if maxdisplay lt 0 then maxdisplay=0.
    ;  Set off-cell quot values to zero

    ;   Adjust min and max for display

    for ix=0,xsize-1 do for iy=0,ysize-1 do if Sfloat[ix,iy,ifrme] le Sthresh or Pfloat[ix,iy,ifrme] le 0. $
    then quot0[ix,iy]=0. else quot0[ix,iy]=Pfloat[ix,iy,ifrme]/Sfloat[ix,iy,ifrme]

;   Display ratio image for iframe=ifirst, with three adjustables:

;  (1) Background threshold:  This zeroes out very dim regions in the original image arrays, e.g., off-cell places.
;     This affects both the display and the actual ratio array quot.
;  (2) Contrast control LO: set value in quot that is displayed as zero brightness; this affects only the display, not the data
;  (3) Contrast control HI: set value in quot that is displayed as full brightness; this affects only the display, not the data

    window,2,xsize=xsize,ysize=ysize,title='P/S FRM # '+STRTRIM(ifrme,2)+' Use keys as below; END w/RT.CLICK'
    bytequot0=BYTSCL(quot0,min=minimum+incmin,max=maximum+incmax)
    tv,bytequot0

    print,'BACKROUND THRESHOLD CONTROL:'
    print,'FIRST, left-click on ratio image, then type one of the following:'
    print,''
    print,' h = HIGHER background threshold'
    print,' l = LOWER background threshold'
    print,''
    print,' b = BRIGHTER by lowering TOP of greyscale'
    print,' d = DARKER by raising TOP of greyscale'
    print,''
    print,' c = INCREASE CONTRAST by raising BOTTOM of greyscale'
    print,' z = DECREASE CONTRAST by lowering BOTTOM of greyscale'
    print,''
    print,' f = FORWARD frame'
    print,' r = REVERSE frame'
    print,''
    print,'SPACEBAR = PROCEED as is to reading actual ratio values'

    A=GET_KBRD(1)
    if A eq 'h' then if alpha le .98 then alpha=alpha+0.02 else alpha=1.
    if A eq 'l' then if alpha ge .02 then alpha=alpha-0.02 else alpha=0.
    if A eq 'c' then nincmin=nincmin+1
    if A eq 'b' then nincmax=nincmax-1
    if A eq 'z' then nincmin=nincmin-1
    if A eq 'd' then nincmax=nincmax+1
    if A eq 'f' and ifrme lt nframes then ifrme=ifrme+1
    if A eq 'r' and ifrme gt 0 then ifrme=ifrme-1

endrep until A eq ' '

print,'RIGHT-CLICK on ratio image to end reading pixels and proceed to making tif files'
rdpix,quot0

; Do ratios for the whole stack, with parameters as selected above

maxsum=max(Pfloat+2*Sfloat)

;   Write tif files for SUM of P and S
multsum=0.001
if maxsum le 6550000. then multsum=0.01
if maxsum le 655000. then multsum=0.1
if maxsum le 65500. then multsum=1.
if maxsum le 6550. then multsum=10.
if maxsum le 655. then multsum=100.
if maxsum le 65.5 then multsum=1000.
if maxsum le 6.55 then multsum=10000.
if maxsum le 0.655 then multsum=100000.

print,'The recorded 16-bit tif SUM files have pixel values multiplied by', multsum
print,'Now writing SUM TIF file #...'
for iframe=ifirst,nframes-1 do begin
    print,iframe
    if iframe ge 0 and iframe le 9 then zeroes='000'
    if iframe ge 10 and iframe le 99 then zeroes='00'
    if iframe ge 100 and iframe le 999 then zeroes='0'
    if iframe ge 1000 and iframe le 9999 then zeroes=''
    sumsingle=REFORM(Pfloat[*,*,iframe]+2*Sfloat[*,*,iframe])
    negatives=where(sumsingle lt 0.,count)
    if count ne 0 then sumsingle[negatives]=0.
    tifnamestack=parentname+'-P2Ssum'+zeroes+STRTRIM(iframe,2)+'.tif'
    sumsingle=REVERSE(sumsingle,2) ; needed to allow Image J default orientation to agree with IDL
    if multwrite eq 0 then WRITE_TIFF,tifnamestack,FIX(multsum*sumsingle),/short,orientation=1
    if multwrite eq 1 then begin
    	tifnamestack=parentname+'-P2Ssum.tif'
    	if iframe eq ifirst then WRITE_TIFF,tifnamestack,FIX(multsum*sumsingle),/short,orientation=1
    	if iframe gt ifirst then WRITE_TIFF,tifnamestack,FIX(multsum*sumsingle),/short,orientation=1,/append
    endif
endfor

;  Write tif files for RATIO of P and S.  Do this for both thresholded ratio (Rfloat) and NONthresholded (RfloatNON)

;  THRESHHOLDED P/S:

;for iframe=ifirst,nframes-1 do for ix=0,xsize-1 do for iy=0,ysize-1 do if Sfloat[ix,iy,iframe] le Sthresh or Pfloat[ix,iy,iframe] lt 0. $
    ;then Rfloat[ix,iy,iframe]=0. else Rfloat[ix,iy,iframe]=Pfloat[ix,iy,iframe]/Sfloat[ix,iy,iframe]

thr=where(Sfloat gt Sthresh and Pfloat gt 0.)
Rfloat[thr]=Pfloat[thr]/Sfloat[thr]     ;  where Sfloat is le to Sthresh or Pfloat is neg, Rfloat will remain at 0.

non=where(Sfloat gt 0. and Pfloat gt 0.)
RfloatNON[non]=Pfloat[non]/Sfloat[non]  ;  where Sfloat is zero or neg or Pfloat is neg, RfloatNON will remain at 0.

stackmean=mean(Rfloat[*,*,[ifirst,nframes-1]])
maxratio=max(Rfloat)

multratio=1.
if maxratio le 6550 then multratio=10.
if maxratio le 655. then multratio=100.
if maxratio le 65.5 then multratio=1000.
if maxratio le 6.55 then multratio=10000.
if maxratio le 0.655 then multratio=100000.

;  The meaningful maxratio is unlikely to be greater than 2 so multratio must be at least 10000.
;  But random pixel dropouts in the denominator or numerator brights might cause a few pixels to have huge ratios,
;  so this next step suppresses that.

if multratio lt 10000. then multratio=10000.

print,'The recorded 16-bit tif RATIO files have pixel values multiplied by', multratio
print,'Now writing RATIO TIF file #...'
for iframe=ifirst,nframes-1 do begin
    print,iframe
    singlemean=mean(Rfloat[*,*,iframe])
    if iframe ge 0 and iframe le 9 then zeroes='000'
    if iframe ge 10 and iframe le 99 then zeroes='00'
    if iframe ge 100 and iframe le 999 then zeroes='0'
    if iframe ge 1000 and iframe le 9999 then zeroes=''
    quotsingle=REFORM(Rfloat[*,*,iframe]);/(singlemean/stackmean); which can be uncommented for samples with high photobleaching
    tifnamestack=parentname+'-PSratioTHR'+zeroes+STRTRIM(iframe,2)+'.tif'
    quotsingle=REVERSE(quotsingle,2); needed to allow Image J default orientation to agree with IDL
    if multwrite eq 0 then WRITE_TIFF,tifnamestack,FIX(multratio*quotsingle),/short,orientation=1
    if multwrite eq 1 then begin
    	tifnamestack=parentname+'-PSratioTHR.tif'
    	if iframe eq ifirst then WRITE_TIFF,tifnamestack,FIX(multratio*quotsingle),/short,orientation=1
    	if iframe gt ifirst then WRITE_TIFF,tifnamestack,FIX(multratio*quotsingle),/short,orientation=1,/append
    endif
endfor

;  NONthreshholded P/S:

non=where(Sfloat gt 0. and Pfloat gt 0.)
RfloatNON[non]=Pfloat[non]/Sfloat[non]  ;  where Sfloat is zero or neg or Pfloat is neg, RfloatNON will remain at 0.
stackmeanNON=mean(RfloatNON[*,*,[ifirst,nframes-1]])
maxratioNON=max(Rfloat)

multratioNON=1.
if maxratioNON le 6550 then multratioNON=10.
if maxratioNON le 655. then multratioNON=100.
if maxratioNON le 65.5 then multratioNON=1000.
if maxratioNON le 6.55 then multratioNON=10000.
if maxratioNON le 0.655 then multratioNON=100000.

;  The meaningful maxratioNON is unlikely to be greater than 2 so multratioNON must be at least 10000.
;  But random pixel dropouts in the denominator or numerator brights might cause a few pixels to have huge ratios,
;  so this next step suppresses that.

if multratioNON lt 10000. then multratioNON=10000.

print,'The recorded 16-bit tif RATIO files have pixel values multiplied by', multratioNON
print,'Now writing RATIO TIF file #...'
for iframe=ifirst,nframes-1 do begin
    print,iframe
    singlemeanNON=mean(RfloatNON[*,*,iframe])
    if iframe ge 0 and iframe le 9 then zeroes='000'
    if iframe ge 10 and iframe le 99 then zeroes='00'
    if iframe ge 100 and iframe le 999 then zeroes='0'
    if iframe ge 1000 and iframe le 9999 then zeroes=''
    quotsingleNON=REFORM(RfloatNON[*,*,iframe]);/(singlemean/stackmean); which can be uncommented for samples with high photobleaching
    tifnamestackNON=parentname+'-PSratioNON'+zeroes+STRTRIM(iframe,2)+'.tif'
    quotsingleNON=REVERSE(quotsingleNON,2); needed to allow Image J default orientation to agree with IDL
    if multwrite eq 0 then WRITE_TIFF,tifnamestack,FIX(multratioNON*quotsingleNON),/short,orientation=1
    if multwrite eq 1 then begin
    	tifnamestack=parentname+'-PSratioNON.tif'
    	if iframe eq ifirst then WRITE_TIFF,tifnamestack,FIX(multratioNON*quotsingleNON),/short,orientation=1
    	if iframe gt ifirst then WRITE_TIFF,tifnamestack,FIX(multratioNON*quotsingleNON),/short,orientation=1,/append
    endif
endfor

print,'Now writing Ceru file #...'
for iframe=ifirst,nframes-1 do begin
    print,iframe
    singlemean=mean(Cerufl[*,*,iframe])
    if iframe ge 0 and iframe le 9 then zeroes='000'
    if iframe ge 10 and iframe le 99 then zeroes='00'
    if iframe ge 100 and iframe le 999 then zeroes='0'
    if iframe ge 1000 and iframe le 9999 then zeroes=''
    Cerusingle=REFORM(Cerufl[*,*,iframe]);/(singlemean/stackmean); which can be uncommented for samples with high photobleaching
    tifnamestack=parentname+'-Cerul'+zeroes+STRTRIM(iframe,2)+'.tif'
    Cerusingle=REVERSE(Cerusingle,2); needed to allow Image J default orientation to agree with IDL
    if multwrite eq 0 then WRITE_TIFF,tifnamestack,FIX(Cerusingle),/short,orientation=1
	if multwrite eq 1 then begin
	  	tifnamestack=parentname+'-Cerul.tif'
		if iframe eq ifirst then WRITE_TIFF,tifnamestack,FIX(Cerusingle),/short,orientation=1
		if iframe gt ifirst then WRITE_TIFF,tifnamestack,FIX(Cerusingle),/short,orientation=1,/append
	endif
endfor

read,'START AGAIN on new images with the SAME NORMALIZATION (1=YES,  0=NO) ?',again
if again eq 1 then goto,startover
print, 'ALL DONE ! Program now ended'

end

;***************************************************************************************************
pro dualalign,parentname,file,nframes,imagedual,xsize,ysize,xSplit,left,right,xpand,newalign,xshift,yshift,$
    nx,ny,width,degtot,magtot,ntie,xout,yout,xin,yin,rpt

;  imagedual is the integer input stack; left and right are the floating point outputs, separated and shifted.

;created 6-21-05 by miriam
;modifed 12-13-05
;modified by Dan from simpleSplit.pro 5-9-06 to test various alignment schemes
;modiified 11-14-06 to read in a normalization file, which should be installed in the root folder c:\

;Take a dual cam image and split the image in a certain way-
;   left will be 0:678 (0:1376/2-1-xcut, where xcut=9)
;   right will be 697:1395 (688+xcut:1375)

;  10-28-08  modified from version appearing in xyratio.pro,
;            to read in files in main program, not here, and to delete options for Sensicam readin
;            Directly imported from FRET.pro
;  2-6-09    variable names generalized to left and right from long442 and short442.  For left and right being two colors,
;            "left is long" wavelength.  For left and right being two polarizations, "left is y" polarization.
;            "new" parameter determines whether:
;            (newalign=1) to find new alignment parameters or
;            (newalign=0) rely on old established typed-in shift parameters with no rot,mag, or warp
;            For newalign=1, the green/red superposition method is used.
;            For newalign=0, parameters typed in to the top of the main program - xshift and yshift - are used
;  2-09-09   Allow successive reuse without re-deriving all the shift, rotate, mag, and warp parameters
;             by passing them all through the calling list and setting rpt=1 in the main program


device,get_screen_size=screen
loadct,0

help,nframes

;  Imagedual can be a multiframe stack
help,imagedual
if nframes eq 1 then imagedual=reform(imagedual); This removes dimensions of size=1
;                               i.e., if only one frame is available this converts the stack to a 2D array.

print,file

xcut=0
xFull=xsize
xSplit=xFull/2-xcut
xpointsleft=indgen(xSplit)
xpointsright=indgen(xSplit)+xFull/2+xcut

if rpt eq 1 then goto,doagain  ;  rpt=1 is used for the second of successive visits to this subroutine

;  Imageleft and imageright are split from the original imagedual (for frame=0 only), but retain the full 12-bit integer range
;  Imageleftbyte and imagerightbyte compress these down to 8-bit integer scale
;  ImageRED and imageGREEN further compress these to 4-bit (only 16 values), indexed along the two axes of the color table

;  Image corresponding to left is to be displayed in 16-color red, which are indicies 0-15 in color table 42
;  So first we compress the range of imageleftbyte (from 0-255) to (0-15) by integer dividing by 16


imageleft=imagedual[xpointsleft,*,0]
imageright=imagedual[xpointsright,*,0]
window,0,title='imageleft',xsize=xSplit,ysize=ysize,xpos=400,ypos=0
tvscl,imageleft
imageleftnew=imageleft;  post-expansion work (even if xpand=1) on image CHx will be called "new"
imagerightnew=imageright

; Now give option for EXPANSION around a chosen point.  The expansion will look like a zoom but the overall image
;   size will be retained so it is effectively a cropping

print,'LEFT CLICK on image, then...'
read,'TYPE in EXPANSION factor (integer: =>1)',xpand
;xpand=1

;The new expanded image will be cropped to the original dimensions (xsplit,ysize) and stored in the cropped, expanded form.

if xpand gt 1 then begin
    print, 'Click cursor LEFT button for new image center'
    cursor,xcent,ycent,/up,/device
    print,'central position',xcent,ycent
    xcrop=xSplit/xpand  ;(xcrop,ycrop) is the size of the cropped image before expansion.
    ycrop=ysize/xpand

    xcroppts=indgen(xcrop)+xcent-xcrop/2; list of the x indicies in the cropped region, as seen from the original imageCHx
    ycroppts=indgen(ycrop)+ycent-ycrop/2; list of the y indicies in the cropped region

    imageleftcroptemp=imageleft[xcroppts,*]
    imagerightcroptemp=imageright[xcroppts,*]
    imageleftcrop=imageleftcroptemp[*,ycroppts]; imageCHxcrop are integer arrays of the cropped areas
    imagerightcrop=imagerightcroptemp[*,ycroppts]

    expand,imageleftcrop,xSplit,ysize,imageleftnew  ;  imageCHxnew is a cropped and expanded version of image CHx
    expand,imagerightcrop,xSplit,ysize,imagerightnew
endif

color42 ;  this generates and activates a new color table appropriate for "overlapping" red and green
;color43;   this generates and activates a new color table appropriate for "overlapping" red and blue

;  Relative shift determination

if newalign eq 1 then begin

    print,'BRIGHTNESS CONTROL for red/green composite image:'
    print,'FIRST, left-click on composite red/green image, then type one of the following:'
    print,''
    print,' b = BRIGHTER'
    print,' z = DIMMER'
    print,''
    print,'SPACEBAR = PROCEED as is'

    maxleft=max(imageleft)
    minleft=min(imageleft)
    maxright=max(imageright)
    minright=min(imageright)

    nbrighten=0
    repeat A=GET_KBRD(0) until A eq ''
    repeat begin

        imageleftbyte=BYTSCL(imageleftnew,max=(1.1^(-nbrighten))*maxleft,min=minleft)
        imageRED=imageleftbyte/16

;  Image corresponding to right is to be displayed in 16-color green, which are the indicies 16*n in color table 42 where n goes from 0 to 15
;  So we first compress the range of imagerightbyte (from 0-255) to (0-15) by integer dividing by 16 (throwing away the remainder)
;  Then we multiply by 16

        imagerightbyte=BYTSCL(imagerightnew,max=(1.1^(-nbrighten))*maxright,min=minright)
        imagegreen=(imagerightbyte/16)*16

        imageSUM=imageRED+imageGREEN

        window,2,xsize=xSplit,ysize=ysize,title='BRIGHTEN (b) or DIM (z)',xpos=400,ypos=0
        tv,imageSUM

        A=GET_KBRD(1)
        if A eq 'b' then nbrighten=nbrighten+1
        if A eq 'z' then nbrighten=nbrighten-1

    endrep until A eq ' '

    window,2,xsize=xSplit,ysize=ysize,title='TRANSLATE (r,l,u,d), ROTATE (c,x), MAGNIFY (m,s)',xpos=400,ypos=0
    tv,imageSUM

    print,'ALIGNMENT CONTROL to move GREEN image'
    print,'FIRST, left-click on composite red/green image, then type one of the following:'
    print,''
    print,' r = shift RIGHT one pixel'
    print,' l = shift LEFT one pixel'
    print,' u = shift UP one pixel'
    print,' d = shift DOWN one pixel'
    print,' c = rotate CLOCKWISE 0,1 deg around middle
    print,' x = rotate COUNTERclockwise 0.1 deg around middle'
    print,' m = MAGNIFY by 0.1% from middle'
    print,' s = SHRINK by 0.1% from middle'
    print,''

    print,' SPACEBAR = QUIT'

    imagenew=FLOAT(imagerightnew)

    nx=0
    ny=0
    nrot=0
    nmag=0

    deg=0.1
    mag=1.001

    repeat begin
        A=GET_KBRD(1)
        if A eq 'r' then begin
            nx=nx+1
            imagenew=shift(imagenew,1,0)
        endif
        if A eq 'l' then begin
            nx=nx-1
            imagenew=shift(imagenew,-1,0)
        endif
        if A eq 'u' then begin
            ny=ny+1
            imagenew=shift(imagenew,0,1)
        endif
        if A eq 'd' then begin
            ny=ny-1
            imagenew=shift(imagenew,0,-1)
        endif
        if A eq 'c' then begin
            nrot=nrot+1
            imagenew=rot(imagenew,deg,1.,cubic=-0.5)
        endif
        if A eq 'x' then begin
            nrot=nrot-1
            imagenew=rot(imagenew,-deg,1.,cubic=-0.5)
        endif
        if A eq 'm' then begin
            nmag=nmag+1
            imagenew=rot(imagenew,0.,mag,cubic=-0.5)
        endif
        if A eq 's' then begin
            nmag=nmag-1
            imagenew=rot(imagenew,0.,1./mag,cubic=-0.5)
        endif

        imageGREEN=(BYTSCL(imagenew,max=(1.1^(-nbrighten))*maxright,min=minright)/16)*16
        imageSUM=imageRED+imageGREEN

        window,2,xsize=xSplit,ysize=ysize,title='TRANSLATE (r,l,u,d), ROTATE (c,x), MAGNIFY (m,s).  Hit SPACEBAR when done',xpos=400,ypos=0
        tv,imageSUM

        print,'ROTATION clockwise = ',deg*nrot
        print,'MAGNIFICATION factor = ',mag^nmag
        print,'X shift = ', nx
        print,'Y shift = ',ny
        print,''
        print,' r = shift RIGHT one pixel'
        print,' l = shift LEFT one pixel'
        print,' u = shift UP one pixel'
        print,' d = shift DOWN one pixel'
        print,' c = rotate CLOCKWISE 0,1 deg around middle
        print,' x = rotate COUNTERclockwise 0.1 deg around middle'
        print,' m = MAGNIFY by 0.1% from middle'
        print,' s = SHRINK by 0.1% from middle'
        print,''
        print,' SPACEBAR = QUIT'

    endrep until A eq ' '
    wshow,2
    window,2,xsize=xSplit,ysize=ysize,title='WARP by clicking tie points. End selection with SPACEBAR ',xpos=400,ypos=0
    tv,imageSUM
    wshow,2

;  Find "tie points" to account for distorted images.  We will correct or "warp" the green image to force- match the red image
;  at tie points and hopefully at intermediate points.

    ntie=0
    print,'Move cursor to composite red/green image'

    repeat begin
        A=GET_KBRD(0)
    endrep until A eq ''

    print,'WARP correction: LEFT click on speckles in composite image'

    again:print,'RIGHT click if done, or LEFT click on GREEN speckle # ',ntie+1
    ntie=ntie+1
    cursor,xi,yi,/device,/up ; these are the coordinates where the speckle appears before correction
    if !mouse.button eq 4 then goto,skipout
    print,'                        Left click on corresponding RED speckle #',ntie
    cursor,xo,yo,/device,/up  ; these are the coordinates where we want the speckle to appear after correction
    print,'speckle # ',STRTRIM(ntie,2),'   (',STRTRIM(xi,2),',',STRTRIM(yi,2),') changed to (',STRTRIM(xo,2),',',STRTRIM(yo,2),')'

    print,''
    print,'Choose at least three tie points for warping to work'

;  Concatenate the clicked points into vectors xin, xout, yin, and yout

    if ntie eq 1 then begin
        xin=xi
        xout=xo
        yin=yi
        yout=yo
    endif else begin
        xin=[xin,xi]
        xout=[xout,xo]
        yin=[yin,yi]
        yout=[yout,yo]
    endelse
    goto,again

;   There is a choice of /quintic, /extrapolate, or /tps in the following warp_tri command
;   Which one works best, with how many tie points needs to be tested on bead images.

    xshift=nx  ; these enable dualalign to output the derived newalign=1 shift parameters to the main program
    yshift=ny
    skipout:print,''

endif else begin; end of newalign=1 shift determination, beginning of newalign=0 shift determination
    nx=xshift
    ny=yshift
    nrot=0  ;  newalign=0 only does x and y shifts, no rotation, magnignification, or warping for now.
    nmag=0.
    ntie=0
endelse

;   Possible smoothing:
read,'SMOOTHING:  # of adjacent pixels to be averaged = ',sm
width=[sm,sm]

print,'Alignment parameters now set.  Nice job!'
if ntie ge 3 then imagenew=warp_tri(xout,yout,xin,yin,imagenew,/extrapolate) else goto,finale
imageGREEN=(BYTSCL(imagenew,max=(1.1^(-nbrighten))*maxright,min=minright)/16)*16
imagetest=imageRED+imageGREEN
window,3,xsize=xSplit,ysize=ysize,title='WARPED image. Program ended',xpos=400,ypos=0
tv,imagetest

finale: degtot=nrot*0.1

magtot=1.001^nmag

;  Do the alignment for ALL the layers in the stack, and export as floating point 3D arrays left and right

doagain: ;  come stright here if rpt=1
left=fltarr(xSplit,ysize,nframes)
right=fltarr(xSplit,ysize,nframes)

if xpand eq 1 then begin
    print,'Alignment on whole stack in progress, now at frame #...'
    for iframe=0,nframes-1 do begin
        print,iframe
        sliceleft=float(imagedual[xpointsleft,*,iframe]); sliceCHx for each frame corresponds to imageCHx for frame=0
        sliceright=float(imagedual[xpointsright,*,iframe])
        if width[0] gt 1 then begin
            sliceleft=SMOOTH(TEMPORARY(sliceleft),width);  sliceleft is a smoothed version replacing slice left
            sliceright=SMOOTH(TEMPORARY(sliceright),width)
        endif
        sliceright=shift(sliceright,nx,ny);  sliceright is the channel that gets all the shifting etc.
        sliceright=rot(sliceright,degtot,magtot,cubic=-0.5)
        if ntie ge 2 then sliceright=warp_tri(xout,yout,xin,yin,sliceright,/extrapolate)
        left[*,*,iframe]=sliceleft[*,*];  CHx is the multiframe stack, floating point
        right[*,*,iframe]=sliceright[*,*]
    endfor
endif

;  For expansion on the whole stack, we have two options, based on the alignment parameters found above:
;    (a)  crop and expand to the original size (xsplit,ysize), then do alignment.  This gives a final zoomed-in view.
;    (b)  Do NOT crop.  Expand to (xsplit*xpand,ysize*xpand), do alignment, then REBIN down to (xsplit,ysize).
;      This should give an alignment made accurate by the expansion, but the original field of view is preserved.

if xpand gt 1 then begin
    read,'CHOOSE:  0 = Original field of view;   1 = Zoomed-in view',zoom
    print,'Alignment on whole stack in progress, now at frame #...'
    xsplitxpand=xsplit*xpand
    ysizexpand=ysize*xpand
    ;leftxpand=fltarr(xsplitxpand,ysizexpand,nframes)
    ;rightxpand=fltarr(xsplitxpand,ysizexpand,nframes)

    for iframe=0,nframes-1 do begin
        print,iframe
        sliceleft=float(imagedual[xpointsleft,*,iframe])
        sliceright=float(imagedual[xpointsright,*,iframe])
        if sm gt 1 then begin
            sliceleft=SMOOTH(TEMPORARY(sliceleft),width);  sliceleft is a smoothed version replacing slice left
            sliceright=SMOOTH(TEMPORARY(sliceright),width)
        endif

        if zoom eq 1 then begin
            sliceleftcroptemp=sliceleft[xcroppts,*]
            slicerightcroptemp=sliceright[xcroppts,*]
            sliceleftcrop=sliceleftcroptemp[*,ycroppts]
            slicerightcrop=slicerightcroptemp[*,ycroppts];    sliceCHx is cropped but not yet expanded

            expand,sliceleftcrop,xSplit,ysize,sliceleft  ;  sliceCHx is a cropped and expanded version of sliceCHx
            expand,slicerightcrop,xSplit,ysize,sliceright

            sliceright=shift(sliceright,nx,ny)
            sliceright=rot(sliceright,degtot,magtot,cubic=-0.5)
            if ntie ge 2 then sliceright=warp_tri(xout,yout,xin,yin,sliceright,/extrapolate)
        endif
        if zoom eq 0 then begin
            expand,sliceleft,xsplitxpand,ysizexpand,sliceleftxpand
            expand,sliceright,xsplitxpand,ysizexpand,slicerightxpand

            slicerightxpand=shift(slicerightxpand,nx,ny)
            slicerightxpand=rot(slicerightxpand,degtot,magtot,cubic=-0.5)
            if ntie ge 2 then slicerightxpand=warp_tri(xout,yout,xin,yin,slicerightxpand,/extrapolate)
            sliceleft=REBIN(sliceleftxpand,xsplit,ysize)
            sliceright=REBIN(slicerightxpand,xsplit,ysize)
        endif
        left[*,*,iframe]=sliceleft[*,*]
        right[*,*,iframe]=sliceright[*,*]
   endfor
endif

end



; ***************************************************************************************************
pro opentfiles,parentname,file,nframes,imageblock,xsize,ysize

;  This program is adapted from openbfiles.pro to read Andor tif files, which have the same numbering scheme as Sensicam b16 files

;  modified 8-29-06 to obtain parentname by truncating off last 8 characters in filename for stacks, and correctly handle single
;     files without truncation.

file = DIALOG_PICKFILE(TITLE='Select files',PATH='E:\Dropbox\TIR\BDNF\BDNF-PHl1103-8RWH fixedBGD', $
       FILTER = '*.tif')

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

for iframe=0,9999 do begin
    if iframe lt 10 then insert='000'
    if iframe ge 10 and iframe lt 100 then insert='00'
    if iframe ge 100 and iframe lt 1000 then insert='0'
    if iframe ge 1000 and iframe lt 10000 then insert=''
    filenew = header + insert + strtrim(iframe, 2) + '.tif'
    print,filenew
    exist=query_tiff(filenew)
    if exist eq 0 then goto,getout
endfor

getout:if iframe eq 0 then begin
    single=1
    nframes=1
    header=STRMID(file,0,length-4)
    print,'single file header = ',header
endif else nframes=iframe


;Set up stack of images
imageblock=UINDGEN(xsize,ysize,nframes)
singleimage=UINDGEN(xsize,ysize)

;Place images in the stack
for iframe=0,nframes-1 do begin
    if iframe lt 10 then insert='000'
    if iframe ge 10 and iframe lt 100 then insert='00'
    if iframe ge 100 and iframe lt 1000 then insert='0'
    if iframe ge 1000 and iframe lt 10000 then insert=''
    filenew = header + insert + strtrim(iframe, 2) + '.tif'
    exist=query_tiff(filenew)
    if exist eq 1 then begin
        singleimage=read_tiff(filenew)
        imageBlock[*,*,iframe]=REVERSE(singleimage,2) ; an inversion of the y-axis to agree with Andor and TIF output; REVERSE could have been used
    endif else begin
        nframes=iframe+1
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
pro opentmulti,parentname,file,nframes,imageblock,xsize,ysize

;  This program is adapted from opentfiles.pro to read Andor multi-image tif files

file = DIALOG_PICKFILE(TITLE='Select files',PATH='E:\Dropbox\TIR\BDNF\BDNF-PHl1103-8RWH fixedBGD', $
       FILTER = '*.tif')

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

;*****************************************************************************************************
pro color42

;  This sets up a color table to allow two image files to be "superimposed"
;    i.e., added, where one is red, the other green, and overlap is yellow.

;  This next command makes it possible to use a 256-color table in "true-color" video driver mode
device,decomposed=0

cmax=255.
fblue=FLTARR(256)
fgreen=FLTARR(256)
fred=FLTARR(256)
blue=BYTARR(256)
green=BYTARR(256)
red=BYTARR(256)
slope=cmax/15.

for i=0,255 do begin
      blue[i]=0
      green[i]=FIX((i/16)*slope)
      red[i]=FIX((i mod 16)*slope)
endfor

;  Convert palette to integers

itab=42
name='redgreenyellow'
MODIFYCT,itab,name,red,green,blue
LOADCT,itab

end


; ***************************************************************************************************
pro color43

;  This sets up a color table to allow two image files to be "superimposed"
;    i.e., added, where one is red, the other blue, and overlap is yellow.

;  This next command makes it possible to use a 256-color table in "true-color" video driver mode
device,decomposed=0

cmax=255.
fblue=FLTARR(256)
fgreen=FLTARR(256)
fred=FLTARR(256)
blue=BYTARR(256)
green=BYTARR(256)
red=BYTARR(256)
slope=cmax/15.

for i=0,255 do begin
      green[i]=0
      blue[i]=FIX((i/16)*slope)
      red[i]=FIX((i mod 16)*slope)
endfor

;  Convert palette to integers

itab=43
name='redbluepurple'
MODIFYCT,itab,name,red,green,blue
LOADCT,itab

end

