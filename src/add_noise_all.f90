
! Program add_noise
!---------------------------------------------------------------------------

Module Main

integer                         :: m_size,n_size
parameter(m_size =4000, n_size = 4000)  ! postage stamp size
real*4, dimension(0:m_size-1,0:n_size-1) :: object
Integer, Dimension(1:2)             :: fpixels,lpixels
Integer, Dimension(1:2)    :: naxes 

character(len=2000) :: filename,outfile

end Module Main

!--------------------------------------------------------

Program add_noise
use main
implicit none

integer :: i,j,ifile
real*4, dimension(0:m_size-1,0:n_size-1) :: noise, obs_image
real*4 :: mean,sigma  ! mean and width of noise dispersion
integer :: iseed
real :: gasdev
external gasdev


mean = 0.0
sigma = 0.001!5.0 star
iseed = -15

fpixels = 1
lpixels = m_size


do ifile = 247,261!1,400

   if(ifile<10)then
      write(filename,101) ifile
101   format('/disk2/ps1/cech/GREAT10/g10_beta/psfs/g10_beta_psfs_',I1,'.fits')
      write(outfile,201) ifile
201   format('/disk2/ps1/cech/GREAT10/g10_beta/psfs/g10_beta_psfs_',I1,'_noise.fits')
   elseif(ifile<100)then
      write(filename,102) ifile
102   format('/disk2/ps1/cech/GREAT10/g10_beta/psfs/g10_beta_psfs_',I2,'.fits')
      write(outfile,202) ifile
202   format('/disk2/ps1/cech/GREAT10/g10_beta/psfs/g10_beta_psfs_',I2,'_noise.fits')
   else
      write(filename,103) ifile
103   format('/disk2/ps1/cech/GREAT10/g10_beta/psfs/g10_beta_psfs_',I3,'.fits')
      write(outfile,203) ifile
203   format('/disk2/ps1/cech/GREAT10/g10_beta/psfs/g10_beta_psfs_',I3,'_noise.fits')
   end if

   do i = 0,m_size-1
      do j = 0,m_size-1
         noise(i,j) = gasdev(iseed)*sigma
      end do
   end do

   call readobj
   
   obs_image = object + noise
   
   call writeimage(obs_image)
end do

end Program add_noise

!--------------------------------------------------------
Subroutine readobj

Use Main
Implicit none

integer                    :: i,j
integer             :: status,unit,readwrite,blocksize
integer             :: group,nfound,naxis
real                :: nullval
logical anyf
integer             :: bitpix
integer             :: xpp,ypp,maxp
integer, Dimension(1:2) :: incs
real, Dimension(0:m_size -1, 0:n_size - 1) :: array


!write(*,*) 'ffll',fpixels(1),fpixels(2),lpixels(1),lpixels(2)
!write(*,*) 'x* y*',Xp,Yp

! The status and file unit parameters initialized.

status = 0      

! Open the read only FITS file 
! Denoted by region, fieldname and camerano from do loop above
! Note blocksize is anobsolete ignorable parameter

readwrite=0

unit = 99
call ftopen(unit,filename,readwrite,blocksize,status)
!call ftnopn(unit,filename,readwrite,status)


! Determine the size of the image.

call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

! Check that it found both NAXIS1 and NAXIS2 keywords.

if (nfound .ne. 2)then
   write(*,*) filename
   print *,'READIMAGE failed to read the NAXISn keywords.',filename
   return
end if

! Initialize variables

group=1
nullval=-999.0
naxis = 2                 ! no of dimensions
incs(1) = 1
incs(2) = 1               ! sampling interval of fits file

! Image info
 
write(*,*) 'dimension size', naxes

!bitpix = -32
! use ftgsve
 
bitpix = -64
! use ftgsvd

! Reading subsection of fits file and creating a 
! (boxsize + m_size/2)^2 data array 

call ftgsvd(unit,group,naxis,naxes,fpixels,lpixels,incs, &
             nullval,array,anyf,status)

! Closing FITS file 
      
call ftflus(unit,status)
call ftclos(unit,status)

! writing real*4 array to NAG routine happy real*8 real_part

write(*,*) 'mm', minval(array), maxval(array)

object = array

! Check for any error, and if so print out error messages.
      
if (status .gt. 0)then
   write(*,*) filename,'error problem, status = ', status
end if  

end subroutine readobj

!-------------------------------------------------------

subroutine writeimage(outarray)
use main
implicit none

integer callno,status,inunit,blocksize,bitpix,naxis
integer i,j,group,arraysize,outunit,readwrite,nkeys,nspace
logical simple,extend
character*(100) record
real, Dimension(0:m_size-1, 0:n_size-1) :: outarray

status=0
simple = .true.
bitpix = -32
naxis = 2
naxes = m_size
extend = .true.

!  Delete the file if it already exists, so we can then recreate it.
call deletefile(outfile,status)


! file unit number

inunit = 88
outunit = 99

!  Create the new empty FITS file.  The blocksize parameter is a
!  historical artifact and the value is ignored by FITSIO.


blocksize=1
call ftinit(outunit,outfile,blocksize,status)

! Open input file to read header

!readwrite=0
!call ftopen(inunit,filename,readwrite,blocksize,status)

! Find number of header keywords

!call ftghsp(inunit,nkeys,nspace,status)      

! read header line from infile
! write header line into outfile

!do i = 1,nkeys
!   call ftgrec(inunit,i,record,status)
!   call ftprec(outunit,record,status)
!end do

call ftphpr(outunit,simple,bitpix,naxis,naxes,0,1,extend,status)

! close infile

!call ftclos(inunit, status)


!  Initialize parameters about the FITS image.
!  BITPIX = -32 means that the image pixels will consist of 32-bit
!  real*4.  The size of the image is given by the NAXES values. 
!  The EXTEND = TRUE parameter indicates that the FITS file
!  may contain extensions following the primary array.
!  BITPIX = -64 , double precision floating point.

simple=.true.
bitpix=-32
naxis=2
extend=.true.


!  Write the array to the FITS file.

group=1
naxes = m_size
arraysize = m_size 

call ftp2de (outunit,group,arraysize,naxes(1),naxes(2),outarray,status)

if (status .gt. 0)write(*,*) 'status = ',status

!  The FITS file must always be closed before exiting the program. 

call ftclos(outunit, status)

if (status .gt. 0)write(*,*) 'status = ',status

end subroutine writeimage


!--------------------------------------------------------

subroutine deletefile(filename,status)
        
!  A simple little routine to delete a FITS file

integer status,unit,blocksize
character*(*) filename
      
!  Simply return if status is greater than zero
if (status .gt. 0)return

!  Get an unused Logical Unit Number to use to open the FITS file
call ftgiou(unit,status)

!  Try to open the file, to see if it exists
call ftopen(unit,filename,1,blocksize,status)

if (status .eq. 0)then
!         file was opened;  so now delete it 
   call ftdelt(unit,status)
else if (status .eq. 103)then
   !         file doesn't exist, so just reset status to zero and clear errors
   status=0
   call ftcmsg
else
   !         there was some other error opening the file; delete the file anyway
   status=0
   call ftcmsg
   call ftdelt(unit,status)
end if

!  Free the unit number for later reuse
call ftfiou(unit, status)

end subroutine deletefile

!--------------------------------------------------------------------
FUNCTION gasdev(idum)
  INTEGER idum
  REAL gasdev
  !U    USES ran1
  INTEGER iset
  REAL fac,gset,rsq,v1,v2,ran1
  SAVE iset,gset
  DATA iset/0/
  if (iset.eq.0) then
1    v1=2.*ran1(idum)-1.
     v2=2.*ran1(idum)-1.
     rsq=v1**2+v2**2
     if(rsq.ge.1..or.rsq.eq.0.)goto 1
     fac=sqrt(-2.*log(rsq)/rsq)
     gset=v1*fac
     gasdev=v2*fac
     iset=1
  else
     gasdev=gset
     iset=0
  endif
  return
END FUNCTION gasdev

!------------------------------------------------------------------
! random number generator from num recipes
! seed with a negative value of idum
! after that use the same value of seed

FUNCTION ran1(idum)
implicit none

INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
REAL*4 ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773)
PARAMETER (IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(1:NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/

if (idum.le.0.or.iy.eq.0) then
   idum=max(-idum,1)
   do j=NTAB+8,1,-1
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      if (j.le.NTAB) iv(j)=idum
   end do
   iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return

end FUNCTION ran1

!-----------------------------------------------------------------
