
! Program to fit a polynomial psf model to the stars
! selected by findstars
! Uses PGPLOT to show diagnostic plots
! Output file contains the polynomial model
!
! Catherine Heymans KSBf90 release for GREAT08
! heymans@roe.ac.uk
!
! Please acknowledge Catherine Heymans and Ludovic Van Waerbeke
! and cite Heymans et al (Mon.Not.Roy.Astron.Soc. 368 (2006) 1323-1339)
! in any publications that make use of this code.
! 
!-------------------------------------------------------------
Module Main

integer:: imax
parameter (imax=1e4)

real*4, dimension(1:imax) :: mag,rg,fr,fwhm
real*4, dimension(1:imax) :: x,y,flux
integer, dimension(1:imax) :: ixchip,iychip

!fitsio parameters
Integer, Dimension(1:2)    :: naxes 
integer                    :: status,unit
integer                    :: group,naxis
real                       :: nullval
logical anyf
integer, Dimension(1:2)    :: incs

integer                         :: m_size,n_size
parameter(m_size = 256, n_size = 256)  ! postage stamp size
!postage stamp array
real, dimension(0:m_size-1,0:n_size-1) :: object

!numrec fitting
integer                       :: nparams
parameter (nparams = 10)
integer, dimension(1:nparams) :: ia
real*4, dimension(1:nparams) :: pfit1,pfit2

! data
integer :: nstars
integer :: irgmax
parameter (irgmax = 100)
real*4, dimension(1:imax,1:4,0:irgmax) :: estar,ecor,pstar,psmstar
real*4, dimension(0:irgmax) :: rgstar,e1ave,e2ave,ec1ave,ec2ave
real*4, dimension(0:irgmax) :: e1ave_err,e2ave_err,ec1ave_err,ec2ave_err
real*4, dimension(0:irgmax) ::shsm

! parameter columns
integer :: ix,iy,imag,ifr,ifwhm,iflux,intot
integer :: rgmax,order,nchipx,nchipy
real*4 :: back

end Module Main
!-------------------------------------------------------------

Program psffit
use main
implicit none

integer:: i,j,k,ii,jj,irg,im, nprev,nim

character*150 arg,opt
integer iargc,narg


! psf fit params
integer :: nsc
real*4 :: xccd,yccd

! in/out files
character*250 filein, fileout,filecrit, filefits, fileparam,imfile,catfile
character*1500 line

! data
real*4, dimension(1:1000)  :: catalogue
real*4 ::xc,yc,xedge,yedge
integer, dimension(1:imax,0:irgmax) :: flag
real*4:: frl,frh,fwhml,fwhmh,magl,magh,frmax
real*4 :: rgbin

! getshapes output
real*4                  :: e(0:1)
real*4                  :: psm(0:1,0:1)
real*4                  :: psh(0:1,0:1)

! psf measurements
real*4 :: det, PsmInv(1:4)
real*4, dimension(1:imax,0:irgmax) :: shsmfact
real*4, dimension(1:imax) :: p1,p2
integer, dimension(1:imax) :: fitflag,chipflag

!postage stamp array
real, dimension(:,:), allocatable :: obj
integer, dimension(1:2) :: fpixels,lpixels
integer :: mpix,npix

character*50 :: toplabel
character*500 :: plot

!--RRG Additions--!
real*4::PSF_q(2,2), PSF_q4(2,2,2,2), RRG_Q(2,2)

! set up default plotting
plot = '/xwin'

narg=iargc()

do i=1,narg
   call getarg(i,opt)
   call getarg(i+1,arg)
   select case (opt)
   case ('-h')
      print *,'psffit.a'
      print *,'Example: '
      print *,'psffit.a -image image.fits -in SEx.cat -crit star_criteria.dat -param KSBf90.param -pgopen /xwin -out PSF.dat'
      print *,''
      print *,'       Defaults: '
      print *,'       -pgopen /xwin  Choose the output of the plots  '
      print *,'       -nim 1  Single image is the default            '
      print *,' If nim>1 then -image and -in should point to lists   '
      STOP

   case ('-in')
      catfile=arg
   case ('-nimage')
      read(arg,*) nim
   case ('-image')
      imfile=arg
   case ('-crit')
      filecrit=arg
   case ('-param')
      fileparam=arg
   case ('-out')
      fileout=arg
   case('-pgopen')
      plot = arg
   end select

enddo

call read_param(fileparam)

! set up the PSF polynomial order to fit

if(order==3)then
   ia =(/1,1,1,1,1,1,1,1,1,1/)   ! 3rd order
elseif(order==0)then
   ia=(/1,0,0,0,0,0,0,0,0,0/)   ! 0th order for testing
elseif(order==1)then
   ia =(/1,1,1,0,0,0,0,0,0,0/)  ! 1st order 
elseif(order==2)then
   ia =(/1,1,1,1,1,1,0,0,0,0/)  ! 2nd order 
else
   write(*,*) 'order restricted to 0,1,2,3'
   stop
end if

! open pgplot window
call pgopen(plot)
call pgsch(1.8)
call pgslw(2)
call pgsubp(2,1)

open(45,file=fileout,status='unknown')
open(33,file=filecrit,status='old')
!read(33,*) nstars,frl,frh,fwhml,fwhmh,magl,magh,frmax
read(33,*) magl,magh,frl,frh,fwhml,fwhmh!,frmax
! set binning scale for the varying weight function
!rgmax set in KSBf90.param
rgbin = frh/rgmax

write(*,*) 'NIM', nim

if(nim==1)then  ! command line files link directly to catalogue and image
   filein = catfile
   filefits = imfile
else ! command line links to files which contain links to the cats and images
   open(66,file=catfile,status='old')
   open(77,file=imfile,status='old')
end if

j = 0  ! star counter - zeroed outwith image loop
nprev = 0  ! cumulative number of stars in previous images

do im = 1,nim

   if(nim>1)then
      read(66,'(a)') filein
      read(77,'(a)') filefits
   end if
   
   open(44,file=filein,status='old')

   ! read in the stars based on the stellar selection
   ! criteria from findstars
   
   do while(j.le.imax)
      read(44,'(a)',end=23) line
      if(line(1:1).ne.'#')then
         read(line(1:1500),*) (catalogue(k),k=1,intot)
         
         write(*,*) catalogue(imag), ifwhm, catalogue(ifwhm)
                  
         if(catalogue(imag)< magh.and.catalogue(imag)> magl.and.&
              catalogue(ifwhm)< fwhmh.and.catalogue(ifwhm)> fwhml.and.&
              catalogue(ifr)< frh.and.catalogue(ifr)> frl)then
            
            j = j + 1
            x(j) = catalogue(ix)
            y(j) = catalogue(iy)
            mag(j) =  catalogue(imag)
            fr(j) = catalogue(ifr)
            fwhm(j) = catalogue(ifwhm)
            flux(j) = catalogue(iflux)
            
         end if
      end if
   end do
   write(*,*) 'WARNING: star list longer than imax in psffit.f90'
   write(*,*) 'WARNING: only using the first 1000 stars - modify imax in code'

23 continue
   nstars = j
   write(*,*) nstars-nprev,' stars selected in image', im

   ! open the image files
   call openfits(filefits)
   
   !calculate CCD chip widths if image is made of multiple CCD chips
   xccd = naxes(1)*1.0/nchipx
   yccd = naxes(2)*1.0/nchipy

   do i = nprev+1,nstars
      
      ! set CCD chip
      ixchip(i) = int(x(i)/xccd) + 1
      iychip(i) = int(y(i)/yccd) + 1
      
      ! select postage stamp from fits image
      ! deal with edges
      
      fpixels(1) = max(1,nint(x(i) - m_size/2)) 
      fpixels(2) = max(1,nint(y(i) - n_size/2)) 
      lpixels(1) = min(nint(x(i) + m_size/2 - 1),naxes(1))
      lpixels(2) = min(nint(y(i) + n_size/2 - 1),naxes(2))
      
      mpix = lpixels(1) - fpixels(1) + 1 
      npix = lpixels(2) - fpixels(2) + 1
      
      allocate (obj(0:mpix-1,0:npix-1))
      
      ! FITSIO routine to extract postage stamp of the star
      
      ! for bitpix = -32 use ftgsve
      ! for bitpix = -64 use ftgsvd 
      
      call ftgsve(unit,group,naxis,naxes,fpixels,lpixels,incs, &
           nullval,obj,anyf,status)
      
      ! Check for any error, and if so print out error messages.
      
      if (status .gt. 0)then
         write(*,*) filefits,'readobj error problem, fitsio status = ', status
         write(*,*) 'Does this file exist?'
         stop
      end if
      
      ! put all objects in a regular sized grid and subtract off 
      ! flat background noise set by back parameter in KSBparam.f90
      
      object = 0.0
      do ii = 0,mpix-1
         do jj = 0,npix-1
            object(ii,jj) = obj(ii,jj) - back   
         end do
      end do
      
      deallocate (obj)
      
      !   uncomment the following to see 2D images of the selected stars
      !   call plotarray(m_size,n_size,object)
      
      !   now loop through different weighting radii to see how the PSF
      !   ellipticity of object (:) varies with size
      
      do irg = 1,rgmax
         
         rgstar(irg) = frh + (irg-1)*rgbin
         rg(i) = rgstar(irg)
         
         xc = x(i) - fpixels(1) 
         yc = y(i) - fpixels(2) 
         xedge = min(xc, naxes(1)-x(i))
         yedge = min(yc, naxes(2)-y(i))
         
         ! getshape is the nuts and bolts of KSB
         
         PSF_Q =0.; PSF_Q4 = 0.
         call getshape(xc,yc,flux(i),rg(i),flag(i,irg),xedge,yedge,e,psm,psh, PSF_q, PSF_q4)
         !--Correct using RRG--!
         call RRG_PSF_Correction(rg(i), PSF_Q, PSF_Q4, RRG_Q)
         
         ! uncomment the following to see raw ellipticities plotted
         !call plotit(100000,x(i),y(i),e(0),e(1))
         
         det = 1.0/(psm(0,0)*psm(1,1) - psm(0,1)*psm(1,0))
         
         Psminv(1) = det * psm(1,1)
         Psminv(2) = -1.0 * det * psm(0,1)
         Psminv(3) = -1.0 * det * psm(1,0)
         Psminv(4) = det * psm(0,0)
         
         pstar(i,1,irg) = Psminv(1) * e(0) + Psminv(2) * e(1)
         pstar(i,2,irg) = Psminv(3) * e(0) + Psminv(4) * e(1)
         
         ! uncomment the following to see p's plotted
         !call plotit(100,x(i),y(i),pstar(i,1,irg),pstar(i,2,irg))
         
         estar(i,1,irg) = e(0)
         estar(i,2,irg) = e(1)
         psmstar(i,1,irg) = psm(0,0)
         psmstar(i,2,irg) = psm(1,0)
         psmstar(i,3,irg) = psm(0,1)
         psmstar(i,4,irg) = psm(1,1)
         
         shsmfact(i,irg) = (psh(0,0) + psh(1,1))/(psm(0,0) + psm(1,1))
         
      end do
   end do
   ! reset star counter and move on to next image
   nprev = nstars
end do

! psffit makes a polynomial fit to the PSF distribution
shsm = 20.0 ! give shsm an average value - properly set later
   
do irg = 1,rgmax
   
   call pgsci(1)
   call pgenv(1.0,naxes(1)*1.0,1.0,naxes(2)*1.0,1,0)
   
   if(irg.ne.0)then
      write(toplabel,102) frh + (irg-1)*rgbin
102   format('Weight radius', f7.2)
   else
      toplabel = 'radius = stellar flux radius'
   end if
   
   call pglab('x','y',toplabel)
   
   do i = 1,nstars
      p1(i) = pstar(i,1,irg)
      p2(i) = pstar(i,2,irg)
      
      if(estar(i,1,irg)**2 + estar(i,2,irg)**2 > 0.3)then
         fitflag(i) = 10
      else
         fitflag(i) = flag(i,irg)
      end if
   end do
   
   do ix = 1,nchipx
      do iy = 1,nchipy
         
         do i = 1,nstars
            if(ixchip(i)==ix.and.iychip(i)==iy)then
               chipflag(i) = fitflag(i)
            else
               chipflag(i) = 1
            end if
         end do
         
         ! uncomment to draw a histogram of p1 and p2
         !call pghist(nstars, p1, -0.1, 0.1, 20, 0)
         !call pgsci(2)
         !call pghist(nstars, p2, -0.1, 0.1, 20, 1)
         
         ! first remove outliers whose noise could bias the first fit
         
         do i = 1,4
            call sigma_clip(Nstars,6.0-i,p1,chipflag)
            call sigma_clip(Nstars,6.0-i,p2,chipflag)
         end do
         
         ! uncomment to see p1 and p2 as a function of stellar magnitude
         !call pgenv(magl,magh,-0.2,0.2,0,0)
         !call pgpt(nstars,mag,p1,4)
         !call pgsci(6)
         !call pgpt(nstars,mag,p2,4)
         
         ! uncomment to see selected e's/p's plotted
         nsc = 0
         do i = 1,nstars
            if(chipflag(i)==0)then
               !call plotit(2*naxes(1),x(i),y(i),estar(i,1,irg),estar(i,2,irg))
               call plotit(nint(naxes(1)/40.0),x(i),y(i),p1(i),p2(i))
               nsc = nsc + 1
            end if
         end do

         write(*,*) nsc, 'stars in chip', ix, iy

         ! fit a polynomial chip by chip
   
         call polyfit(p1,p2,chipflag)

!         call fitplot(nint(naxes(1)/40.0),chipflag)
         
         ! correct the PSF ellipticities with the above model

         call ecorrect(irg,chipflag)

         ! calculate average shsmfact(rg) in each CCD
         
         j = 0
         shsm(irg) = 0
         do i = 1,nstars
            if(chipflag(i)==0)then
               j = j + 1
               shsm(irg) = shsm(irg) +  shsmfact(i,irg)
            end if
            if(ixchip(i)==ix.and.iychip(i)==iy)then
               fitflag(i) = chipflag(i)
            end if
         end do
         
         shsm(irg) = shsm(irg)/(j*1.0)

         ! for each CCD chip write out the PSF parameters shsm and p

         write(45,108) irg,ix,iy,shsm(irg),&
                       (pfit1(j),j=1,nparams),(pfit2(j),j=1,nparams)
         
108      format(3I7, 21e15.5)

      end do
   end do
   call e_scatter_plot(irg,fitflag)
end do

! this plots how the PSF ellipticity varies with a different weighting
! radius rg

call plot_rg_variation

! Closing FITS file 
      
call ftflus(unit,status)
call ftclos(unit,status)

! Close pgplot window

call pgend
end Program psffit

!****************************************************************
! RRG PSF Correction model
!****************************************************************
!Uses 2nd and 4th order moments output by the KSB method [getshape()].
!Authored by CAJD on 13th November 2013.
!
!
!****************************************************************

subroutine RRG_PSF_Correction(rwindow, Q2, Q4, RRG_Q2)
  use Main
  !-Implements Equation (52) of Rhodes, Refridgier and Groth '10, to return the unweighted quadropole moment.
  !--Summation is assumed in the kk indices
  real*4::rwindow
  real*4::Q2(2,2), Q4(2,2,2,2)
  real*4::RRG_Q2(2,2)

  integer::ii,jj
  
  if(size(Q2) /= 4) then
     print *, size(Q2), size(Q2,1), size(q2,2)
     STOP 'RRG_PSF_Correction - FATAL ERROR - input Q2 is not of the correct size (2x2=4)'
  end if
  if(size(Q4) /= 16) STOP 'RRG_PSF_Correction - FATAL ERROR - input Q4 is not of the correct size (2x2x2x2=16)'

  do ii = 1, 2
     do jj = 1, 2
        RRG_Q2(ii,jj) = Q2(ii,jj) -(1.e0_4/(2.e0_4*rwindow*rwindow))*(Q2(ii,jj)*(Q2(1,1)+Q2(2,2)) - (Q4(ii,jj,1,1) + Q4(ii,jj,2,2)))
     end do
  end do

end subroutine RRG_PSF_Correction


!****************************************************************
! KSB nuts and bolts
!****************************************************************
!
! KSBf90 version of getshape
! With thanks to Nick Kaiser for the original imcat C version 
! of the getshape subroutine.
!
!Edited on the 12th Nov 2013 by Christopher Duncan, to add the PSF measurement Eqaution (50) of RRG. 
!****************************************************************

Subroutine getshape(xs,ys,fluxs,rwindow,gsflag,xedge,yedge,e,psm,psh, q, q4)
Use Main; use Moments, only:Weighted_Intensity_Moment_0, Weighted_Intensity_Moment_2, Weighted_Intensity_Moment_4
Implicit none

real*4         :: xs,ys,rwindow,fluxs
integer        :: istar,gsflag
real*4         :: xedge,yedge
real*4,intent(out)         :: q(1:2,1:2), q4(2,2,2,2)

real*4                        :: denom
integer                       :: i0, j0, i, j, di, dj, rmax, l, m
real*4                        :: r, dx, dy
real*4                        :: W, Wp, Wpp, fc, DD, DD1, DD2
real*4                        :: Xsm(0:1,0:1), Xsh(0:1,0:1)
real*4                        :: em(0:1), eh(0:1)
integer                       :: R_MAX_FACTOR
integer                       :: negflux
integer                       :: N1, N2

!--RRG Internal Declarations---!
integer::ii,jj,kk,ll
real*4::q0, q0_test, q_test(2,2)
real*4, allocatable, dimension(:)::obj_Grid_x, Obj_Grid_y, dx_test, dy_test
real*4,allocatable::W_test(:,:)
real*4::dxy(2)

real*4                  :: xim(0:1)
real*4                  :: d(0:1)
real*4                  :: e(0:1)
real*4                  :: psm(0:1,0:1)
real*4                  :: psh(0:1,0:1)
real*4                  :: fb0
real*4, dimension(0:1)     :: dfb

if(rwindow.le.0.0)then
   write(*,*) 'problem rwindow = 0 in getshapes'
end if

xim(0) = xs
xim(1) = ys

R_MAX_FACTOR = 4

N1 = n_size
N2 = m_size  

rmax = rwindow * R_MAX_FACTOR + 1       ! round up integer
i0 = xim(1)                          ! round down integer
j0 = xim(0)                          !  ""

! zero arrays

eh = 0.0
em = 0.0
e = 0.0
Xsh = 0.0
Xsm = 0.0
psm = 0.0
psh = 0.0
q = 0.0
d = 0.0
   
q0 = 0.0
q4 = 0.0

q_test = 0.e0_4; q0_test = 0.e0_4


!!$!--Implementation of CAJD methods for calculation of moments--!
!!$!-Unused as requires Intensity Map to be defined over a grid, howeve this grid could be easily populated-!
!!$!--Define Object Grid, used in my subroutines. To keep the Guassian Width defined in the same way as below, Grids should have dx = dy = 1--!
!!$!---NOTE, THIS IS NOT NORMALISED--!
!!$print *, 'Seting up grid and calling weights'
!!$allocate(Obj_Grid_x(size(Object,1))); allocate(Obj_Grid_y(size(Object,2))); Obj_Grid_x = 0.; Obj_Grid_y = 0.
!!$do i = 1, maxval( (/ size(Obj_Grid_x), size(Obj_Grid_y) /) )
!!$   if( i <= size(Obj_Grid_x) ) Obj_Grid_x(i) = i-1
!!$   if( i <= size(Obj_Grid_y) ) Obj_Grid_y(i) = i-1
!!$end do
!!$!--TESTING--!
!!$allocate(dx_test(size(Obj_Grid_x))); allocate(dy_test(size(Obj_Grid_y))); allocate(W_test(size(Obj_Grid_x),size(Obj_Grid_y)))
!!$dx_test = 0.; dy_test = 0.; W_test = 0.
!!$
!!$print *, 'Calling weighted intensity moments with width:', rwindow, xim
!!$call Weighted_Intensity_Moment_0(Intensity_Map = object, Grid_1 = Obj_Grid_x, Grid_2 = Obj_Grid_y, Centroid_Position = (/xim(0),xim(1)/), Gauss_Width = rwindow, I = q0)
!!$print *, 'Calling 2nd order:'
!!$call Weighted_Intensity_Moment_2(Intensity_Map = object, Grid_1 = Obj_Grid_x, Grid_2 = Obj_Grid_y, Centroid_Position = (/xim(0),xim(1)/), Gauss_Width =rwindow, IM = q, Normalisation = 1.e0_4, dx =  dx_test, dy =  dy_test, Wt = W_test)
!!$call Weighted_Intensity_Moment_4(Intensity_Map = transpose(object), Grid_1 = Obj_Grid_x, Grid_2 = Obj_Grid_y, Centroid_Position = (/xim(1),xim(0)/), Gauss_Width =rwindow, IM = q4, Normalisation = 1.e0_4)
!!$deallocate(Obj_Grid_x, Obj_Grid_y)


if(rmax>xedge.or.rmax>yedge)then
   gsflag = 1
   !go to 88
else
   do i = max(i0-rmax,0),min(i0+rmax,N2)
      do j = max(j0-rmax,0),min(j0+rmax,N1)
         
         di = i - i0
         dj = j - j0
         dx = j - xim(0)
         dy = i - xim(1)
         
         dxy = (/dx,dy/)

         r = sqrt(dx**2 + dy**2)
         
         if (r <= rmax)then
            W = exp(-0.5 * r**2 / rwindow**2)
            Wp = -0.5 * W / rwindow**2
            Wpp = 0.25 * W / rwindow**4
            
            fc = object(j,i)

!!$ DELETE         
!!$            !--TESTING--!
!!$            print *, 'fc*W:', W_test(j+1,j+1)/(fc*W)
!!$            print *, 'dx:', dx_test(j+1)/j; 
!!$            print *, 'dy:', dy_test(i+1)/j

   
!            if(gsflag.ge.2.and.object(j,i)<0.0)then
!               fc = (object(j-1,i-1) + object(j-1,i) + &
!                    object(j-1,i+1) + object(j,i-1) + &
!                    object(j,i+1)+ object(j+1,i-1)+ &
!                    object(j+1,i)+ object(j+1,i+1)) / 8.0
!            end if
            
            d(0) = d(0) + (W * fc * dx)
            d(1) = d(1) + (W * fc * dy)

            !--Set up quadropole moments-!
            q0 = q0 + (fc * W)
            do ii = 1, 2
               do jj = 1, 2
                  q(ii,jj) =  q(ii,jj) + (fc * W * dxy(ii) * dxy(jj))
                  do kk = 1,2
                     do ll = 1,2
                        q4(ii,jj,kk,ll) = q4(ii,jj,kk,ll) + (fc * W * dxy(ii) * dxy(jj) * dxy(kk) * dxy(ll))
                     end do
                  end do
               end do
            end do
!!$
!!$            q_test(1,1) = q_test(1,1) + (fc * W * dx * dx)
!!$            q_test(2,2) = q_test(2,2) + (fc * W * dy * dy)
!!$            q_test(1,2) = q_test(1,2) + (fc * W * dx * dy)
!!$            q_test(2,1) = q_test(2,1) + (fc * W * dx * dy)
!!$
!!$            !--Chris changes--!
!!$            q0_test = q0_test + (fc * W)
!!$            !-----------------!


            DD = di * di + dj * dj
            DD1 = dj * dj - di * di
            DD2 = 2 * di * dj
            Xsm(0,0) = Xsm(0,0) + &
                       (2 * W + 4 * Wp * DD + 2 * Wpp * DD1 * DD1) * fc 
            Xsm(1,1) = Xsm(1,1) + &
                       (2 * W + 4 * Wp * DD + 2 * Wpp * DD2 * DD2) * fc 
            Xsm(0,1) = Xsm(0,1) + (2 * Wpp * DD1 * DD2 * fc) 
            Xsm(1,0) = Xsm(1,0) + (2 * Wpp * DD1 * DD2 * fc)
            em(0) = em(0) + (4 * Wp + 2 * Wpp * DD) * DD1 * fc
            em(1) = em(1) + (4 * Wp + 2 * Wpp * DD) * DD2 * fc
            Xsh(0,0) = Xsh(0,0) + (2 * W * DD + 2 * Wp * DD1 * DD1) * fc
            Xsh(1,1) = Xsh(1,1) + (2 * W * DD + 2 * Wp * DD2 * DD2) * fc
            Xsh(0,1) = Xsh(0,1) + (2 * Wp * DD1 * DD2 * fc)
            Xsh(1,0) = Xsh(1,0) + (2 * Wp * DD1 * DD2 * fc)
            eh(0) = eh(0) + (2 * Wp * DD * DD1 * fc)
            eh(1) = eh(1) + (2 * Wp * DD * DD2 * fc)

         end if
      end do
   end do

!!$   !---TESTING DELETE---!
!!$   print *, 'q0:', q0/q0_test
!!$   print *, 'q:', q/q_test
!!$   read(*,*)

   !normalise d
   
   if (fluxs > 0)then 
      do l = 0,1
         d(l) = d(l)/fluxs
      end do
      negflux = 0
   else 
      negflux = 1
   end if
   
   !--RRG PSF MODEL--!
   

   
   ! calculate ellipticities
   denom = q(1,1) + q(2,2)                 
   
   if (denom > 0)then
      e(0) = (q(1,1) - q(2,2)) / denom
      e(1) = (q(1,2) + q(2,1)) / denom
      em = em / denom
      eh = eh / denom
      eh = eh + (2 * e)

      do l = 0,1
         do m = 0,1
            psm(l,m) = 0.5 * (Xsm(l,m) / denom - e(l) * em(m))
            psh(l,m) = Xsh(l,m) / denom - e(l) * eh(m)
         end do
      end do
   else
      gsflag = 4
   end if
end if

!--Renormalise--!
q = q/q0
q4 = q4/q0

!88 continue



!write(*,*) 'e',(e(i),i=0,1),gsflag
!write(*,*) (d(i),i=0,1)
!write(*,*) 'psm',(psm(i,0),i=0,1),(psm(i,1),i=0,1)
!write(*,*) 'psh',(psh(i,0),i=0,1),(psh(i,1),i=0,1)

End subroutine getshape

!-------------------------------------------------------------------------

Subroutine ecorrect(irg,rej)

Use Main
Implicit none

Integer :: i,j,irg
Real*4, dimension(1:2) :: pfitted
integer, dimension(1:imax) :: rej
Real*4, dimension(1:2) :: eopsm


do i = 1,nstars

   if(rej(i)==0)then
   
      call calcp(i,pfitted)
      
      eopsm(1) = Psmstar(i,1,irg)*pfitted(1)  + Psmstar(i,2,irg)*pfitted(2)
      eopsm(2) = Psmstar(i,3,irg)*pfitted(1)  + Psmstar(i,4,irg)*pfitted(2) 
      
      do j = 1,2
         ecor(i,j,irg) = estar(i,j,irg) - eopsm(j)
      end do
      
   end if
end do
   
End Subroutine ecorrect

!------------------------------------------------------
Subroutine e_scatter_plot(irg,rej)
use main
implicit none

Integer :: i,j,nd,irg
integer, dimension(1:imax) :: rej
Real*4, dimension(1:2) :: esum, esumsq, ecsum, ecsumsq
Real*4, dimension(1:2) :: e_mean, e_sd
Real*4, dimension(1:2) :: ec_mean, ec_sd
real*4                 :: fwhm_mean,fr_mean

nd = 0
esum = 0.0
esumsq = 0.0
ecsum = 0.0
ecsumsq = 0.0
fwhm_mean = 0.0

do i = 1,nstars

   if(rej(i)==0)then
      do j = 1,2
         ecsum(j) = ecsum(j) + ecor(i,j,irg)
         ecsumsq(j) = ecsumsq(j) + ecor(i,j,irg)**2
         esum(j) = esum(j) + estar(i,j,irg)
         esumsq(j) = esumsq(j) + estar(i,j,irg)**2
      end do
      
      nd = nd + 1
      fwhm_mean = fwhm_mean + fr(i)!FWHM(i)

   end if
end do
   
ec_mean = ecsum / (nd*1.0)
ec_sd = sqrt(ecsumsq/(nd*1.0) - ec_mean**2 )
e_mean = esum / (nd*1.0)
e_sd = sqrt(esumsq/(nd*1.0) - e_mean**2 )
fr_mean = fwhm_mean/nd

write(*,*) 'uncorrected',e_mean, e_sd
write(*,*) 'corrected',ec_mean, ec_sd
write(*,*) nd, e_mean, e_sd, fwhm_mean/nd 

! save values for estar(rg) test

e1ave(irg) = e_mean(1)
e1ave_err(irg) = e_sd(1)/sqrt(nd-1.0)
e2ave(irg) = e_mean(2)
e2ave_err(irg) = e_sd(2)/sqrt(nd-1.0)
ec1ave(irg) = ec_mean(1)
ec1ave_err(irg) = ec_sd(1)/sqrt(nd-1.0)
ec2ave(irg) = ec_mean(2)
ec2ave_err(irg) = ec_sd(2)/sqrt(nd-1.0)

! draw plots
call pgsci(1)
call pgenv(-0.1,0.1,-0.1,0.1,1,0)
call pglab('e1','e2','PSF ellipticity')
do i = 1,nstars
   if(rej(i)==0)then
      call pgpt(1,estar(i,1,irg),estar(i,2,irg),4)
   end if
end do
call pgsci(6)
do i = 1,nstars
   if(rej(i)==0)then
      call pgpt(1,ecor(i,1,irg),ecor(i,2,irg),3)
   end if
end do

End Subroutine e_scatter_plot

!********************************************************************
! Fitting and Statistics routines
!********************************************************************

Subroutine polyfit(p1,p2,flag)
use main
implicit none

integer :: i,gdstr
real*4, dimension(1:nstars) ::p1,p2
integer, dimension(1:nstars) :: flag

real*4, dimension(1:nparams,1:nparams) :: covar
real*4, dimension(1:nstars) :: xpos,ypos
real*4, dimension(1:nstars) :: tofit1,tofit2,sig
real*4 :: chisq

gdstr = 0
sig = 0.1
! set starting guess equal to zeroth order average
pfit1 = 0.0
pfit2 = 0.0
do i = 1,nstars
   if(flag(i)==0)then
      gdstr = gdstr + 1
      xpos(gdstr) = x(i)
      ypos(gdstr) = y(i)
      tofit1(gdstr) = p1(i)
      tofit2(gdstr) = p2(i)
      pfit1(1) = pfit1(1) + p1(i)
      pfit2(1) = pfit2(1) + p2(i)
   end if
end do

pfit1(1) = pfit1(1) / gdstr
pfit2(1) = pfit2(1) / gdstr

if(gdstr>10)then   

   call lfit(xpos,ypos,tofit1,sig,gdstr,pfit1,ia,nparams,covar,nparams,chisq)
   call lfit(xpos,ypos,tofit2,sig,gdstr,pfit2,ia,nparams,covar,nparams,chisq)
else
   write(*,*) 'Less than 10 stars in fit, reduce nchipx, nchipy'
end if

!write(*,*) 'p1', pfit1
!write(*,*) 'p2', pfit2

end Subroutine polyfit
!---------------------------------------------------------------

Subroutine calcp(i,pfitted)
Use Main
Implicit none

integer,intent(in)         :: i
integer :: j,k
real*4          :: afunc(1:nparams)
Real*4, dimension(1:2) ::pfitted

call funcs(x(i),y(i),afunc,nparams)

pfitted = 0.0

do j = 1,nparams
   pfitted(1) = pfitted(1) + ia(j)*pfit1(j)*afunc(j)
   pfitted(2) = pfitted(2) + ia(j)*pfit2(j)*afunc(j)
end do

end Subroutine calcp
!-----------------------------------------------------------------------

Subroutine sigma_clip(N,clip,cliparray,flag)
implicit none
integer :: i,ct,rct,N
real*4, dimension(1:N) :: cliparray
integer, dimension(1:N) :: flag
real*4 :: clipsum, sumosq
real*4 :: mean,sd,clip

! calculate mean and variance

clipsum = 0.0
ct = 0
rct = 0
sumosq = 0

do i = 1,N

   if(flag(i)==0)then
      ct = ct + 1
      clipsum = clipsum + cliparray(i)
      sumosq = sumosq + cliparray(i)**2
   end if

end do

mean = clipsum/(ct*1.0)
sd = sqrt( (sumosq/(ct*1.0)) - mean**2 )

! remove > clip*sigma outliers

do i = 1,N
   if(flag(i)==0)then
      if(abs(cliparray(i)-mean)>clip*sd)then
         flag(i)= flag(i) + 10
      else
         rct = rct + 1
      end if
   end if
end do

!if(clip==3.0)then
!   write(*,107) N, ct, rct
!107 format('Nstars selected = ', I4, ' Good shapes =', I4, ' After sigma clip', I4) 
!end if

end Subroutine sigma_clip

!********************************************************************
! PGPLOT plotting routines
!********************************************************************

Subroutine plot_rg_variation
use main
implicit none
integer :: i,j

call pgsubp(2,2)
call pgsci(1)
call pgenv(rgstar(1)-0.2,rgstar(rgmax)+0.2,minval(e1ave)-0.002,maxval(e1ave)+0.002,0,0)
call pglab('rg','<e1_star>','PSF ellipticity variation with rg')
call pgsls(4)
call pgmove(0.0,0.0)
call pgdraw(15.0,0.0)
call pgsls(1)
call pgline(rgmax,rgstar,e1ave)
call pgerry(rgmax,rgstar,e1ave+e1ave_err,e1ave-e1ave_err,1.0)
call pgsci(6)
call pgline(rgmax,rgstar+0.05,ec1ave)
call pgerry(rgmax,rgstar+0.05,ec1ave+ec1ave_err,ec1ave-ec1ave_err,1.0)

call pgsci(1)
call pgenv(rgstar(1)-0.2,rgstar(rgmax)+0.2,minval(e2ave)-0.002,maxval(e2ave)+0.002,0,0)
call pglab('rg','<e2_star>','')
call pgsls(4)
call pgmove(0.0,0.0)
call pgdraw(15.0,0.0)
call pgsls(1)
call pgline(rgmax,rgstar,e2ave)
call pgerry(rgmax,rgstar,e2ave+e2ave_err,e2ave-e2ave_err,1.0)
call pgsci(6)
call pgline(rgmax,rgstar+0.05,ec2ave)
call pgerry(rgmax,rgstar+0.05,ec2ave+ec2ave_err,ec2ave-ec2ave_err,1.0)

call pgsci(1)
call pgenv(rgstar(1)-0.2,rgstar(rgmax)+0.2,minval(shsm),maxval(shsm),0,0)
call pglab('rg','shsmfact','')
call pgline(rgmax,rgstar,shsm)


end Subroutine plot_rg_variation

!-----------------------------------------------------------------------

Subroutine plotarray(N1,N2,array)
implicit none

integer :: N1,N2
real*4, dimension(1:N1,1:N2) :: array
real*4, dimension(1:6) :: Tr

Tr = 0
Tr(2) = 1.0
Tr(6) = 1.0

call pgenv(1.0,N1*1.0,1.0,1.0*N2,1,0)
call pggray(array, N1, N2, 1, N1, 1, N2, minval(array), maxval(array), TR)

end Subroutine plotarray
!--------------------------------------------------------------------

Subroutine plotit(magfact,xstar,ystar,e1,e2)
Implicit none

real*8            :: eta,axis,theta,rad_to_deg
real*8, dimension(1:2)     :: thetasin,thetacos
real*4, dimension(1:2)     :: Xvect,Yvect
real*4 :: e1,e2,xstar,ystar
integer,intent(in)                     :: magfact
!parameter(magfact = 500)
integer                    :: i,j,k

rad_to_deg = 180.0/3.1415926
         
eta= sqrt(e1**2 + e2**2)
axis= sqrt((eta+ 1.0)/(1.0-eta))

if(eta.ne.0.0)then
   
   theta = acos(e1/eta) * rad_to_deg
            
   thetacos(1) = theta/2.0
   thetacos(2) = (360.0 - theta)/2.0
   
   theta= asin(e2/eta)* rad_to_deg
            
   if(theta>0)then
      thetasin(1) = theta/2.0
      thetasin(2) = (180.0 - theta)/2.0
   else
      thetasin(1) = (180.0 - theta)/2.0
      thetasin(2) = (360.0 + theta)/2.0
      
   end if
   
   do k = 1,2
      do j = 1,2
         if(abs(thetasin(k) - thetacos(j))<1e-8)then
            theta = thetasin(k)
         end if
      end do
   end do
   
else
   theta = 0.1
end if
         
Xvect(1) = xstar - (magfact*eta/2. * cos(theta/rad_to_deg))
Xvect(2) = xstar + (magfact*eta/2. * cos(theta/rad_to_deg))
Yvect(1) = ystar - (magfact*eta/2. * sin(theta/rad_to_deg))
Yvect(2) = ystar + (magfact*eta/2. * sin(theta/rad_to_deg))

call pgline(2,Xvect,Yvect)
Yvect(1) = ystar 
Xvect(1) = xstar
!call pgpt(1,Xvect,Yvect,4)
      
end subroutine plotit

!--------------------------------------------------------------------------

Subroutine fitplot(magfact,flag)
use main
Implicit none

integer         :: i,j,k,magfact
real*4          :: afunc(1:nparams)
integer, dimension(1:nstars) :: flag
real*4 :: p1,p2,xstar,ystar


do i = 1,nstars
   if(flag(i)==0)then

      xstar = x(i)
      ystar = y(i)

      call funcs(xstar,ystar,afunc,nparams)

      p1 = 0.0
      p2 = 0.0

      do j = 1,nparams
         p1 = p1 + ia(j)*pfit1(j)*afunc(j)
         p2 = p2 + ia(j)*pfit2(j)*afunc(j)
      end do

      call pgsci(6)
      call plotit(magfact,xstar,ystar,p1,p2)

   end if
end do

call pgsci(1)

end Subroutine fitplot

!**************************************************************
! Numerical Recipes
! You should have an onsite Numerical Recipes
! licence to use the following 
!**************************************************************

SUBROUTINE funcs(x1val,x2val,afunc,ma)
  integer :: ma
  real*4  :: x1val,x2val,afunc(1:ma)

  afunc(1)=1.
  afunc(2)=x1val
  afunc(3)=x2val
  afunc(4)=x1val**2
  afunc(5)=x1val*x2val
  afunc(6)=x2val**2
  afunc(7)=x1val**3
  afunc(8)=x2val**3
  afunc(9)=x1val**2 * x2val
  afunc(10)=x1val* x2val**2

END SUBROUTINE funcs

!-----------------------------------------------------------------------
      SUBROUTINE LFIT(Xa,Xb,Y,SIG,NDAT,A,IA,MA,COVAR,npc,CHISQ)
      integer ma,ia(ma),npc,ndat,mmax
      PARAMETER (MMAX=50)
      real*4 chisq,Xa(NDAT),Xb(NDAT),Y(NDAT),SIG(NDAT),A(MA)
      real*4 COVAR(npc,npc)

      integer  i,j,k,l,m,mfit
      real*4 sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)

      mfit = 0
      do j = 1,ma
         if(ia(j).ne.0) mfit = mfit + 1
      end do

      if(mfit.eq.0)then
         write(*,*) 'lfit: no parameters to be fitted'
         stop
      end if
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE


      DO I=1,ndat
        CALL FUNCS(Xa(I),Xb(I),AFUNC,MA)
        YM=Y(I)
        IF(MFIT.LT.MA) THEN
          DO J=1,MA
             if(ia(j).eq.0) ym = ym - a(j)*afunc(j)
          end do
        ENDIF

        SIG2I=1./SIG(I)**2
        j = 0

        do l = 1,ma
           if(ia(l).ne.0)then
              j = j + 1
              WT=AFUNC(l)*SIG2I

              k = 0
              do m = 1,l
                 IF(IA(M).NE.0)THEN
                    k = k + 1
                    covar(j,k) = covar(j,k) + wt*afunc(m)
                 end if
              end do

              BETA(J)=BETA(J)+YM*WT
           end if
        end do
      end do

      DO 21 J=2,MFIT
         DO 19 K=1,J-1
            COVAR(K,J)=COVAR(J,K)
 19      CONTINUE
 21   CONTINUE

      CALL GAUSSJ(COVAR,MFIT,npc,BETA,1,1)

      j = 0
      do l = 1,ma
         if(ia(l).ne.0)then
            j = j + 1
            a(l) = beta(j)
         end if
      end do

      CHISQ=0.
      DO 24 I=1,ndat
        CALL FUNCS(Xa(I),Xb(I),AFUNC,MA)
        SUM=0.
        DO 23 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
23      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
24    CONTINUE
      CALL COVSRT(COVAR,npc,MA,ia,MFIT)
      RETURN
      END

!------------------------------------------------------------------------
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      PARAMETER (NMAX=50)
      integer np,K
      real*4 A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)

      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                write(*,*) 'Singular matrix'
                stop
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) then
           write(*,*) 'Singular matrix'
           stop
        end IF
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,nint(INDXR(L)))
            A(K,nint(INDXR(L)))=A(K,nint(INDXC(L)))
            A(K,nint(INDXC(L)))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
!----------------------------------------------------------------------
      SUBROUTINE COVSRT(COVAR,npc,MA,ia,MFIT)
      real*4 COVAR(npc,npc)
      integer ma,ia(ma),mfit,npc
      integer  i,j,k
      real*4 swap

      do i = mfit + 1, ma
         do j = 1,i
            covar(i,j) = 0
            covar(j,i) = 0
         end do
      end do

      k = mfit
      do j = ma, 1, -1
         if(ia(j).ne.0)then
            do i = 1,ma
               swap = covar(i,k)
               covar(i,k) = covar(i,j)
               covar(i,j) = swap
            end do
            do i = 1,ma
               swap = covar(k,i)
               covar(k,i) = covar(j,i)
               covar(j,i) = swap
            end do
            k = k - 1
         end if
      end do
      
      RETURN
      END

!***************************************************************
! FITSIO ROUTINES 
!***************************************************************

Subroutine openfits(filename)
use main
implicit none

character(len=200) filename
integer             :: readwrite,blocksize
integer             :: nfound

! The status and file unit parameters initialized.

status = 0      

! Open the read only FITS file 
! Denoted by region, fieldname and camerano from do loop above
! Note blocksize is anobsolete ignorable parameter

readwrite=0

unit = 99
write(*,*) filename
call ftopen(unit,filename,readwrite,blocksize,status)

! Determine the size of the image.

call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

! Check that it found both NAXIS1 and NAXIS2 keywords.

if (nfound .ne. 2)then
   write(*,*) filename
   print *,'READIMAGE failed to read the fits file',filename
   return
end if

! Initialize variables

group=1
nullval=-999.0
naxis = 2                 ! no of dimensions
incs(1) = 1
incs(2) = 1               ! sampling interval of fits file

end Subroutine openfits
!--------------------------------------------------------
!-----------------------------------------------------------------------

Subroutine read_param(fileparam)
use main
implicit none

character*150 arg,opt
character*250 fileparam
character*1500 line
integer :: i

!read parameter file

open(88,file=trim(fileparam),status='old')

do i = 1,1000
   read(88,'(a)',end=24) line
   if(line(1:1).ne.'#')then
      read(line,*) opt,arg
      select case (opt)
      case ('-x')
         read(arg,*) ix
      case ('-y')
         read(arg,*) iy
      case ('-mag')
         read(arg,*) imag
      case('-fr')
         read(arg,*) ifr
      case('-fwhm')
         read(arg,*) ifwhm
      case('-flux')
         read(arg,*) iflux
      case('-ntot')
         read(arg,*) intot
      case('-nchipx')
         read(arg,*) nchipx
      case('-nchipy')
         read(arg,*) nchipy
      case('-order')
         read(arg,*) order
      case('-rgmax')
         read(arg,*) rgmax
      case('-back')
         read(arg,*) back
      end select
   end if
enddo
24 continue

if(rgmax>100)then
   write(*,*) 'WARNING: rgmax is limited to 100'
   write(*,*) 'Change irgmax in code to increase above 100'
   stop
end if

if(intot>1000)then
   write(*,*) 'WARNING:  Number of columns in SExtractor exceeds 1000 limit'
   stop
end if


end Subroutine read_param
