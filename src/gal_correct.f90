! Program to measure galaxy shapes
! and correct galaxies with the polynomial psf model 
! created by psffit
!
! Can use PGPLOT to show diagnostic plots
! Output file contains the input SExtractor catalogue with additional
! shape and shear parameters appended
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
parameter (imax=1e7)

real*4, dimension(1:imax) :: mag,rg,fr,fwhm,SN
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
real*4 :: xccd,yccd

!numrec fitting
integer                       :: nparams, irgmax
parameter (nparams = 10, irgmax = 100)
integer, dimension(1:nparams) :: ia
real*4, dimension(1:nparams,1:irgmax,1:10,1:10) :: pfit1,pfit2
real*4, dimension(1:irgmax,1:10,1:10) :: shsmfact

real*4, dimension(1:imax) :: x,y,flux
integer, dimension(1:imax) :: irg

! parameter columns
integer :: ix,iy,imag,ifr,ifwhm,iflux,intot
integer :: rgmax,order,nchipx,nchipy
real*4 :: back

! fits image
real, dimension(:,:), allocatable :: image

end Module Main
!-------------------------------------------------------------

Program gal_correct
use main; use Calculate_Moments, only:getshape
implicit none

integer:: i,j,k,ii,jj,kk
integer :: iix,iiy,icx,icy
integer :: header
integer :: ngals

character*150 arg,opt
integer iargc,narg

! in/out files
character*250 filein, fileout,filepsf, filefits,filecrit, fileparam
character*1000 line,string

! data
real*4, dimension(1:1000)  :: catalogue
real*4 ::xc,yc,xedge,yedge
integer, dimension(1:imax) :: flag
real*4:: frl,frh,fwhml,fwhmh,magl,magh,frmax
real*4 :: rgbin

! getshapes output
real*4                  :: e(0:1),ecor(0:1),shear(0:1), Pgamma
real*4                  :: psm(0:1,0:1)
real*4                  :: psh(0:1,0:1)
real*4, dimension(0:1) :: aveshear, sumosq, sd, mean


! psf correction
Real*4, dimension(1:2) :: pfitted
Real*4, dimension(0:1) :: eopsm

!postage stamp array
integer, dimension(1:2) :: fpixels,lpixels
integer :: mpix,npix

character*50 :: toplabel


!--RRG ADDITIONS--!
real*4::KSB_Q2(2,2), KSB_Q4(2,2,2,2)

narg=iargc()

do i=1,narg
   call getarg(i,opt)
   call getarg(i+1,arg)
   select case (opt)
   case ('-h')
      print *,'gal_correct.a'
      print *,'Example: '
      print *,'galcorrect.a -image image.fits -in SEx.cat -psf PSF.dat -crit starcriteria.dat -param KSBf90.param -out shape.cat'
      print *,''
      STOP

   case ('-in')
      filein=arg
   case ('-image')
      filefits=arg
   case ('-crit')
      filecrit=arg
   case ('-param')
      fileparam=arg
   case ('-psf')
      filepsf=arg
   case ('-out')
      fileout=arg
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

open(44,file=filein,status='old')
open(33,file=filecrit,status='old')
open(55,file=filepsf,status='old')
open(45,file=fileout,status='unknown')

aveshear = 0.0
sumosq = 0.0
kk = 0

! open the image files
call openfits(filefits)

!calculate CCD chip widths
xccd = naxes(1)*1.0/nchipx
yccd = naxes(2)*1.0/nchipy

!uncomment for any diagnostic plots
!! open pgplot window
!call pgopen('/xwin')
!call pgsch(1.8)
!call pgslw(2)
!call pgsubp(6,4)

! read psf file containing the chip by chip polynomial model
do i = 1,rgmax
   do icx = 1,nchipx
      do icy = 1,nchipy
         read(55,*) ii,iix,iiy,shsmfact(ii,iix,iiy),&
                    (pfit1(j,ii,iix,iiy),j=1,nparams),&
                    (pfit2(j,ii,iix,iiy),j=1,nparams)
         if(i.ne.ii.or.iix.ne.icx.or.iiy.ne.icy)then
            write(*,*) 'error reading psf file'
            write(*,*) 'check nchipx and nchipy are correct'
            write(*,*) i,icx,icy
            write(*,*) ii,iix,iiy
            stop
         end if
      end do
   end do
end do

! read criteria file
read(33,*) i,frl,frh,fwhml,fwhmh,magl,magh,frmax

! set binning scale for the varying weight function
rgbin = frh/rgmax

header = 0
do i = 1,imax
   read(44,'(a)',end=23) line
   if(line(1:1).ne.'#')then

      if(header == 0)then  ! add to header
         write(45,'(a)') '#  +1 rg          KSBf90 output'
         write(45,'(a)') '#  +2 e1_obs      KSBf90 output'
         write(45,'(a)') '#  +3 e2_obs      KSBf90 output'
         write(45,'(a)') '#  +4 e1_cor      KSBf90 output'
         write(45,'(a)') '#  +5 e2_cor      KSBf90 output'
         write(45,'(a)') '#  +6 pgamma      KSBf90 output'
         header = 1
      end if

      read(line(1:500),*) (catalogue(k),k=1,intot)

      x(i) = catalogue(ix)
      y(i) = catalogue(iy)
      mag(i) =  catalogue(imag)
      fr(i) = catalogue(ifr)
      fwhm(i) = catalogue(ifwhm)
      flux(i) = catalogue(iflux)
      SN(i) = catalogue(iflux)/catalogue(iflux + 1)
      rg(i) = max(min(2.0*frh,fr(i)),frh)
      irg(i) = max(nint((rg(i) - frh)/rgbin),1)

      ! set CCD chip
      ixchip(i) = int(x(i)/xccd) + 1
      iychip(i) = int(y(i)/yccd) + 1
            
      if(fr(i)>frh.and.mag(i)>magl)then  ! galaxy

         ! select postage stamp from fits image
         ! deal with edges
         
         fpixels(1) = max(1,nint(x(i) - m_size/2)) 
         fpixels(2) = max(1,nint(y(i) - n_size/2)) 
         lpixels(1) = min(nint(x(i) + m_size/2 - 1),naxes(1))
         lpixels(2) = min(nint(y(i) + n_size/2 - 1),naxes(2))
         
         mpix = lpixels(1) - fpixels(1) + 1 
         npix = lpixels(2) - fpixels(2) + 1
          
         ! put all objects in a regular sized grid and subtract off 
         ! flat background noise set by back param in KSBf90.param
         
         object = 0.0
         do ii = 0,mpix-1
            do jj = 0,npix-1
               object(ii,jj) = image(fpixels(1)+ii,fpixels(2)+jj) - back
             end do
         end do
         
         xc = x(i) - fpixels(1) 
         yc = y(i) - fpixels(2) 
         xedge = min(xc, naxes(1)-x(i))
         yedge = min(yc, naxes(2)-y(i))

         !   uncomment to see 2D images of the selected object
         !   call plotarray(m_size,n_size,object)
         !   call pgsci(6)
         !   call pgpt(1,xc,yc,4)
         
         ! getshape is the nuts and bolts of KSB
                  
         KSB_Q2 = 0.0; KSB_Q4 = 0.0
         call getshape(xc,yc,flux(i),rg(i),flag(i),xedge,yedge,e,psm,psh, KSB_Q2, KSB_Q4)

         call calcp(i,irg(i),pfitted)
         
         eopsm(0) = psm(0,0)*pfitted(1)  + psm(0,1)*pfitted(2)
         eopsm(1) = psm(1,0)*pfitted(1)  + psm(1,1)*pfitted(2) 
         
         do j = 0,1
            ecor(j) = e(j) - eopsm(j)
         end do
         
         !Hoekstra correction can be included
         !psh = psh*(1.0- (ecor(0)**2 + ecor(1)**2)/2.0)

         Pgamma = (psh(0,0) - psm(0,0)*shsmfact(irg(i),ixchip(i),iychip(i))) &
                + (psh(1,1) - psm(1,1)*shsmfact(irg(i),ixchip(i),iychip(i))) 
         
         do j = 0,1
            shear(j) = 2.0*ecor(j)/Pgamma
         end do

         !Here are my selection criteria of what I consider to be
         !a good galaxy for the shear to be calculated for the 
         !GREAT08 challenge.  For astronomy in general I would
         !use a SN>10 cut
         
         j = 1
         if(flag(i)==0.and.rg(i)>1.15*frh.and.&
            abs(e(0))<0.5.and.abs(e(1))<0.5.and.&
            abs(shear(0))<1.0.and.abs(shear(1))<1.0.and.SN(i)>5.0)then
            kk = kk + 1
            aveshear = aveshear + shear
            sumosq = sumosq + shear**2
!            write(45,'(4e15.5,I4)') x(i),y(i),shear(0),shear(1),j

         end if
      
         if(flag(i)==0.and.abs(e(0))<0.5.and.abs(e(1))<0.5)then
            write(string,102) rg(i), e(0), e(1), ecor(0), ecor(1), Pgamma
102         format(f8.2,5f10.4)
            write(45,'(a)') trim(line)// trim(string)
         end if
      end if
   else
!      write(45,'(a)') trim(line)
   end if
end do

write(*,*) 'WARNING: catalogue longer than imax in gal_correct.f90'

23 continue

ngals = i - 1

! Close pgplot window
! uncomment for any diagnostic plots
!call pgend

mean = aveshear/kk
sd = sqrt(sumosq/(kk*1.0) - mean**2)/sqrt(kk*1.0)

!write(*,109) (mean(i),sd(i),i=0,1), kk
!109 format(4f15.7, I8)

!GREAT08 Submission format
write(*,110) (mean(i),i=0,1),(sd(i),i=0,1), kk
110 format(4f15.7,I7)

end Program 

!******************************************************************
!RRG nuts and bolts
!******************************************************************
! Authored by CAJD on 13th November 2013.
!******************************************************************

integer function factorial(n)
  integer::n
  
  integer::count
  
  factorial = 1
  do count = n, 1, -1
     factorial = factorial*count
  end do
  return
end function factorial

subroutine RRG_Isotropic_Correction(PSF_rwindow, rwindow, Q2, RRG_Q2)
  !--Implements equation (49) of RRG, corrects for the isotropic (guassian) part of the PSF. Q2 is the galaxy quadropole--!
  !--Note: Does this need implemented before/after the anisotropic - probably after--!
  real*4,intent(in)::PSF_rwindow, rwindow
  real*4, intent(in)::Q2(2,2)
  real*4, intent(out):: RRG_Q2(2,2)

  integer::i,j, Delta(2,2)
  
  Delta = 0.0; Delta(1,1) = 1.0; Delta(2,2) = 1.0
  
  do i = 1, 2
     do j = 1, 2
        RRG_Q2(i,j) = ((PSF_rwindow/rwindow)**4.0)*(Q2(i,j) - rwindow*rwindow*Delta(i,j))
     end do
  end do

end subroutine RRG_Isotropic_Correction

subroutine RRG_Anisotropic_Correction(rwindow, Q2, Q4, PSF_Q2, PSF_Q4, RRG_Q2)
  !--Implements equation (46) of RRG to correct the galaxy quadropole moment--!
  real*4::rwindow
  real*4::Q2(2,2), PSF_Q2(2,2), PSF_Q4(2,2,2,2), Q4(2,2,2,2)
  real*4::RRG_Q2(2,2)

  integer::i,j,k,l, ll
  real*4::C(2,2,2,2), Delta(2,2), CQ4(2,2,2,2),  Permutation_Set(4)
  real*4::  Permutations(24,4)

!  allocate(Permutations(factorial(size(Permutation_Set)), size(Permutation_Set))); Permutations = 0

  !--Calculate Corrected 4th order moment (CQ4), given in equation [50]--!
  Delta = 0.0
  Delta(1,1) = 1.0; Delta(2,2) = 1.0
  CQ4 = 0.0
  do i = 1, 2
     do j = 1, 2
        do k = 1, 2
           do l= 1, 2
              Permutation_Set = (/i,j,k,l/)
              call permutate(Permutation_Set, Permutations)
              CQ4(i,j,k,l) = Q4(i,j,k,l) - PSF_Q4(i,j,k,l)
              do ll = 1, size(Permutations,1)
                 CQ4(i,j,k,l) = CQ4(i,j,k,l) - 6.e0_4*PSF_Q2(Permutations(ll,1), Permutations(ll,2))*Q2(Permutations(ll,3),Permutations(ll,4)) + 6.e0_4*PSF_Q2(Permutations(ll,1), Permutations(ll,2))*PSF_Q2(Permutations(ll,3),Permutations(ll,4))
              end do
           end do
        end do
     end do
  end do

  !--Calculate Matrix C [Equation (40)]--!
  C = 0.0
  do i = 1, 2
     do j = 1, 2
        do k = 1, 2
           do l = 1, 2
              C(i,j,k,l) = Delta(i,k)*Delta(j,k) -(1.e0_4/(rwindow*rwindow))*(Q2(k,i)*Delta(j,l) + Q2(k,j)*Delta(i,l)) + (1.e0_4/(2.e0_4*rwindow*rwindow))*(CQ4(i,j,k,l) - Q2(i,j)*Q2(k,l))
           end do
        end do
     end do
  end do

  !--Equation (46)
  do i = 1, 2
     do j = 1, 2
        RRG_Q2(i,j) = Q2(i,j)
        !--Summation over k,l assumed--!
        do k = 1, 2
           do l = 1, 2
              RRG_Q2(i,j) = RRG_Q2(i,j) - C(i,j,k,l)*PSF_Q2(k,l)
           end do
        end do
     end do
  end do

end subroutine RRG_Anisotropic_Correction

  recursive subroutine permutate(E, P)
    integer, intent(in)  :: E(:)       ! array of objects                       
    integer, intent(out) :: P(:,:)     ! permutations of E                      
    integer  :: N, Nfac, i, k, S(size(P,1)/size(E), size(E)-1) 
    N = size(E); Nfac = size(P,1);
    do i=1,N                           ! cases with E(i) in front               
      if( N>1 ) call permutate((/E(:i-1), E(i+1:)/), S)
      forall(k=1:Nfac/N) P((i-1)*Nfac/N+k,:) = (/E(i), S(k,:)/)
    end do
  end subroutine permutate

!!$!****************************************************************
!!$! KSB nuts and bolts
!!$!****************************************************************
!!$!
!!$! KSBf90 version of getshape
!!$! With thanks to Nick Kaiser for the original imcat C version 
!!$! of the getshape subroutine.
!!$!
!!$!****************************************************************
!!$
!!$Subroutine getshape(xs,ys,fluxs,rwindow,gsflag,xedge,yedge,e,psm,psh)
!!$Use Main
!!$Implicit none
!!$
!!$real*4         :: xs,ys,rwindow,fluxs
!!$integer        :: istar,gsflag
!!$real*4         :: xedge,yedge
!!$
!!$real*4                        :: q(0:1,0:1), denom
!!$integer                       :: i0, j0, i, j, di, dj, rmax, l, m
!!$real*4                        :: r, dx, dy
!!$real*4                        :: W, Wp, Wpp, fc, DD, DD1, DD2
!!$real*4                        :: Xsm(0:1,0:1), Xsh(0:1,0:1)
!!$real*4                        :: em(0:1), eh(0:1)
!!$integer                       :: R_MAX_FACTOR
!!$integer                       :: negflux
!!$integer                       :: N1, N2
!!$
!!$real*4                  :: xim(0:1)
!!$real*4                  :: d(0:1)
!!$real*4                  :: e(0:1)
!!$real*4                  :: psm(0:1,0:1)
!!$real*4                  :: psh(0:1,0:1)
!!$real*4                  :: fb0
!!$real*4, dimension(0:1)     :: dfb
!!$
!!$if(rwindow.le.0.0)then
!!$   write(*,*) 'problem rwindow = 0 in getshapes'
!!$end if
!!$
!!$xim(0) = xs
!!$xim(1) = ys
!!$
!!$R_MAX_FACTOR = 4
!!$
!!$N1 = n_size
!!$N2 = m_size  
!!$
!!$rmax = rwindow * R_MAX_FACTOR + 1       ! round up integer
!!$i0 = xim(1)                          ! round down integer
!!$j0 = xim(0)                          !  ""
!!$
!!$! zero arrays
!!$
!!$eh = 0.0
!!$em = 0.0
!!$e = 0.0
!!$Xsh = 0.0
!!$Xsm = 0.0
!!$psm = 0.0
!!$psh = 0.0
!!$q = 0.0
!!$d = 0.0
!!$   
!!$if(rmax>xedge.or.rmax>yedge)then
!!$   gsflag = 1
!!$   go to 88
!!$else
!!$   do i = max(i0-rmax,0),min(i0+rmax,N2)
!!$      do j = max(j0-rmax,0),min(j0+rmax,N1)
!!$         
!!$         di = i - i0
!!$         dj = j - j0
!!$         dx = j - xim(0)
!!$         dy = i - xim(1)
!!$         
!!$         r = sqrt(dx**2 + dy**2)
!!$         
!!$         if (r <= rmax)then
!!$            W = exp(-0.5 * r**2 / rwindow**2)
!!$            Wp = -0.5 * W / rwindow**2
!!$            Wpp = 0.25 * W / rwindow**4
!!$            
!!$            fc = object(j,i)
!!$            
!!$!            if(gsflag.ge.2.and.object(j,i)<0.0)then
!!$!               fc = (object(j-1,i-1) + object(j-1,i) + &
!!$!                    object(j-1,i+1) + object(j,i-1) + &
!!$!                    object(j,i+1)+ object(j+1,i-1)+ &
!!$!                    object(j+1,i)+ object(j+1,i+1)) / 8.0
!!$!            end if
!!$            
!!$            d(0) = d(0) + (W * fc * dx)
!!$            d(1) = d(1) + (W * fc * dy)
!!$            q(0,0) = q(0,0) + (fc * W * dx * dx)
!!$            q(1,1) = q(1,1) + (fc * W * dy * dy)
!!$            q(0,1) = q(0,1) + (fc * W * dx * dy)
!!$            q(1,0) = q(1,0) + (fc * W * dx * dy)
!!$            DD = di * di + dj * dj
!!$            DD1 = dj * dj - di * di
!!$            DD2 = 2 * di * dj
!!$            Xsm(0,0) = Xsm(0,0) + &
!!$                       (2 * W + 4 * Wp * DD + 2 * Wpp * DD1 * DD1) * fc 
!!$            Xsm(1,1) = Xsm(1,1) + &
!!$                       (2 * W + 4 * Wp * DD + 2 * Wpp * DD2 * DD2) * fc 
!!$            Xsm(0,1) = Xsm(0,1) + (2 * Wpp * DD1 * DD2 * fc) 
!!$            Xsm(1,0) = Xsm(1,0) + (2 * Wpp * DD1 * DD2 * fc)
!!$            em(0) = em(0) + (4 * Wp + 2 * Wpp * DD) * DD1 * fc
!!$            em(1) = em(1) + (4 * Wp + 2 * Wpp * DD) * DD2 * fc
!!$            Xsh(0,0) = Xsh(0,0) + (2 * W * DD + 2 * Wp * DD1 * DD1) * fc
!!$            Xsh(1,1) = Xsh(1,1) + (2 * W * DD + 2 * Wp * DD2 * DD2) * fc
!!$            Xsh(0,1) = Xsh(0,1) + (2 * Wp * DD1 * DD2 * fc)
!!$            Xsh(1,0) = Xsh(1,0) + (2 * Wp * DD1 * DD2 * fc)
!!$            eh(0) = eh(0) + (2 * Wp * DD * DD1 * fc)
!!$            eh(1) = eh(1) + (2 * Wp * DD * DD2 * fc)
!!$
!!$         end if
!!$      end do
!!$   end do
!!$
!!$   !normalise d
!!$   
!!$   if (fluxs > 0)then 
!!$      do l = 0,1
!!$         d(l) = d(l)/fluxs
!!$      end do
!!$      negflux = 0
!!$   else 
!!$      negflux = 1
!!$   end if
!!$   
!!$   
!!$   ! calculate ellipticities
!!$   denom = q(0,0) + q(1,1)                 
!!$   
!!$   if (denom > 0)then
!!$      e(0) = (q(0,0) - q(1,1)) / denom
!!$      e(1) = (q(0,1) + q(1,0)) / denom
!!$      em = em / denom
!!$      eh = eh / denom
!!$      eh = eh + (2 * e)
!!$
!!$      do l = 0,1
!!$         do m = 0,1
!!$            psm(l,m) = 0.5 * (Xsm(l,m) / denom - e(l) * em(m))
!!$            psh(l,m) = Xsh(l,m) / denom - e(l) * eh(m)
!!$         end do
!!$      end do
!!$   else
!!$      gsflag = 4
!!$   end if
!!$end if
!!$
!!$88 continue
!!$
!!$!write(*,*) 'e',(e(i),i=0,1),gsflag
!!$!write(*,*) (d(i),i=0,1)
!!$!write(*,*) 'psm',(psm(i,0),i=0,1),(psm(i,1),i=0,1)
!!$!write(*,*) 'psh',(psh(i,0),i=0,1),(psh(i,1),i=0,1)
!!$
!!$End subroutine getshape

!-------------------------------------------------------------------------

!********************************************************************
! Fitting and Statistics routines
!********************************************************************

!---------------------------------------------------------------

Subroutine calcp(i,k,pfitted)
Use Main
Implicit none

integer,intent(in)         :: i
integer :: j,k
real*4          :: afunc(1:nparams)
Real*4, dimension(1:2) ::pfitted

call funcs(x(i),y(i),afunc,nparams)

pfitted = 0.0

do j = 1,nparams
   pfitted(1) = pfitted(1) + ia(j)*pfit1(j,k,ixchip(i),iychip(i))*afunc(j)
   pfitted(2) = pfitted(2) + ia(j)*pfit2(j,k,ixchip(i),iychip(i))*afunc(j)
end do

end Subroutine calcp
!-----------------------------------------------------------------------


!********************************************************************
! PGPLOT plotting routines
!********************************************************************


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

!**************************************************************
! Numerical Recipes
! You should have a on site numerical recipes licence
! to use these subroutines
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


!***************************************************************
! FITSIO ROUTINES 
!***************************************************************

Subroutine openfits(filename)
use main
implicit none

character(len=200) filename
integer             :: readwrite,blocksize
integer             :: nfound
integer, dimension(1:2) :: fpixels,lpixels

! The status and file unit parameters initialized.

status = 0      

! Open the read only FITS file 
! Denoted by region, fieldname and camerano from do loop above
! Note blocksize is anobsolete ignorable parameter

readwrite=0

unit = 99
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

allocate (image(1:naxes(1),1:naxes(2)))

! FITSIO routine

fpixels = 1
lpixels = naxes

call ftgsve(unit,group,naxis,naxes,fpixels,lpixels,incs, &
     nullval,image,anyf,status)

! Check for any error, and if so print out error messages.

if (status .gt. 0)then
   write(*,*) filename,&
        'readobj error problem, fitsio status = ', status
end if

! Closing FITS file 
      
call ftflus(unit,status)
call ftclos(unit,status)
end Subroutine openfits
!--------------------------------------------------------
!-----------------------------------------------------------------------

Subroutine read_param(fileparam)
use main
implicit none

character*150 arg,opt
character*250 fileparam
character*500 line
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
