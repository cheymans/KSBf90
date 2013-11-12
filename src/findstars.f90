! Program to select stars from the mag/fwhm and mag/fr locii
! Uses PGPLOT for the user to select the appropriate region
! Output file contains the selection criteria
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

integer:: imax,Ngals
parameter (imax=2e6)
real*4, dimension(1:imax) :: mag,objsize,fr,fwhm
real*8, dimension(1:imax) :: ra,dec
! parameter columns
integer :: ix,iy,imag,ifr,ifwhm,iflux,intot
integer :: rgmax,order,nchipx,nchipy
real*4 :: back

end Module Main
!-------------------------------------------------------------

Program findstars
use main
implicit none

integer:: i,j,k,ii

character*150 arg,opt
integer iargc,narg

! in/out files
character*250 filein, fileout,fileparam
character*500 line

! data
real*4, dimension(1:1000)  :: catalogue
real*4, dimension(1:imax) :: x,y
integer :: nstars
real*4 :: frmax

! plotting 
character*1 ch
character*20 ylabel
real*4 :: xl,xh,yl,yh
real*4 :: frl,frh,magl,magh
real*4 :: xx,yy

narg=iargc()

do i=1,narg
   call getarg(i,opt)
   call getarg(i+1,arg)
   select case (opt)
   case ('-h')
      print *,'findstars.a'
      print *,'Example: '
      print *,'findstars.a -in SEx.cat -param KSBf90.param -out star_criteria.dat'
      print *,'Update any changes to the SEXtractor catalogue in KSBf90.param'
      STOP

   case ('-in')
      filein=arg
   case ('-out')
      fileout=arg
   case ('-param')
      fileparam = arg

   end select
end do

! read parameter file
call read_param(fileparam)

open(44,file=filein,status='old')
open(45,file=fileout,status='unknown')

j = 0
frmax = 0.0
do i = 1,imax
   read(44,'(a)',end=23) line
   if(line(1:1).ne.'#')then
      j = j + 1
      read(line(1:500),*) (catalogue(k),k=1,intot)
      x(j) = catalogue(ix)
      y(j) = catalogue(iy)
      mag(j) =  catalogue(imag)
      fr(j) = catalogue(ifr)
      fwhm(j) = catalogue(ifwhm)
      ra(j) = catalogue(9)
      dec(j) = catalogue(10)

      ! find maximum galaxy size for PSF modelling with a limit
      ! to remove spurious size results that SExtractor can produce

      if(fr(j)>frmax.and.fr(j)<15.0)then
         frmax = fr(j)
      end if

   end if
end do
write(*,*) 'WARNING: Input file longer than imax in findstars.f90'
write(*,*) 'WARNING: Only inspecting the first imax=1e6 SExtractor objects'

23 continue
Ngals = j

write(*,*) 'Magnitude: min,max',ngals,minval(mag),maxval(mag)
write(*,*) 'FR: min,max',ngals,minval(fr(1:ngals)),maxval(fr(1:ngals))


call pgbegin(0,'/xwin',1,1)
call pgask(.false.)
call PGSFS(2)

21 continue

do ii = 1,2
 
! In the first pass this code displays the magnitude flux radius relation
! In the second pass the magnitude fhwm relation is drawn

   if(ii==1)then
      objsize = fr
      ylabel = 'Flux radius'
   else
      objsize = fwhm
      ylabel = 'FWHM'
   end if

! The following code structure is thanks to Ed Olding and Gavin Dalton

   ch = 'r'

   do while (ch.ne.'q'.and.ch.ne.'x'.and.ch.ne.'n')
      if (ch.eq.'r') then  ! redraw everything
         xl = minval(mag)
         xh = maxval(mag)
         yl = minval(objsize)
         yh = maxval(objsize)
         call pgsci(1)
         call newplot(xl,xh,yl,yh,mag,objsize,Ngals,ylabel)
         if(ii==2)then
            call plot_selected_stars(magl,magh,frl,frh)
         end if
      else if (ch.eq.'d') then
         go to 21
      else if (ch.eq.'l') then
         xl = xx
         yl = yy
         write(6,*) 'Hit l again to define zoom region'
         call pgsci(5)
         call pgband(2,1,xl,yl,xx,yy,ch)
         if (ch.eq.'l') then
            xh = max(xl,xx)
            yh = max(yl,yy)
            xl = min(xl,xx)
            yl = min(yl,yy)
            call pgsci(1)
            call newplot(xl,xh,yl,yh,mag,objsize,Ngals,ylabel)
            if(ii==2)then
               call plot_selected_stars(magl,magh,frl,frh)
            end if
         end if
         call pgsci(1)
      else if (ch.eq.'s') then
         xl = xx
         yl = yy
         write(6,*) 'Hit s again to define selection region'
         call pgsci(5)
         call pgband(2,1,xl,yl,xx,yy,ch)
         xh = max(xl,xx)
         yh = max(yl,yy)
         xl = min(xl,xx)
         yl = min(yl,yy)
         call pgsci(2)
         call pgrect(xl,xh,yl,yh)

         if(ii==1)then
            frl = yl
            frh = yh
            magl = xl
            magh = xh
         end if
         call plot_selected_stars(magl,magh,frl,frh)         
         call pgsci(1)

      end if
      
      call menu(xx,yy,ch)
   end do
   
   if(ch.eq.'q')then
      stop
   elseif (ch.eq.'x') then

! check inputs
      if(frl<0.5.or.frh>10.0.or. yl<0.5.or.yh>10.0)then
         write(*,*) 'size restrictions look too low or too high'
         write(*,*) 'flux radius', frl,frh
         write(*,*) 'fwhm', yl,yh
         write(*,*) 'Selection rejected'
         stop
      end if
      if(magl<0.5.or.magh>30.0.or. xl<0.5.or.xh>30.0)then
         write(*,*) 'magnitude restrictions look too low or too high'
         write(*,*) 'magnitude from flux radius selection', magl,magh
         write(*,*) 'magnitude from fwhm selection', xl,xh
         write(*,*) 'Selection rejected'
         stop
      end if
      write(*,*)' Findstars Size/Mag selection'
      write(*,*) 'flux radius', frl,frh
      write(*,*) 'fwhm', yl,yh
      write(*,*) 'magnitude', max(xl,magl), min(xh,magh)

      ! Count stars
      nstars = 0
      do i = 1,Ngals
         if(fr(i)<frh.and.fr(i)>frl.and.&
            fwhm(i)<yh.and.fwhm(i)>yl.and.&
            mag(i)>max(xl,magl).and.mag(i)<min(xh,magh)) then
            nstars = nstars + 1
            write(55,'(2f15.5,f7.2)') ra(i),dec(i),mag(i)
         end if
      end do
      write(*,*) 'Total star count', nstars
      if(nstars<10)then
         write(*,*) 'WARNING:  You have selected less than 10 stars'
         write(*,*) 'Selection rejected'
         stop
      end if

      write(45,107) nstars,frl, frh, yl, yh, max(xl,magl), min(xh,magh),frmax
   end if
end do
call pgend

107 format(I8, 7f12.5)
end Program findstars

!---------------------------------------------

subroutine newplot(xl,xh,yl,yh,x,y,Npt,ylabel)
implicit none
integer :: Npt
real*4:: xl,xh,yl,yh
real*4, dimension(1:Npt)::x,y
character*20 ylabel

call pgadvance
call pgenv(xl,xh,yl,yh,0,0)
call pglabel('Mag',ylabel,'Findstars')
call pgpt(Npt,x,y,1)

end subroutine newplot

!---------------------------------------------
subroutine menu(xx,yy,ch)
implicit none

  character(len=1) ch
  real*4 xx,yy
  write(6,*)
  write(6,*) '----------------------------------------------'
  write(6,*) '                 Menu Options'
  write(6,*) '----------------------------------------------'
  write(6,*) 'd -- delete everything and start again'
  write(6,*) 'l -- zoom in to user defined box'
  write(6,*) 'q -- exit without saving changes'
  write(6,*) 'r -- redraw field and delete current size/mag selection'
  write(6,*) 's -- select region'
  write(6,*) 'n -- next selection'
  write(6,*) 'x -- write current stellar selection and exit'
  write(6,*) 'z -- zoom in about cursor position'
  write(6,*) '----------------------------------------------'
  call pgcurse(xx,yy,ch)

end subroutine menu
!-----------------------------------------------------------------------

Subroutine plot_selected_stars(xl,xh,yl,yh)  
use main
implicit none
real*4 :: xl,xh,yl,yh
integer :: i

call pgsci(4)

do i = 1,Ngals
   if(fr(i)<yh.and.fr(i)>yl.and.mag(i)>xl.and.mag(i)<xh)then
      call pgpt(1,mag(i),objsize(i),1)
   end if
end do

end Subroutine plot_selected_stars
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

if(intot>1000)then
   write(*,*) 'WARNING:  Number of columns in SExtractor exceeds 1000 limit'
   stop
end if

end Subroutine read_param
