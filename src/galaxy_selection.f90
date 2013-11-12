! 
! Program to test average shear measurements
! as a function of galaxy magnitude, size, S/N etc
!
! Uses PGPLOT to make selection criteria
! Output file contains the galaxy selection criteria
!
! Catherine Heymans IPM School 2008
!
! Please acknowledge Catherine Heymans and Ludovic Van Waerbeke
! in any publications that make use of this code.
! 
!-------------------------------------------------------------
Module Main

integer:: imax
parameter (imax=1e7)

real*4, dimension(1:imax) :: mag,rg,fr,fwhm,SN
real*4, dimension(1:imax) :: x,y,flux,e1,e2,g1,g2,pgamma_fit
real*4, dimension(1:imax) :: wt,pgamma
real*8, dimension(1:imax) :: ra,dec
! parameter columns
integer :: ix,iy,imag,ifr,ifwhm,iflux,intot
integer :: rgmax,order,nchipx,nchipy
real*4 :: back

integer  :: nbins,ngals
end Module Main
!-------------------------------------------------------------

Program gal_selection
use main
implicit none

integer:: i,j,k,ii,jj,kk,icat

character*150 arg,opt
integer iargc,narg

! in/out files
character*500 filein, fileout, fileshear,fileparam
character*500 line

integer :: pgfitflag

! data
real*4, dimension(1:100)  :: catalogue
real*4 :: g1temp,g2temp,pgtemp
real*4 :: rgmin,magmax, SNmin,magmin

real*4, dimension(1:50) :: x_bin,y_bin,sd_bin
real*4  :: g1mean,g1sd,g2mean,g2sd
!
! set up defaults
nbins = 10

narg=iargc()

do i=1,narg
   call getarg(i,opt)
   call getarg(i+1,arg)
   select case (opt)
   case ('-h')
      print *,'galaxy_selection.a'
      print *,'Example: '
      print *,'galaxy_selection.a -x 7 -y 8 -flux 2 -mag 4 -fr 11 -fwhm 12 -ntot 14 -in shearcat.list -fit 1 -out selection.cat'
      print *,''
      print *,'       -in shearcat.list contains all the files that you wish to be evaluated together eg 1 STEP set'
      print *,'       Defaults: '
      print *,'       -nbins 10'
      print *,'       -fit 1 means fit Pgamma'
      print *,'       -fit 0 means use raw pgammas'
      STOP

   case ('-in')
      filein=arg
   case ('-out')
      fileout=arg
   case ('-param')
      fileparam=arg
   case('-nbins')
      read(arg,*) nbins
   case('-fit')
      read(arg,*) pgfitflag

   end select

enddo

call read_param(fileparam)

open(67,file=fileout,status='unknown')

i = 0
!open(44,file=filein,status='old')
do icat = 1,1!00
!   read(44,'(a)',end=22) fileshear

!   open(55, file=fileshear, status='old')
   open(55, file=filein, status='old')
   
   do while (i.le.imax)
      read(55,'(a)',end=23) line
      if(line(1:1).ne.'#')then
         read(line(1:500),*) (catalogue(k),k=1,intot+6)
         g1temp = 2.0*catalogue(intot+4)/catalogue(intot+6)
         g2temp = 2.0*catalogue(intot+5)/catalogue(intot+6)
         pgtemp = catalogue(intot+6)
         
         if(abs(g1temp)<1.0.and.abs(g2temp)<1.0.and.pgtemp>0.0)then
            i = i + 1
            ra(i)= catalogue(9)
            dec(i)= catalogue(10)
            x(i) = catalogue(ix)
            y(i) = catalogue(iy)
            mag(i) =  catalogue(imag)
            fr(i) = catalogue(ifr)
            fwhm(i) = catalogue(ifwhm)
            flux(i) = catalogue(iflux)
            rg(i) = catalogue(intot+1)
            pgamma(i) = catalogue(intot+6)
            SN(i) = min(catalogue(iflux)/catalogue(iflux + 1),100.0)
            e1(i) = catalogue(intot+4)
            e2(i) = catalogue(intot+5)
            wt(i) = 1.0
         end if
      end if
   end do
   write(*,*) icat, 'reached imax limit in galaxy_selection'
   23 continue
end do
22 continue

ngals = i
write(*,*) 'Total galaxies in systematics check', ngals

call pgopen('/xwin')
call pgsch(1.8)
call pgslw(2)
call pgask(.true.)
If(pgfitflag==1)then
   !Fit pgamma
   call fit_Pgamma(ngals,Pgamma,rg,wt)
   !pause
   ! calculate shear using fitted Pgamma
   do i = 1,ngals
      g1(i) = 2.0 * e1(i) / pgamma_fit(i)
      g2(i) = 2.0 * e2(i) / pgamma_fit(i)
   end do
else
   ! If not Pgamma fit use raw pgammas
   do i = 1,ngals
      g1(i) = 2.0 * e1(i) / pgamma(i)
      g2(i) = 2.0 * e2(i) / pgamma(i)
   end do
end If

call pgask(.false.)
call pgsubp(3,2)

do ii = 1,4

   call test_trends

   if(ii==1)then
      write(*,*) 'low size cut rg > ?'
      read(*,*) rgmin
      do i = 1,ngals
         if(rg(i)<rgmin)then
            wt(i) = 0.0
         end if
      end do
   elseif(ii==2)then
      write(*,*) 'faint mag cut mag < ?'
      read(*,*) magmax
      do i = 1,ngals
         if(mag(i)>magmax)then
            wt(i) = 0.0
         end if
      end do
   elseif(ii==3)then
      write(*,*) 'low SN cut SN > ?'
      read(*,*) SNmin
      do i = 1,ngals
         if(SN(i)<SNmin)then
            wt(i) = 0.0
         end if
      end do
   elseif(ii==4)then
      write(*,*) 'bright mag cut mag < ?'
      read(*,*) magmin
      do i = 1,ngals
         if(mag(i)<magmin)then
            wt(i) = 0.0
         end if
      end do
   end if
end do

call test_trends

call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g1,g2,ngals)
call binplot('g2','<g1>',x_bin,y_bin,sd_bin,nbins)
call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g2,g1,ngals)
call binplot('g1','<g2>',x_bin,y_bin,sd_bin,nbins)

call stats(ngals,g1,wt,g1mean,g1sd)
call stats(ngals,g2,wt,g2mean,g2sd)

write(*,*) 'Final galaxy selection gives'
write(*,107) g1mean,g1sd
107 format('g1 = ', f8.4, '+/-', f8.4)
write(*,108) g2mean,g2sd
108 format('g2 = ', f8.4, '+/-', f8.4)

! write out
close(55)
open(55, file=filein, status='old')
i = 0   
do while (i.le.imax)
   read(55,'(a)',end=24) line
   if(line(1:1).ne.'#')then
      read(line(1:500),*) (catalogue(k),k=1,intot+6)
      g1temp = 2.0*catalogue(intot+4)/catalogue(intot+6)
      g2temp = 2.0*catalogue(intot+5)/catalogue(intot+6)
      pgtemp = catalogue(intot+6)
         
      if(abs(g1temp)<1.0.and.abs(g2temp)<1.0.and.pgtemp>0.0)then
         i = i + 1


         if(wt(i)>0.0)then
            !  write(67,'(4f15.5,f7.3)') x(i),y(i),g1(i),g2(i),1.0
            !write(67,'(6e15.5)') ra(i),dec(i),mag(i),g1(i),g2(i)
            write(67,'(a,1x,2e15.5)')trim(line),g1(i),g2(i)
         end if
      end if
   end if
end do
24 continue
call pgend



end Program gal_selection
!--------------------------------------------------------------------

Subroutine bin_and_average(nbins,binned,sd_bin,xbin,wt,ydata,xdata,N)
implicit none

integer  :: N,nbins
real*4, dimension(1:nbins) :: xbin,binned,sd_bin
real*4, dimension(1:N) :: wt, ydata, xdata
real*4, dimension(1:nbins) :: sum_o_weights,ave, sumosq,ct
real*4   :: datamin,datamax
real*4   :: max,min,binsize,ymin,ymax

integer  :: i,j,k

datamax = maxval(xdata)
datamin = minval(xdata)
write(*,*)'in sub', datamax,datamin
binsize = (datamax - datamin)/nbins
ave = 0.0
sum_o_weights = 0.0
sumosq = 0.0
ct = 0.0

do i = 1,nbins

   min = datamin + (i-1)*binsize
   max = datamin + i*binsize

   xbin(i) = min + (max-min)/2.0

   do j = 1,N
      if(xdata(j)<max.and.xdata(j).ge.min.and.wt(j)>0.0)then
         sum_o_weights(i) = sum_o_weights(i) + wt(j)
         ave(i) = ave(i) + ydata(j)*wt(j)
         sumosq(i) = sumosq(i) + ydata(j)**2
         ct(i) = ct(i) + 1
      end if
   end do
end do

ymin = 1000000.0
ymax = -1000000.0

! normalise

do i = 1,nbins
   if(sum_o_weights(i)>0.0.and.ct(i)>25.0)then
      binned(i) = ave(i)/sum_o_weights(i)
      if(abs(ct(i)-sum_o_weights(i))<0.001)then
         ! no weights
         sd_bin(i) = sqrt(sumosq(i)/(ct(i)*1.0) - (ave(i)/ct(i))**2)/sqrt(ct(i)*1.0)
      else
         sd_bin(i) = sqrt(1.0 / sum_o_weights(i))
      end if
   else
      binned(i) = -99.99
      sd_bin(i) = 0.99
   end if
   write(*,*) i,ct(i),binned(i),sd_bin(i)
end do

end Subroutine bin_and_average
!-------------------------------------------------------------

Subroutine binplot(xlabel,ylabel,xpt,ypt,sd,N)
implicit none
integer ::i,N
character(len=4) :: xlabel,ylabel
real*4, dimension(1:N) ::xpt,ypt,sd
real*4 :: ymin,ymax
! calculate ymin and ymax
ymin = 1000
ymax = -1000
do i = 1,N
   if(ypt(i)>-99.99.and.sd(i)<0.02.and.ypt(i)-sd(i)<ymin)then
      ymin = ypt(i) - sd(i)
   end if
   if(sd(i)<0.02.and.ypt(i)+sd(i)>ymax)then
      ymax = ypt(i) + sd(i)
   end if
end do

call pgenv(0.9*minval(xpt),1.1*maxval(xpt),ymin,ymax, 0,0)
call pglab(xlabel,ylabel,'')
call pgpt(N,xpt,ypt,4)
call pgerry(N,xpt,ypt+sd,ypt-sd,1.0)

end Subroutine binplot
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
end Subroutine read_param
!-----------------------------------------------------------------------

Subroutine test_trends
use main
implicit none

! bindata

real*4, dimension(1:50) :: x_bin,y_bin,sd_bin

call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g1,mag,ngals)
call binplot('mag ','<g1>',x_bin,y_bin,sd_bin,nbins)

call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g1,rg,ngals)
call binplot('size','<g1>',x_bin,y_bin,sd_bin,nbins)

call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g1,SN,ngals)
call binplot('StoN','<g1>',x_bin,y_bin,sd_bin,nbins)

call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g2,mag,ngals)
call binplot('mag ','<g2>',x_bin,y_bin,sd_bin,nbins)

call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g2,rg,ngals)
call binplot('size','<g2>',x_bin,y_bin,sd_bin,nbins)

call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,g2,SN,ngals)
call binplot('StoN','<g2>',x_bin,y_bin,sd_bin,nbins)


!call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,pgamma,mag,ngals)
!call binplot('mag','pgamma',x_bin,y_bin,sd_bin,nbins)

!call bin_and_average(nbins,y_bin,sd_bin,x_bin,wt,pgamma,rg,ngals)
!call binplot('size','pgamma',x_bin,y_bin,sd_bin,nbins)

end Subroutine test_trends
!------------------------------------------
Subroutine stats(N,ydata,wt,mean,sd)
implicit none
integer  :: i,N,ct
real*4, dimension(1:N) :: wt, ydata
real*4 :: mean,sd,sumosq,ave

sumosq = 0.0
ave = 0.0
ct = 0

do i = 1,N
   if(wt(i)>0.0)then
      ave = ave + ydata(i)*wt(i)
      sumosq = sumosq + ydata(i)**2
      ct = ct + 1
   end if
end do
if(ct>5)then
   mean = ave/(ct*1.0)
   sd = sqrt(sumosq/(ct*1.0) - (ave/(ct*1.0))**2)/sqrt(ct*1.0)
  ! write(*,*) mean,sd, ct
else
   write(*,*) 'less than 5 objects in the statistics test'
   stop
end if

end Subroutine stats

!-----------------------------------------

Subroutine fit_Pgamma(N,Pg,rsize,weight)
use main
implicit none
integer :: i,j,k,N
real*4, dimension(1:N) :: Pg, rsize,weight,binweight
integer :: nobin
parameter(nobin = 10)
real*4, dimension(1:nobin) :: pg_bin,rg_bin,sd_bin
real*4 :: datamax, datamin,binsize
!numrec fitting
integer :: nparam
parameter(nparam = 4)
integer, dimension(1:nparam) :: ia
real*4, dimension(1:nparam,1:nparam) :: covar
real*4, dimension(1:nparam) :: fitparam
real*4 :: chisq
real*4, dimension(1:1000) :: xpt,ypt

datamax = maxval(rsize)
datamin = minval(rsize)
write(*,*) datamax, datamin
binsize = (datamax - datamin)/nobin
call bin_and_average(nobin,pg_bin,sd_bin,rg_bin,weight,Pg,rsize,N)

! 3 sigma clip - some very high pgamma values
binweight = weight
do i = 1,N
   j = nint((rsize(i) - datamin)/binsize)
   if(abs(pgamma(i)-pg_bin(j))>2.5*sd_bin(j))then
      binweight(i)=0.0
   end if
end do

call bin_and_average(nobin,pg_bin,sd_bin,rg_bin,binweight,pgamma,rsize,N)

fitparam = 0
fitparam(1) = 1.0
ia =(/1,1,0,0/)  ! 2nd order 

do i = 1,nobin
   write(*,*) rg_bin(i),pg_bin(i), sd_bin(i)
end do
write(*,*) nparam, ia, nobin
call LFIT(rg_bin,pg_bin,sd_bin,nobin,fitparam,IA,nparam,COVAR,nparam,CHISQ)

call pgenv(datamin,datamax,0.0,2.0,1,0)
call pglab('rg','pgamma','Fitting Pgamma')
call pgpt(N,rsize,Pg,1)
call pgsci(6)
call pgpt(nobin,rg_bin,pg_bin,4)
call pgerry(nobin,rg_bin,pg_bin+sd_bin,pg_bin-sd_bin,1.0)

do i = 1,100
   xpt(i) = i*0.1
   ypt(i) = fitparam(1) + fitparam(2)*xpt(i) + &
        fitparam(3)*xpt(i)**2 + fitparam(4)*xpt(i)**3 
end do

call pgsci(4)
call pgline(100,xpt,ypt)

do i = 1,N
   pgamma_fit(i) = fitparam(1) + fitparam(2)*rsize(i) + &
                   fitparam(3)*rsize(i)**2 + fitparam(4)*rsize(i)**3 
end do

pause
call pgsci(1)
end Subroutine fit_Pgamma

!**************************************************************
! Numerical Recipes
!**************************************************************
!-----------------------------------------------------------------------

SUBROUTINE funcs(xval,afunc,ma)
  integer :: ma
  real*4  :: xval,afunc(1:ma)

  afunc(1)=1.
  afunc(2)=xval
  afunc(3)=xval**2
  afunc(4)=xval**3

END SUBROUTINE funcs

!-----------------------------------------------------------------------
     SUBROUTINE LFIT(Xa,Y,SIG,NDAT,A,IA,MA,COVAR,npc,CHISQ)
      integer ma,ia(1:ma),npc,ndat,mmax
      PARAMETER (MMAX=50)
      real*4 chisq,Xa(1:NDAT),Y(1:NDAT),SIG(1:NDAT),A(1:MA)
      real*4 COVAR(1:npc,1:npc)

      integer  i,j,k,l,m,mfit
      real*4 sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)


      mfit = 0
      do j = 1,ma
         if(ia(j).ne.0) mfit = mfit + 1
      end do

      if(mfit.eq.0)then
         write(*,*) ma,'lfit: no parameters to be fitted - test', ia
         stop
      end if
      DO 14 J=1,MFIT
        DO 13 K=1,MFIT
          COVAR(J,K)=0.
13      CONTINUE
        BETA(J)=0.
14    CONTINUE


      DO I=1,ndat
        CALL FUNCS(Xa(I),AFUNC,MA)
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
        CALL FUNCS(Xa(I),AFUNC,MA)
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

