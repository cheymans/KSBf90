module Calculate_Moments
  implicit none
!****************************************
!Contains the nuts and bolts KSB shape measurement routines
!****************************************

contains

Subroutine getshape(xs,ys,fluxs,rwindow,gsflag,xedge,yedge,e,psm,psh, q, q4, Renormalise_Moments)
Use Main; use Moments, only:Weighted_Intensity_Moment_0, Weighted_Intensity_Moment_2, Weighted_Intensity_Moment_4
Implicit none

real*4         :: xs,ys,rwindow,fluxs
integer        :: istar,gsflag
real*4         :: xedge,yedge
real*4,intent(out)         :: q(1:2,1:2), q4(2,2,2,2)
logical,intent(in)::Renormalise_Moments

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
if(Renormalise_Moments) then
   q = q/q0
   q4 = q4/q0
end if

End Subroutine getshape




end module Calculate_Moments
