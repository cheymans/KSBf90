Program set_obs_group
implicit none
integer :: i,j,k,ntile,obs_group,days

! see Figure 1 of Heymans et al 2008 for the PSF groupings

open(45,file='dates.dat',status='old')
open(46,file='tile_group.dat',status='unknown')

do i = 1,80
   read(45,*) ntile,days

   if(days<4)then  
      obs_group = 1
   elseif(days<7)then  
      obs_group = 2
   elseif(days<10)then
      obs_group = 3
   elseif(days<12)then  
      obs_group = 4
   elseif(days<15)then
      obs_group = 5
   elseif(days<20)then
      obs_group =6
   else  ! 2 gyro mode
      obs_group = 7
   end if

   write(46,*) ntile, obs_group

end do

end Program set_obs_group
