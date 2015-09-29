! A cellular automaton that tries to generate a BIBD
! for the given parameters
program bibd_ca

  use incidence_structure
  use utils

  implicit none

  type(IncidenceStructure) is
  integer opt_steps
  integer v,k,lmbd
  integer i,j
  real r_r,b_r
  character(100) f

  ! If the number of command line params is invalid, print instructions
  if (command_argument_count() /= 3 .and. command_argument_count() /= 4) then
     write (*,*) "Usage:"
     call get_command_argument(0,f)
     write (*,*) f ,"v k λ [opt_steps]"
     stop
  else
     call get_command_argument(1,f)
     read (f,"(I10)") v
     call get_command_argument(2,f)
     read (f,"(I10)") k
     call get_command_argument(3,f)
     read (f,"(I10)") lmbd
     if(command_argument_count() == 4) then
         call get_command_argument(4,f)
         read (f,"(I10)") opt_steps
     else
         opt_steps=2
     endif
  endif

  call seedRandomGenerator()
  call construct(is,v,k,lmbd)

  if(opt_steps > is%v * is%b) opt_steps = is%v*is%b

  call randomCA_BIBD(is,opt_steps)

  print *, "An incidence matrix for the given parameters found:"
  call writeMatrix(is)

  ! Memory cleanup
  call deconstruct(is)
end program bibd_ca

subroutine randomCA_BIBD(is, opt_steps)
  use mtmod
  use incidence_structure
  use n_opt

  implicit none

  type(IncidenceStructure) is
  integer opt_steps

  integer maxChangeFactor, changeFactor
  integer row, active_col, dormant_col
  integer i, point
  integer inc_ratio, min_dim, max_dim
  integer cfa, cfd

  integer v,r,b,k,lambda

  v=is%v
  r=is%r
  b=is%b
  k=is%k
  lambda=is%lambda

  min_dim = min(is%v, is%b)
  max_dim = max(is%v, is%b)

  maxChangeFactor = (v-1)*(abs(r-lambda)) + b

  ! Rince and repeat until BIBD
  do while(.true.)
     call increment(is%generations,1)
     if (isBIBD(is)) return
     if(nOpt(opt_steps, opt_steps, is)) return

     ! Pick a random row, and an active and a dormant cell in it
     row = randomInt(is%v)+1
     active_col = randomInt(is%b)+1
     do while(.not.active(is,row,active_col))
        active_col = randomInt(is%b)+1
     enddo

     dormant_col = randomInt(is%b)+1
     do while(.not.dormant(is,row,dormant_col))
        dormant_col = randomInt(is%b)+1
     enddo

     cfD = changeFactor(is,row,dormant_col)
     cfA = changeFactor(is,row,active_col)

     i=randomInt(maxChangefactor)
     if(i < cfD&
          .and. i < cfA) then
        ! Exchange places
        call flip(is,row,active_col)
        call flip(is,row,dormant_col)
     endif
     
  enddo
end subroutine randomCA_BIBD

integer function changeFactor(is,row,col)  
  use incidence_structure
  implicit none
  
  type(IncidenceStructure) is
  integer row,col
  
  integer delta
  integer i

  changefactor = 0

  ! v * (abs(lambda - r))
  do i=1, is%v
     if(i .eq. row) cycle
     delta = is%lambda - is%row_intersection(row, i)

     if(active(is,row,col).and.active(is,i,col).and.delta<0) then
        call decrement(changefactor, delta)
     else if(active(is,row,col).and.dormant(is,i,col).and.delta<0) then
        call increment(changefactor, delta)
     else if(dormant(is,row,col).and.active(is,i,col).and.delta>0) then
        call increment(changefactor, delta)
     else if(dormant(is,row,col).and.dormant(is,i,col).and.delta>0) then
        call decrement(changefactor, delta)
     endif
  enddo

  !do i=1, is%b
  !   if(i .eq. col) cycle
  !   delta = is%lambda - is%col_intersection(col, i)
  !   if(active(is,row,col).and.active(is,row,i).and.delta<0) then
  !      call decrement(changefactor, delta)
  !   else if(dormant(is,row,col).and.active(is,col,i).and.delta>0) then
  !      call increment(changefactor, delta)
  !   endif
  !enddo

  delta = is%k - is%sum_in_col(col)
  if (active(is,row,col)) then
     call decrement(changefactor, delta)
  else if (dormant(is,row,col)) then
     call increment(changefactor, delta)
  endif

  return
end function changeFactor
