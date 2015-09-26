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

  integer maxChangeFactorActive, maxChangeFactorDormant
  integer row, col
  integer i, point
  integer inc_ratio, min_dim, max_dim

  min_dim = min(is%v, is%b)
  max_dim = max(is%v, is%b)
  inc_ratio = 1

  maxChangeFactorDormant = inc_ratio*(is%v + is%b - 2) + is%sum_ideal + is%r + is%k
  maxChangeFactorActive = inc_ratio*(is%v + is%b - 2) + is%sum_total + is%v + is%b

  ! Rince and repeat until BIBD
  do while(.true.)
     call increment(is%generations,1)
     if (isBIBD(is)) return
     if(nOpt(opt_steps, opt_steps, is)) return

     row = randomInt(is%v)+1
     col = randomInt(is%b)+1

     if(dormant(is,row,col)&
         .and. randomInt(maxChangeFactorDormant) < changeFactor(is,row,col,inc_ratio)) then
        call flip(is,row,col)
     endif

     if(active(is,row,col)&
         .and. randomInt(maxChangeFactorActive) < changeFactor(is,row,col,inc_ratio)) then
        call flip(is,row,col)
     endif

  enddo
end subroutine randomCA_BIBD
