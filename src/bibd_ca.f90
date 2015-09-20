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

  integer changeFactor, maxChangeFactorActive, maxChangeFactorDormant
  integer maxChangeFactor
  integer row, col, i, point
  integer inc_ratio, min_dim, max_dim
  integer lambda_row_delta
  integer lambda_col_delta

  min_dim = min(is%v, is%b)
  max_dim = max(is%v, is%b)
  inc_ratio = max_dim / min_dim

  maxChangeFactorDormant = inc_ratio*(is%v + is%b - 2) + is%sum_ideal + is%r + is%k
  maxChangeFactorActive = inc_ratio*(is%v + is%b - 2) + is%sum_total + is%v + is%b


  ! Rince and repeat until BIBD
  do while(.true.)
     call increment(is%generations,1)
     if (isBIBD(is)) return
     if(nOpt(opt_steps, opt_steps, is)) return

     changeFactor = 0

     row = randomInt(is%v)+1
     col = randomInt(is%b)+1

     do point=1, is%v
        lambda_row_delta = is%lambda - is%row_intersection(row, point)
        if(lambda_row_delta /= 0) then
            call increment(changefactor, inc_ratio)
        else
            if(is%v >= is%b) call decrement(changeFactor, inc_ratio)
        endif
     enddo
     do point=1, is%b
        lambda_col_delta = is%lambda - is%col_intersection(col, point)
        if(lambda_col_delta /= 0) then
            call increment(changefactor, inc_ratio)
        else
            if(is%v < is%b) call decrement(changeFactor, inc_ratio)
        endif
     enddo

     if (dormant(is,row,col)) then
        call increment(changefactor, is%sum_ideal - is%sum_total)
        call increment(changefactor, is%r - is%sum_in_row(row))
        call increment(changefactor, is%k - is%sum_in_col(col))

        if(randomInt(maxChangeFactorDormant) < changeFactor) call flip(is,row,col)
     else if (active(is,row,col)) then
        call decrement(changefactor, is%sum_ideal - is%sum_total)
        call decrement(changefactor, is%r - is%sum_in_row(row))
        call decrement(changefactor, is%k - is%sum_in_col(col))

        if(randomInt(maxChangeFactorActive) < changeFactor) call flip(is,row,col)
     endif
  enddo
end subroutine randomCA_BIBD
