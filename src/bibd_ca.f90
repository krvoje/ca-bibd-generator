! A cellular automaton that tries to generate a BIBD
! for the given parameters  
program bibd_ca

  use incidence_structure
  use utils
  
  implicit none
  
  type(IncidenceStructure) is
  integer optSteps
  integer v,k,lmbd
  integer i,j
  real r_r,b_r
  character(100) f

  ! If the number of command line params is invalid, print instructions
  if (command_argument_count() /= 3 .and. command_argument_count() /= 4) then
     write (*,*) "Usage:"
     call get_command_argument(0,f)
     write (*,*) f ,"v k λ [optSteps]"
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
        read (f,"(I10)") optSteps
     else
        optSteps=2
     endif
  endif

  call seedRandomGenerator()
  call construct(is,v,k,lmbd)

  if(optSteps > is%v * is%b) optSteps = is%v*is%b

  call randomCA_BIBD(is,optSteps)
  print *, "An incidence matrix for the given parameters found:"
  call writeMatrix(is)  

  ! Memory cleanup
  call deconstruct(is)
end program bibd_ca

subroutine randomCA_BIBD(is, optSteps)
  use mtmod
  use incidence_structure
  use n_opt
  
  implicit none

  type(IncidenceStructure) is
  integer optSteps

  integer changeFactor, maxChangeFactorActive, maxChangeFactorDormant
  integer maxChangeFactor
  integer row, col, i, point
  integer incRatio

  maxChangeFactorDormant = (is%b/is%v)*(is%v + is%b - 2) + is%sum_ideal + is%r + is%k
  maxChangeFactorActive = (is%b/is%v)*(is%v + is%b - 2) + is%sum_total + is%v + is%b

  incRatio=max(is%v, is%b) / min(is%v, is%b)
  
  ! Rince and repeat until BIBD
  do while(.true.)
     call increment(is%generations,1)
     if (isBIBD(is))&
        return     

     do i = optSteps, optSteps, -1
        if(nOpt(i, i, is)) return
     enddo

     changeFactor = 0

     row = randomInt(is%v)+1
     col = randomInt(is%b)+1

     ! Not really sure why, but this block below really speeds up the convergence     
     do point=1,max(is%v,is%b)
        if(point <= is%v .and. is%row_intersection(row,point) /= is%lambda) then
           call increment(changefactor, incRatio)
        else
            if(is%v >= is%b) call decrement(changeFactor, incRatio)
        endif

        if(point <= is%b .and. is%col_intersection(col,point) /= is%lambda) then
           call increment(changefactor, incRatio)
        else
           if(is%b > is%v) call decrement(changeFactor, incRatio)
        endif
     enddo
     
     if (dormant(is,row,col)) then
        if (is%sum_total /= is%sum_ideal) call increment(changefactor, is%sum_ideal - is%sum_total)
        if (is%sum_in_row(row) /= is%r) call increment(changefactor, is%r - is%sum_in_row(row))
        if (is%sum_in_col(col) /= is%k) call increment(changefactor, is%k - is%sum_in_col(col))
        
        if(randomInt(maxChangeFactorDormant) < changeFactor) call flip(is,row,col)
     else if (active(is,row,col)) then
        if (is%sum_total /= is%sum_ideal) call increment(changefactor, is%sum_total - is%sum_ideal)
        if (is%sum_in_row(row) /= is%r) call increment(changefactor, is%sum_in_row(row) - is%r)
        if (is%sum_in_col(col) /= is%k) call increment(changefactor, is%sum_in_col(col) - is%k)

        if(randomInt(maxChangeFactorActive) < changeFactor) call flip(is,row,col)
     endif
  enddo
end subroutine randomCA_BIBD
