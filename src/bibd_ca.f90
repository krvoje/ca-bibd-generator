! A cellular automaton that tries to generate a BIBD
! for the given parameters
program bibd_ca

  use incidence_structure
  use utils
  use tabu_list

  implicit none

  type(IncidenceStructure) is
  type(TabuList) rows_tabu, cols_tabu
  integer optSteps
  integer v,k,lmbd
  integer i,j
  real r_r,b_r
  integer tabu_size
  character(100) f

  ! If the number of command line params is invalid, print instructions
  if (command_argument_count() /= 5) then
     write (*,*) "Usage:"
     call get_command_argument(0,f)
     write (*,*) f ,"v k λ optSteps tabuSize"
     stop
  else
     call get_command_argument(1,f)
     read (f,"(I10)") v
     call get_command_argument(2,f)
     read (f,"(I10)") k
     call get_command_argument(3,f)
     read (f,"(I10)") lmbd
     call get_command_argument(4,f)
     read (f,"(I10)") optSteps
     call get_command_argument(5,f)
     read (f,"(I10)") tabu_size
  endif

  call seedRandomGenerator()
  call construct(is,v,k,lmbd)

  if(optSteps > is%v * is%b) optSteps = is%v*is%b

  call mkTabuList(rows_tabu, tabu_size)
  call mkTabuList(cols_tabu, tabu_size)

  call randomCA_BIBD(is, optSteps, rows_tabu, cols_tabu)
  print *, "An incidence matrix for the given parameters found:"
  call writeMatrix(is)

  ! Memory cleanup
  call deconstruct(is)
  call delTabuList(rows_tabu)
  call delTabuList(cols_tabu)
end program bibd_ca

subroutine randomCA_BIBD(is, optSteps, rows_tabu, cols_tabu)
  use mtmod
  use incidence_structure
  use tabu_list
  use n_opt

  implicit none

  type(IncidenceStructure) is
  type(TabuList) rows_tabu, cols_tabu
  integer optSteps

  integer changeFactor, maxChangeFactorActive, maxChangeFactorDormant
  integer maxChangeFactor
  integer row, col, i, point
  integer inc_ratio, min_dim, max_dim, max_point
  integer lambda_row_delta
  integer lambda_col_delta
  
  min_dim = min(is%v, is%b)
  max_dim = max(is%v, is%b)
  max_point = max(is%r, is%k)
  inc_ratio = is%v / is%k

  maxChangeFactorDormant = inc_ratio*(is%v + is%b - 2) + is%sum_ideal + is%r + is%k
  maxChangeFactorActive = inc_ratio*(is%v + is%b - 2) + is%sum_total + is%v + is%b


  ! Rince and repeat until BIBD
  do while(.true.)
     call increment(is%generations,1)
     
     if (isBIBD(is)) return

     do i = optSteps, optSteps, -1
        if(nOpt(i, i, is)) return
     enddo

     row = randomInt(is%v)+1
     col = randomInt(is%b)+1

     changeFactor = 0

     do while(is_in_list(rows_tabu, row))
        row = randomInt(is%v)+1    
     enddo
     call push(rows_tabu, row)

     do while(is_in_list(cols_tabu, col))
        col = randomInt(is%b)+1
     enddo
     call push(cols_tabu, col)

     ! Combined loop for 1:is%v and 1:is%b
     do point=1,max_dim
        lambda_row_delta = is%lambda - is%row_intersection(row, point)
        if(point <= is%v) then
            if(lambda_row_delta /= 0) then
                call increment(changefactor, inc_ratio)
            else
                if(is%v >= is%b) call decrement(changeFactor, inc_ratio)
            endif
        endif

        lambda_col_delta = is%lambda - is%col_intersection(col, point)
        if(point <= is%b) then
            if(lambda_col_delta /= 0) then
                call increment(changefactor, inc_ratio)
            else
                if(is%v < is%b) call decrement(changeFactor, inc_ratio)
            endif
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
