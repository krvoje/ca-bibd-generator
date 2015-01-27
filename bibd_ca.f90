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
  if (command_argument_count() /= 4) then
     write (*,*) "Usage:"
     call get_command_argument(0,f)
     write (*,*) f ,"v k λ optSteps"
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
  endif

  call seedRandomGenerator()
  call construct(is,v,k,lmbd)

  call randomCA_BIBD(is,optSteps)
  print *, "An incidence matrix for the given parameters found:"
  call writeMatrix(is)  

  ! Memory cleanup
  call deconstruct(is)
end program bibd_ca

subroutine randomCA_BIBD(is, optSteps)
  use mtmod
  use incidence_structure
  
  implicit none

  type(IncidenceStructure) is
  integer optSteps

  integer changeFactor, maxChangeFactor
  integer row, col
  integer vertex, blok
  logical nOpt

  maxChangeFactor = is%v -1 + 3
  ! Rince and repeat until BIBD
  do while(.true.)

     if (isBIBD(is)) return
     if(nOpt(optSteps, optSteps, is)) return

     changeFactor = 0

     row = randomInt(is%v)+1
     col = randomInt(is%b)+1
     
     if (dormant(is,row,col)) then
        do vertex=1,is%v
           if(vertex==row) cycle
           if (is%ROW_INTERSECTION(row,vertex) < is%LAMBDA) call increment(changefactor, 1)
        enddo
        if (is%SUM_TOTAL < is%SUM_IDEAL) call increment(changefactor, 1)
        if (is%SUM_IN_COL(col) < is%k) call increment(changefactor, 1)
        if (is%SUM_IN_ROW(row) < is%r) call increment(changefactor, 1)

        if (is%SUM_TOTAL > is%SUM_IDEAL) call decrement(changefactor, 1)
        if (is%SUM_IN_COL(col) > is%k) call decrement(changefactor, 1)
        if (is%SUM_IN_ROW(row) > is%r) call decrement(changefactor, 1)
     else if (active(is,row,col)) then
        do vertex=1,is%v
           if(vertex==row) cycle
           if (is%ROW_INTERSECTION(row,vertex) > is%LAMBDA) call increment(changefactor, 1)
        enddo
        if (is%SUM_TOTAL > is%SUM_IDEAL) call increment(changefactor, 1)
        if (is%SUM_IN_COL(col) > is%k) call increment(changefactor, 1)
        if (is%SUM_IN_ROW(row) > is%r) call increment(changefactor, 1)

        if (is%SUM_TOTAL < is%SUM_IDEAL) call decrement(changefactor, 1)
        if (is%SUM_IN_COL(col) < is%k) call decrement(changefactor, 1)
        if (is%SUM_IN_ROW(row) < is%r) call decrement(changefactor, 1)
     endif
     if(randomInt(maxChangeFactor) < changeFactor) call flip(is,row,col)
  enddo
end subroutine randomCA_BIBD

recursive logical function nOpt(n,topOpt,is) result (successfulOpt)
  use incidence_structure
  implicit none

  type(IncidenceStructure) is
  integer n,row,col,topOpt  

  logical sumTotal_lt, sumTotal_gt
  logical sumInRow_lt, sumInRow_gt

  successfulOpt=.false.
  if(n<1) return
  if (abs(is%heuristic_distance) > is%max_hd_coef * n) return

  if (abs(is%SUM_IDEAL - is%SUM_TOTAL)/=n) then
     successfulOpt = nOpt(n-1,n-1,is)
     return
  endif
  
  sumTotal_lt = (is%SUM_TOTAL <= is%SUM_IDEAL - n)
  sumTotal_gt = (is%SUM_TOTAL >= is%SUM_IDEAL + n)
    
  do row=1,is%v
     sumInRow_lt = (is%SUM_IN_ROW(row) <= is%r - n)
     sumInRow_gt = (is%SUM_IN_ROW(row) >= is%r + n)
     do col=1,is%b
        if(&
             (dormant(is,row,col) &
             .and. sumTotal_lt &
             .and. (sumInRow_lt .or. is%SUM_IN_COL(col) <= is%k - n) &
             )&
             .or.&
             (active(is,row,col) &
             .and. sumTotal_gt &
             .and. (sumInRow_gt .or. is%SUM_IN_COL(col) >= is%k + n))&
             ) then
           call flip(is,row,col)
           sumInRow_lt = (is%SUM_IN_ROW(row) <= is%r-n)
           sumInRow_gt = (is%SUM_IN_ROW(row) >= is%r+n)
           if(n == 1) then
              if (isBIBD(is)) then
                 print *, topOpt, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else if(nOpt(n-1,topOpt,is)) then
              successfulOpt=.true.
              return
           else
              call flip(is,row,col)
              sumInRow_lt = (is%SUM_IN_ROW(row) <= is%r-n)
              sumInRow_gt = (is%SUM_IN_ROW(row) >= is%r+n)
              successfulOpt=.false.
              return
           endif
        endif
     enddo
  enddo
  return
end function nOpt
