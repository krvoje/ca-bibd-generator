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

  integer changeFactor, maxChangeFactorActive, maxChangeFactorDormant
  integer row, col
  integer vertex, blok
  logical nOpt

  maxChangeFactorDormant = (is%v - 1) * is%LAMBDA + is%SUM_IDEAL + is%r + is%k
  maxChangeFactorActive = (is%v - 1) * (is%v - is%LAMBDA) + is%SUM_TOTAL + is%v + is%b
  
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
           if (is%ROW_INTERSECTION(row,vertex) < is%LAMBDA)&
                call increment(changefactor, is%LAMBDA - is%ROW_INTERSECTION(row,vertex))
           !if (is%ROW_INTERSECTION(row,vertex) > is%LAMBDA) call decrement(changefactor, 1)
        enddo
        if (is%SUM_TOTAL < is%SUM_IDEAL) call increment(changefactor, is%SUM_IDEAL - is%SUM_TOTAL)
        if (is%SUM_IN_ROW(row) < is%r) call increment(changefactor, is%r - is%SUM_IN_ROW(row))
        if (is%SUM_IN_COL(col) < is%k) call increment(changefactor, is%k - is%SUM_IN_COL(col))

        if (is%SUM_TOTAL > is%SUM_IDEAL) call decrement(changefactor, is%SUM_TOTAL - is%SUM_IDEAL)
        if (is%SUM_IN_ROW(row) > is%r) call decrement(changefactor, is%SUM_IN_ROW(row) - is%r)
        if (is%SUM_IN_COL(col) > is%k) call decrement(changefactor, is%SUM_IN_COL(col) - is%k)
        if(randomInt(maxChangeFactorDormant) < changeFactor) call flip(is,row,col)
     else if (active(is,row,col)) then
        do vertex=1,is%v
           if(vertex==row) cycle
           if (is%ROW_INTERSECTION(row,vertex) > is%LAMBDA)&
                call increment(changefactor, is%ROW_INTERSECTION(row,vertex) - is%LAMBDA)              
           !if (is%ROW_INTERSECTION(row,vertex) < is%LAMBDA) call decrement(changefactor, 1)
        enddo
        if (is%SUM_TOTAL > is%SUM_IDEAL) call increment(changefactor, is%SUM_TOTAL - is%SUM_IDEAL)
        if (is%SUM_IN_ROW(row) > is%r) call increment(changefactor, is%SUM_IN_ROW(row) - is%r)
        if (is%SUM_IN_COL(col) > is%k) call increment(changefactor, is%SUM_IN_COL(col) - is%k)

        if (is%SUM_TOTAL < is%SUM_IDEAL) call decrement(changefactor, is%SUM_IDEAL - is%SUM_TOTAL)
        if (is%SUM_IN_ROW(row) < is%r) call decrement(changefactor, is%r - is%SUM_IN_ROW(row))
        if (is%SUM_IN_COL(col) < is%k) call decrement(changefactor, is%k - is%SUM_IN_COL(col))
        if(randomInt(maxChangeFactorActive) < changeFactor) call flip(is,row,col)        
     endif
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
