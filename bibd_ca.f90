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
  logical bibdFound, nOpt

  integer nothingChanged

  nothingChanged=0

  bibdFound=.false.
  ! Rince and repeat until BIBD
  do while(bibdFound .eqv. .false.)
     if (isBIBD(is)) then
        bibdFound=.true.
        return
     endif     
     if(nOpt(optSteps, optSteps, is)) then
        bibdFound=.true.
        return        
     endif

     changeFactor=0

     row=randomInt(is%v)+1
     col=randomInt(is%b)+1
        
     if (dormant(is,row,col)) then
        do vertex=1,is%v
           if(vertex==row) cycle
           if (is%ROW_INTERSECTION(row,vertex) < is%LAMBDA) then
              call increment(changefactor, is%LAMBDA - is%ROW_INTERSECTION(row,vertex))
           endif
        enddo
        if (is%SUM_TOTAL < is%SUM_IDEAL) then
           call increment(changefactor, is%SUM_IDEAL - is%SUM_TOTAL)
        endif
        if (is%SUM_IN_COL(col)<is%k) then
           call increment(changefactor, is%k)
        endif        
        if (is%SUM_IN_ROW(row)<is%r) then
           call increment(changefactor, is%r)
        endif
        maxChangeFactor = is%LAMBDA*(is%v-1) &
             + is%SUM_IDEAL &
             + is%r &
             + is%k
        if(randomInt(maxChangeFactor) < changeFactor) then
           call flip(is,row,col)
        endif
     else if (active(is,row,col)) then
        do vertex=1,is%v
           if(vertex==row) cycle
           if (is%ROW_INTERSECTION(row,vertex) > is%LAMBDA) then
              call increment(changefactor, is%ROW_INTERSECTION(row,vertex) - is%LAMBDA)
           endif
        enddo
        if (is%SUM_TOTAL > is%SUM_IDEAL) then
           call increment(changefactor, is%SUM_TOTAL - is%SUM_IDEAL)
        endif
        if (is%SUM_IN_COL(col)>is%k) then
           call increment(changefactor, is%v)
        endif
        if (is%SUM_IN_ROW(row)>is%r) then
           call increment(changefactor, is%b)
        endif
        maxChangeFactor = is%b*(is%v-1) &
             + is%SUM_TOTAL &
             + is%b &
             + is%v
        if(randomInt(maxChangeFactor) < changeFactor) then
           call flip(is,row,col)
        endif
     else
        print *, "Severe fatal error"
        call exit(1)
     endif
  enddo
end subroutine randomCA_BIBD

recursive logical function nOpt(n,topOpt,is) result (successfulOpt)
  use incidence_structure
  implicit none

  type(IncidenceStructure) is
  integer n,i,j,topOpt  

  logical sumTotal_lt, sumTotal_gt
  logical sumInRow_lt, sumInRow_gt

  successfulOpt=.false.
  if(n<1) then
     successfulOpt=.false.
     return
  endif

  if (abs(is%heuristic_distance)>is%max_hd_coef * n) then
     successfulOpt=.false.
     return
  endif

  ! print *, is%heuristic_distance, is%heuristic_distance/is%max_hd_coef
  ! call writeMatrix(is)
  
  if (abs(is%SUM_IDEAL - is%SUM_TOTAL)/=n) then
     successfulOpt = nOpt(n-1,n-1,is)
     return
  endif
  
  !call writeMatrix(is)

  sumTotal_lt = is%SUM_TOTAL<=is%SUM_IDEAL - n
  sumTotal_gt = is%SUM_TOTAL>=is%SUM_IDEAL + n
    
  do i=1,is%v
     sumInRow_lt=is%SUM_IN_ROW(i)<=is%r-n
     sumInRow_gt=is%SUM_IN_ROW(i)>=is%r+n
     do j=1,is%b
        if(&
             (dormant(is,i,j) &
             .and. sumTotal_lt &
             .and. (sumInRow_lt .or. is%SUM_IN_COL(j)<=is%k-n) &
             )&
             .or.&
             (active(is,i,j) &
             .and. sumTotal_gt &
             .and. (sumInRow_gt .or. is%SUM_IN_COL(j)>=is%k+n))&
             ) then
           call flip(is,i,j)
           sumInRow_lt=is%SUM_IN_ROW(i)<=is%r-n
           sumInRow_gt=is%SUM_IN_ROW(i)>=is%r+n
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
              call flip(is,i,j)
              sumInRow_lt=is%SUM_IN_ROW(i)<=is%r-n
              sumInRow_gt=is%SUM_IN_ROW(i)>=is%r+n
              successfulOpt=.false.
              return
           endif
        endif
     enddo
  enddo
  !write (*,*) "DEBUG: end of nOpt", n
  return
end function nOpt
