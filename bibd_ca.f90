! A cellular automaton that tries to generate a BIBD
! for the given parameters  
program bibd_ca

  use incidence_structure
  use utils
  
  implicit none
  
  type(IncidenceStructure) IS
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
  call construct(IS,v,k,lmbd)

  call randomCA_BIBD(IS,optSteps)
  print *, "An incidence matrix for the given parameters found:"
  call writeMatrix(IS)  

  ! Memory cleanup
  call deconstruct(IS)
end program bibd_ca

subroutine randomCA_BIBD(IS, optSteps)
  use mtmod
  use incidence_structure
  
  implicit none

  type(IncidenceStructure) IS
  integer optSteps

  integer changeFactor, maxChangeFactor, i, j, vertex
  logical bibdFound, nOpt

  bibdFound=.false.
  ! Rince and repeat until BIBD
  do while(bibdFound .eqv. .false.)
     if (isBIBD(IS) .or. nOpt(optSteps, optSteps, IS)) then
        bibdFound=.true.
        return
     endif

     changeFactor=0

     i=randomInt(IS%VERTICES)+1
     j=randomInt(IS%BLOCKS)+1

     if (dormant(IS,i,j)) then
        do vertex=1,IS%VERTICES
           if (IS%BLOCK_INTERSECTION(i,vertex)<IS%LAMBDA) then
              ! Can increase at most for v
              call increment(changefactor, 1)
           endif
        enddo
        if (IS%SUM_TOTAL_LESS_THAN_IDEAL) then
           ! Can increase at most for (v*r-1)
           call increment(changefactor, abs(IS%IDEAL_SUM - IS%SUM_TOTAL))
        endif
        if (IS%SUM_IN_COL(j)<IS%INCIDENCES_PER_VERTICE) then
           ! Can increase at most for (k-1)
           call increment(changefactor, (IS%INCIDENCES_PER_VERTICE-IS%SUM_IN_COL(j)))
        endif
        if (IS%SUM_IN_ROW(i)<IS%VERTICES_PER_BLOCK) then
           ! Can increase at most for (r-1)
           call increment(changefactor, (IS%VERTICES_PER_BLOCK-IS%SUM_IN_ROW(i)))
        endif
        maxChangeFactor=IS%VERTICES &
             + (IS%IDEAL_SUM-1) &
             + (IS%INCIDENCES_PER_VERTICE-1) &
             + (IS%VERTICES_PER_BLOCK-1)
        if(randomInt(maxChangeFactor) < changeFactor) then
           call flip(IS,i,j)
        endif
     else if (active(IS,i,j)) then
        do vertex=1,IS%VERTICES
           if (IS%BLOCK_INTERSECTION(i,vertex)>IS%LAMBDA) then
              ! Can increase at most for v              
              call increment(changefactor, 1)
           endif
        enddo
        if (IS%SUM_TOTAL_MORE_THAN_IDEAL) then
           ! Can increase at most for (v*b-(v*r-1))
           call increment(changefactor, abs(IS%IDEAL_SUM - IS%SUM_TOTAL))
        endif
        if (IS%SUM_IN_COL(j)>IS%INCIDENCES_PER_VERTICE) then
           ! Can increase at most for (v-(k-1))
           call increment(changefactor, (IS%SUM_IN_COL(j)-IS%INCIDENCES_PER_VERTICE))
        endif
        if (IS%SUM_IN_ROW(i)>IS%VERTICES_PER_BLOCK) then
           ! Can increase at most for (b-(r-1))
           call increment(changefactor, (IS%SUM_IN_ROW(i)-IS%VERTICES_PER_BLOCK))
        endif
        maxChangeFactor=IS%VERTICES &
             + (IS%VERTICES_X_BLOCKS-(IS%IDEAL_SUM-1)) &
             + (IS%VERTICES-(IS%INCIDENCES_PER_VERTICE-1)) &
             + (IS%BLOCKS-(IS%VERTICES_PER_BLOCK-1))
        if(randomInt(maxChangeFactor) < changeFactor) then
           call flip(IS,i,j)
        endif
     else
        print *, "Severe fatal error"
        call exit(1)
     endif
  enddo
end subroutine randomCA_BIBD

recursive logical function nOpt(n,topOpt,IS) result (successfulOpt)
  use incidence_structure
  implicit none

  type(IncidenceStructure) IS
  integer n,i,j,topOpt  

  logical sumTotal_lt, sumTotal_gt
  logical sumInRow_lt, sumInRow_gt

  successfulOpt=.false.
  if(n<1) then
     successfulOpt=.false.
     return
  endif
  
  if (abs(IS%IDEAL_SUM - IS%SUM_TOTAL)/=n) then
     successfulOpt = nOpt(n-1,n-1,IS)
     return
  endif

  !call writeMatrix(IS)

  sumTotal_lt = IS%SUM_TOTAL<=IS%IDEAL_SUM - n
  sumTotal_gt = IS%SUM_TOTAL>=IS%IDEAL_SUM + n
    
  do i=1,IS%VERTICES
     sumInRow_lt=IS%SUM_IN_ROW(i)<=IS%VERTICES_PER_BLOCK-n
     sumInRow_gt=IS%SUM_IN_ROW(i)>=IS%VERTICES_PER_BLOCK+n
     do j=1,IS%BLOCKS
        if(&
             (dormant(IS,i,j) &
             .and. sumTotal_lt &
             .and. (sumInRow_lt .or. IS%SUM_IN_COL(j)<=IS%INCIDENCES_PER_VERTICE-n) &
             )&
             .or.&
             (active(IS,i,j) &
             .and. sumTotal_gt &
             .and. (sumInRow_gt .or. IS%SUM_IN_COL(j)>=IS%INCIDENCES_PER_VERTICE+n))&
             ) then
           call flip(IS,i,j)
           sumInRow_lt=IS%SUM_IN_ROW(i)<=IS%VERTICES_PER_BLOCK-n
           sumInRow_gt=IS%SUM_IN_ROW(i)>=IS%VERTICES_PER_BLOCK+n
           if(n == 1) then
              if (isBIBD(IS)) then
                 print *, topOpt, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else if(nOpt(n-1,topOpt,IS)) then
              successfulOpt=.true.
              return
           else
              call flip(IS,i,j)
              sumInRow_lt=IS%SUM_IN_ROW(i)<=IS%VERTICES_PER_BLOCK-n
              sumInRow_gt=IS%SUM_IN_ROW(i)>=IS%VERTICES_PER_BLOCK+n
              successfulOpt=.false.
              return
           endif
        endif
     enddo
  enddo
  !write (*,*) "DEBUG: end of nOpt", n
  return
end function nOpt
