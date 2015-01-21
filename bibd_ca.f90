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

  integer changeFactor, maxChangeFactor
  integer row, col
  integer vertex, blok
  logical bibdFound, nOpt

  integer nothingChanged

  nothingChanged=0

  bibdFound=.false.
  ! Rince and repeat until BIBD
  do while(bibdFound .eqv. .false.)
     if (isBIBD(IS) .or. nOpt(optSteps, optSteps, IS)) then
        bibdFound=.true.
        return
     endif

     changeFactor=0

     row=randomInt(IS%VERTICES)+1
     col=randomInt(IS%BLOCKS)+1
        
     if (dormant(IS,row,col)) then
        do vertex=1,IS%VERTICES
           if(vertex==row) cycle
           if (IS%ROW_INTERSECTION(row,vertex) < IS%LAMBDA) then
              call increment(changefactor, IS%LAMBDA - IS%ROW_INTERSECTION(row,vertex))
           endif
        enddo
        if (IS%SUM_TOTAL < IS%SUM_IDEAL) then
           call increment(changefactor, IS%SUM_IDEAL - IS%SUM_TOTAL)
        endif
        if (IS%SUM_IN_COL(col)<IS%INCIDENCES_PER_VERTICE) then
           call increment(changefactor, IS%INCIDENCES_PER_VERTICE)
        endif        
        if (IS%SUM_IN_ROW(row)<IS%VERTICES_PER_BLOCK) then
           call increment(changefactor, IS%VERTICES_PER_BLOCK)
        endif
        maxChangeFactor = IS%LAMBDA*(IS%VERTICES-1) &
             + IS%SUM_IDEAL &
             + IS%INCIDENCES_PER_VERTICE &
             + IS%VERTICES_PER_BLOCK
        if(randomInt(maxChangeFactor) < changeFactor) then
           call flip(IS,row,col)
        endif
     else if (active(IS,row,col)) then
        do vertex=1,IS%VERTICES
           if(vertex==row) cycle
           if (IS%ROW_INTERSECTION(row,vertex) > IS%LAMBDA) then
              call increment(changefactor, IS%ROW_INTERSECTION(row,vertex) - IS%LAMBDA)
           endif
        enddo
        if (IS%SUM_TOTAL > IS%SUM_IDEAL) then
           call increment(changefactor, IS%SUM_TOTAL - IS%SUM_IDEAL)
        endif
        if (IS%SUM_IN_COL(col)>IS%INCIDENCES_PER_VERTICE) then
           call increment(changefactor, IS%BLOCKS)
        endif
        if (IS%SUM_IN_ROW(row)>IS%VERTICES_PER_BLOCK) then
           call increment(changefactor, IS%VERTICES)
        endif
        maxChangeFactor = IS%BLOCKS*(IS%VERTICES-1) &
             + IS%SUM_TOTAL &
             + IS%BLOCKS &
             + IS%VERTICES
        if(randomInt(maxChangeFactor) < changeFactor) then
           call flip(IS,row,col)
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
  
  if (abs(IS%SUM_IDEAL - IS%SUM_TOTAL)/=n) then
     successfulOpt = nOpt(n-1,n-1,IS)
     return
  endif

  !call writeMatrix(IS)

  sumTotal_lt = IS%SUM_TOTAL<=IS%SUM_IDEAL - n
  sumTotal_gt = IS%SUM_TOTAL>=IS%SUM_IDEAL + n
    
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
