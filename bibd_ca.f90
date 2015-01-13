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
     if (nOpt(optSteps, IS) .eqv. .true.) then
        bibdFound=.true.
        return
     endif

     changeFactor=0

     i=randomInt((IS%VERTICES))+1
     j=randomInt((IS%BLOCKS))+1   

     if (dormant(IS,i,j)) then
        do vertex=1,IS%VERTICES
           if (IS%BLOCK_INTERSECTION(i,vertex)<IS%LAMBDA) then
              ! Can increase at most for v
              call increment(changefactor, 1)
           endif
        enddo
        if (IS%SUM_TOTAL<IS%VERTICES*IS%VERTICES_PER_BLOCK) then
           ! Can increase at most for (v*r-1)
           call increment(changefactor, (IS%VERTICES*IS%VERTICES_PER_BLOCK-IS%SUM_TOTAL))
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
             + (IS%VERTICES*IS%VERTICES_PER_BLOCK-1) &
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
        if (IS%SUM_TOTAL>IS%VERTICES*IS%VERTICES_PER_BLOCK) then
           ! Can increase at most for (v*b-(v*r-1))
           call increment(changefactor, (IS%SUM_TOTAL-IS%VERTICES*IS%VERTICES_PER_BLOCK))
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
             + (IS%VERTICES*IS%BLOCKS-(IS%VERTICES*IS%VERTICES_PER_BLOCK-1)) &
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

recursive logical function nOpt(n,IS) result (successfulOpt)
  use incidence_structure
  implicit none

  type(IncidenceStructure) IS
  integer n,i,j

  successfulOpt=.false.
  if(n<1) then
     successfulOpt=.false.
     return
  else if (abs(IS%SUM_TOTAL-IS%VERTICES*IS%VERTICES_PER_BLOCK)/=n) then
     if(n==1) then
        successfulOpt=.false.
        return
     else
        successfulOpt = nOpt(n-1,IS)
        return
     endif
  endif

  !call writeMatrix(IS)

  do i=1,IS%VERTICES
     do j=1,IS%BLOCKS
        if(& ! If it makes sence to do an n-opt step
             (dormant(IS,i,j) &
             .and. IS%SUM_TOTAL<=IS%VERTICES*IS%VERTICES_PER_BLOCK-n &
             .and. (IS%SUM_IN_ROW(i)<=IS%VERTICES_PER_BLOCK-n .or. IS%SUM_IN_COL(j)<=IS%INCIDENCES_PER_VERTICE-n) &
             )&
             .or.&
             (active(IS,i,j) &
             .and. IS%SUM_TOTAL>=IS%VERTICES*IS%VERTICES_PER_BLOCK+n &
             .and. (IS%SUM_IN_ROW(i)>=IS%VERTICES_PER_BLOCK+n .or. IS%SUM_IN_COL(j)>=IS%INCIDENCES_PER_VERTICE+n))&
             ) then
           call flip(IS,i,j)
           if(n == 1) then
              if (isBIBD(IS).eqv..True.) then
                 !print *, n, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else if (abs(IS%SUM_TOTAL-IS%VERTICES*IS%VERTICES_PER_BLOCK)==n - 1) then
              !write (*,*) "DEBUG: Calling n-1 opt", n
              if(nOpt(n-1,IS)) then
                 print *, n, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else
              call flip(IS,i,j)
              successfulOpt=.false.
              return
           endif
        endif
     enddo
  enddo
  !write (*,*) "DEBUG: end of nOpt", n
  return
end function nOpt
