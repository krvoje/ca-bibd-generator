! A cellular automaton that tries to generate a BIBD
! for the given parameters  
program bibd_ca

  use incidence_structure
  use utils
  
  implicit none
  
  type(IncidenceStructure) incs
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

  call construct(incs,v,k,lmbd)

  call seedRandomGenerator()
  
  ! Rince and repeat until BIBD
  do while(isBIBD(incs) .eqv. .false.)     
     call randomCA_BIBD(incs,optSteps)
     print *, "An incidence matrix for the given parameters found:"
     call writeMatrix(incs)
     stop
  enddo

  ! Memory cleanup
  call deconstruct(incs)
end program bibd_ca

subroutine randomCA_BIBD(incs, optSteps)
  use mtmod
  use incidence_structure
  
  implicit none

  type(IncidenceStructure) incs

  integer rowSum, colSum
  integer changeFactor ! Whether the cell will come dormant, or stay dormant
  integer i,j,point
  integer opt_row, opt_col
  integer optSteps
  integer distanceTilGoal, nothing, dot_prod
  logical theEnd, isIn, nOpt

  integer opt_count
  !write (*,*) "DEBUG: randomCA_BIBD"
  theEnd=.false.
  incs%sumTotal=0
  changeFactor=0

  opt_count=0

  nothing=0

  do while(theEnd.eqv..false.)
     !call writeMatrix(incs)
     if (nOpt(optSteps, incs) .eqv. .true.) then
        theEnd=.true.
        print *, "Total count of opt steps:", opt_count
        return
     endif

     !if(nothing>100) then
     !   print *, "Oh deer: ", nothing
     !endif

     changeFactor=0

     i=int(generateRandomNumber()*(incs%v))+1
     j=int(generateRandomNumber()*(incs%b))+1
     colSum=incs%sumInCol(j)
     rowSum=incs%sumInRow(i)

     if (dormant(incs,i,j)) then
        do point=1,incs%v
           if (point==i) cycle
           dot_prod=incs%dp(i,point)
           if (dot_prod<incs%lmbd) changeFactor=changeFactor+1 ! At most v-1           
        enddo
        if (incs%sumTotal<incs%v*incs%r) changeFactor=changeFactor+(incs%v*incs%r-incs%sumTotal)
        if (colSum<incs%k) changeFactor=changeFactor+(incs%k-colSum) ! At most k
        if (rowSum<incs%r) changeFactor=changeFactor+(incs%r-rowSum) ! At most r
        if(int(generateRandomNumber()*(incs%v+incs%k-1+incs%r-1))<changeFactor) then
           call flip(incs,i,j)
           !if(nothing/=0) print *, "Inactive for", nothing, "iterations"
           nothing=0
        else
           nothing=nothing+1
        endif
     else if (active(incs,i,j)) then
        do point=1,incs%v
           if (point==i) cycle
           dot_prod=incs%dp(i,point)
           if (dot_prod>incs%lmbd) changeFactor=changeFactor+1 ! Najviše v-1
        enddo
        if (incs%sumTotal>incs%v*incs%r) changeFactor=changeFactor+(incs%sumTotal-incs%v*incs%r)
        if (colSum>incs%k) changeFactor=changeFactor+(colSum-incs%k) ! Najviše v-k+1?
        if (rowSum>incs%r) changeFactor=changeFactor+(rowSum-incs%r) ! Najviše b-r+1?
        if(int(generateRandomNumber()*(incs%v-incs%k+1+incs%b-incs%r+1))<changeFactor) then
           call flip(incs,i,j)
           !if(nothing/=0) print *, "Inactive for", nothing, "iterations"
           nothing=0
        else
           nothing=nothing+1           
        endif
     else
        print *, "Severe fatal error"
        call exit(1)
     endif
  enddo
end subroutine randomCA_BIBD

recursive logical function nOpt(n,incs) result (successfulOpt)
  use incidence_structure
  implicit none

  type(IncidenceStructure) incs
  integer n,i,j

  successfulOpt=.false.
  if(n<1) then
     successfulOpt=.false.
     return
  else if (abs(incs%sumTotal-incs%v*incs%r)/=n) then
     if(n==1) then
        successfulOpt=.false.
        return
     else
        successfulOpt = nOpt(n-1,incs)
        return
     endif
  endif

  !call writeMatrix(incs)

  do i=1,incs%v
     do j=1,incs%b
        if(& ! If it makes sence to do an n-opt step
             (dormant(incs,i,j) &
             .and. incs%sumTotal==incs%v*incs%r-n &
             .and. (incs%sumInRow(i)==incs%r-n .or. incs%sumInCol(j)==incs%k-n) &
             )&
             .or.&
             (active(incs,i,j) &
             .and. incs%sumTotal==incs%v*incs%r+n &
             .and. (incs%sumInRow(i)==incs%r+n .or. incs%sumInCol(j)==incs%k+n))&
             ) then
           call flip(incs,i,j)
           if(n == 1) then
              if (isBIBD(incs).eqv..True.) then
                 !print *, n, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else if (abs(incs%sumTotal-incs%v*incs%r)==n - 1) then
              !write (*,*) "DEBUG: Calling n-1 opt", n
              if(nOpt(n-1,incs)) then
                 print *, n, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else
              call flip(incs,i,j)
              print *, n, "-opt fail"
              successfulOpt=.false.
              return
           endif
        endif
     enddo
  enddo
  !write (*,*) "DEBUG: end of nOpt", n
  return
end function nOpt
