! A cellular automaton that tries to generate a BIBD
! for the given parameters

program bibd

  implicit none

  logical isBIBD
  real generateRandomNumber

  integer, dimension(:,:), allocatable:: incidences
  integer, dimension(:), allocatable:: sumInRow
  integer, dimension(:), allocatable:: sumInCol
  integer vertices,k,blocks,lmbd,r,i,j,optSteps,sumTotal
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
     read (f,"(I10)") vertices
     call get_command_argument(2,f)
     read (f,"(I10)") k
     call get_command_argument(3,f)
     read (f,"(I10)") lmbd
     call get_command_argument(4,f)
     read (f,"(I10)") optSteps
  endif

  ! Test for neccessary BIBD conditions
  r_r=lmbd*(vertices-1)/(k-1)
  b_r=r_r*vertices/k

  if (r_r/=int(r_r).or.b_r/=int(b_r)) then
     write (*,*) "Invalid BIBD parameters."
     stop
  else
     r=int(r_r)
     blocks=int(b_r)
     print *, "v,k,λ,b,r=", vertices,k,lmbd,blocks,r
  endif

  allocate(incidences(1:vertices,1:blocks))
  allocate(sumInRow(1:vertices))
  allocate(sumInCol(1:blocks))

  call seedRandomGenerator()
  incidences(1:vertices,1:blocks)=int(generateRandomNumber())
  sumInRow(1:vertices)=0
  sumInCol(1:blocks)=0

  ! Check if the random generator generated correctly TODO: This is debug, remove
  do i=1,vertices
     do j=1,blocks
        if(incidences(i,j)==1) then
           sumInRow(i)=sumInRow(i)+1
           sumInCol(j)=sumInCol(j)+1
           sumTotal=sumTotal+1
        endif
        if (incidences(i,j)/=1 .and. incidences(i,j)/=0) then
           print *, "Random generator error, the incidence matrix is:"
           print *, incidences(i,j)
        endif
     enddo
  enddo

  ! Rince and repeat until BIBD
  do while(isBIBD(incidences,vertices,k,lmbd,blocks,r,sumInRow,sumInCol) .eqv. .false.)     
     call randomCA_BIBD(incidences,vertices,k,lmbd,blocks,r,optSteps,sumTotal,sumInRow,sumInCol)
     print *, "An incidence matrix for the given parameters found:"
     call writeMatrix(incidences, vertices, blocks)
     stop
  enddo

  ! Memory cleanup
  deallocate(incidences)
  deallocate(sumInRow)
  deallocate(sumInCol)
end program bibd

subroutine writeMatrix(incidences, v, b)
  integer v, b
  integer incidences(1:v, 1:b)
  write (*,*) "writeMatrix"
  do i=1,v
     do j=1,b
        write (*,"(I1)",advance='no') incidences(i,j)
     enddo
     print *
  enddo  
end subroutine writeMatrix

!!***
!! Checks if the incidence matrix with the other parameters is a BIBD
!! Returns .True. if this is a 2-(v,k,λ) BIBD
!! .False. otherwise
logical function isBIBD(incidences,v,k,lmbd,b,r,sumInRow,sumInCol)

  integer r,b,v,k,lmbd,i,j
  integer incidences(1:v,1:b)
  integer sumInRow(1:v)
  integer sumInCol(1:b)
  !write (*,*) "DEBUG: isBIBD"
 !  call writeMatrix(incidences,v,k)
   
  ! If the incidence matrix has a row with sum non-equal to r, this is not a BIBD
  do i=1,v
     if(sumInRow(i)/=r) then
        isBIBD=.False.
        return
     endif
  enddo

  ! If the incidence matrix has a row with sum non-equal to k, this is not a BIBD
  do i=1,b
     if(sumInCol(i)/=r) then
        isBIBD=.False.
        return
     endif
  enddo

  ! If there are two blocks (rows) whose intersection does not contain lambda vertices,
  ! this is not a BIBD
  do i=1,(v-1)
     do j=(i+1),v
        if (dot_product(incidences(i,:),incidences(j,:))/=lmbd) then
           isBIBD=.False.
           return
        endif
     enddo
  end do

  ! Otherwise, we got a BIBD
  isBIBD=.True.
  write (*,*) "Is a BIBD!"
  return
end function isBIBD

subroutine randomCA_BIBD(incidences,v,k,lmbd,b,r,optSteps,sumTotal,sumInRow,sumInCol)
  use mtmod

  implicit none

  ! Sum of all incidences, computed while running
  integer sumTotal, rowSum, colSum
  integer vitality ! Whether the cell will come alive, or stay dormant
  integer r,b,v,k,lmbd,i,j,point
  integer incidences(1:v,1:b)
  integer sumInRow(1:v)
  integer sumInCol(1:b)  
  integer opt_row, opt_col
  integer optSteps
  integer distanceTilGoal, nothing
  logical theEnd, isBIBD, isIn, nOpt

  real generateRandomNumber

  integer opt_count
  !write (*,*) "DEBUG: randomCA_BIBD"
  theEnd=.false.
  sumTotal=0
  vitality=0

  opt_count=0

  nothing=0

  do while(theEnd.eqv..false.)
     !call writeMatrix(incidences,v,b)
     if (nOpt(optSteps, incidences,v,k,lmbd,b,r,sumTotal,sumInRow,sumInCol) .eqv. .true.) then
        theEnd=.true.
        print *, "Total count of opt steps:", opt_count
        return
     endif

     vitality=0

     i=int(generateRandomNumber()*(v))+1
     j=int(generateRandomNumber()*(b))+1
     colSum=sumInCol(j)
     rowSum=sumInRow(i)

     if (incidences(i,j)==0) then
        do point=1,v
           if (point==i) cycle
           if (dot_product(incidences(i,:),incidences(point,:))<lmbd) vitality=vitality+1 ! At most v-1           
        enddo
        if (sumTotal<v*r) vitality=vitality+1 ! At most 1
        if (colSum<k) vitality=vitality+(k-colSum) ! At most k
        if (rowSum<r) vitality=vitality+(r-rowSum) ! At most r
        if(int(generateRandomNumber()*(v+k-1+r-1))<vitality) then
           call flip(i,j,incidences,v,b,sumTotal,sumInRow,sumInCol)
           nothing=0
        else
           nothing=nothing+1
           !if(nothing>optSteps)
           !print *, "nothing done, vitality: ", nothing, vitality
        endif
     else if (incidences(i,j)==1) then
        do point=1,v
           if (point==i) cycle
           if (dot_product(incidences(i,:),incidences(point,:))>lmbd) vitality=vitality+1 ! Najviše v-1
        enddo
        if (sumTotal>v*r) vitality=vitality+1 ! Najviše 1
        if (colSum>k) vitality=vitality+(colSum-k) ! Najviše v-k+1?
        if (rowSum>r) vitality=vitality+(rowSum-r) ! Najviše b-r+1?
        if(int(generateRandomNumber()*(v+v-k+1+b-r+1))<vitality) then
           call flip(i,j,incidences,v,b,sumTotal,sumInRow,sumInCol)
           nothing=0
        else
           nothing=nothing+1
           !if(nothing>optSteps)
           !print *, "nothing done, vitality: ", nothing, vitality
        endif
     else
        print *, "Severe fatal error"
        call exit(1)
     endif     
  enddo
end subroutine randomCA_BIBD

subroutine seedRandomGenerator()
  use mtmod

  external ZBQLINI
  !!write (*,*) "DEBUG: seedRandomGenerator"
  call ZBQLINI(time())
  call sgrnd(time())
  call srand(time())
end subroutine seedRandomGenerator

real function generateRandomNumber()
  use mtmod

  double precision ZBQLU01,ZBQLUAB, ZBQLEXP
  integer i
  !!write (*,*) "DEBUG: generateRandomNumber"
  !generateRandomNumber=ZBQLU01(i)
  !generateRandomNumber=grnd() !<- Ova implementacija Mersenne twistera je koma za ovu svrhu
  generateRandomNumber=rand()
  if (generateRandomNumber==1) generateRandomNumber=0 ! Računamo da se nikad ne dobije 1, jer bi to potrgalo sve
  return

end function generateRandomNumber

recursive logical function nOpt(n,incidences,v,k,lmbd,b,r,sumTotal,sumInRow,sumInCol) result (successfulOpt)
  implicit none

  integer n
  integer incidences(1:v,1:b)
  integer sumInRow(1:v)
  integer sumInCol(1:b)
  integer v,k,lmbd,b,r,sumTotal,i,j
  logical isBIBD

  successfulOpt=.false.
  if(n<1) then
     successfulOpt=.false.
     return
  else if (abs(sumTotal-v*r)/=n) then
     if(n==1) then
        successfulOpt=.false.
        return
     else
        successfulOpt = nOpt(n-1,incidences,v,k,lmbd,b,r,sumTotal,sumInRow,sumInCol)
        return
     endif
  endif

  !call writeMatrix(incidences,v,b)
  
  do i=1,v
     do j=1,b
        if(& ! If it makes sence to do an n-opt step
             (incidences(i,j)==0 .and. sumTotal==v*r-n .and. (sumInRow(i)==r-n .or. sumInCol(j)==k-n) )&
             .or.&
             (incidences(i,j)==1 .and. sumTotal==v*r+n .and. (sumInRow(i)==r+n .or. sumInCol(j)==k+n) )&
             ) then
           call flip(i,j,incidences,v,b,sumTotal,sumInRow,sumInCol)
           if(n == 1) then
              if (isBIBD(incidences,v,k,lmbd,b,r,sumInRow,sumInCol).eqv..True.) then
                 !print *, n, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else if (abs(sumTotal-v*r)==n - 1) then
              !write (*,*) "DEBUG: Calling n-1 opt", n
              if(nOpt(n-1,incidences,v,k,lmbd,b,r,sumTotal,sumInRow,sumInCol)) then
                 print *, n, "-opt yielded a BIBD"
                 successfulOpt=.true.
                 return
              endif
           else
              call flip(i,j,incidences,v,b,sumTotal,sumInRow,sumInCol)
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

! Flips a CA's cell value, and updates the sum
subroutine flip (i,j,incidences,v,b,sumTotal,sumInRow,sumInCol)
  integer sumTotal, v, b, i, j
  integer incidences(1:v,1:b)
  integer sumInRow(1:v)
  integer sumInCol(1:b)
  !!write (*,*) "DEBUG: flip"
  incidences(i,j)=abs(incidences(i,j)-1)
  if (incidences(i,j)==0) then
     sumTotal=sumTotal-1
     sumInRow(i)=sumInRow(i)-1
     sumInCol(j)=sumInCol(j)-1
  endif
  if (incidences(i,j)==1) then
     sumTotal=sumTotal+1
     sumInRow(i)=sumInRow(i)+1
     sumInCol(j)=sumInCol(j)+1
  endif
end subroutine flip
