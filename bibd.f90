! A cellular automaton that tries to generate a BIBD
! for the given parameters

program bibd

  implicit none

  logical isBIBD
  real generateRandomNumber

  integer, dimension(:,:), allocatable:: incidences
  integer vertices,k,blocks,lmbd,r,i,j
  real r_r,b_r
  character(100) f

  ! If the number of command line params is invalid, print instructions
  if (command_argument_count() /= 3) then
     write (*,*) "Usage:"
     call get_command_argument(0,f)
     write (*,*) f ,"v k λ"
     stop        
  else
     call get_command_argument(1,f)
     read (f,"(I10)") vertices
     call get_command_argument(2,f)
     read (f,"(I10)") k
     call get_command_argument(3,f)
     read (f,"(I10)") lmbd
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

  call seedRandomGenerator()
  incidences(1:vertices,1:blocks)=int(generateRandomNumber())

  ! Check if the random generator generated correctly TODO: This is debug, remove
  do i=1,vertices
     do j=1,blocks
        if (incidences(i,j)/=1 .and. incidences(i,j)/=0) then
           print *, "Random generator error, the incidence matrix is:"
           print *, incidences(i,j)
        endif
     enddo
  enddo

  ! Rince and repeat until BIBD
  do while(isBIBD(incidences,vertices,k,lmbd,blocks,r) .eqv. .false.)
     call randomCA_BIBD(incidences,vertices,k,lmbd,blocks,r)
     print *, "An incidence matrix for the given parameters found:"
     do i=1,vertices
        do j=1,blocks
           write (*,"(I1)",advance='no') incidences(i,j)
        enddo
        print *
     enddo
     stop
  enddo

  ! Memory cleanup
  deallocate(incidences)
end program bibd

!!***
!! Checks if the incidence matrix with the other parameters is a BIBD
!! Returns .True. if this is a 2-(v,k,λ) BIBD
!! .False. otherwise
logical function isBIBD(incidences,v,k,lmbd,b,r)
  integer r,b,v,k,lmbd,i,j
  integer incidences(1:v,1:b)

  ! If the incidence matrix has a row with sum non-equal to r, this is not a BIBD
  do i=1,v
     if (sum(incidences(i,:))/=r) then
        isBIBD=.False.
        return
     endif
  enddo

  ! If the incidence matrix has a row with sum non-equal to k, this is not a BIBD
  do i=1,b
     if (sum(incidences(:,i))/=k) then
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

subroutine randomCA_BIBD(incidences,v,k,lmbd,b,r)        
  use mtmod

  implicit none

  ! Sum of all incidences, computed while running
  integer sum_total, sum_rows, sum_cols
  integer cell_liveness_factor ! Whether the cell will come alive, or stay dormant
  integer r,b,v,k,lmbd,i,j,point
  integer incidences(1:v,1:b)
  integer opt_row, opt_col
  logical theEnd, isBIBD, oneOpt, twoOpt, isIn

  real generateRandomNumber

  integer opt_count

  theEnd=.false.
  sum_total=0
  cell_liveness_factor=0

  opt_count=0

  do i=1,v
     sum_total = sum_total + sum(incidences(i,:))
  enddo

  i=int(generateRandomNumber()*(v-1))+1
  j=int(generateRandomNumber()*(b-1))+1

  do while(theEnd.eqv..false.)
     ! If it makes sence to try a 2-opt
     if (abs(sum_total-v*r)==2) then
        opt_count=opt_count+v*b
        if (twoOpt(incidences,v,k,lmbd,b,r,sum_total) .eqv. .true.) then
           theEnd=.true.
           print *, "Total count of 2-opt steps:", opt_count
           return
        endif
     endif

     ! If it makes sence to try a 1-opt
     if (abs(sum_total-v*r)==1) then  !ako opće ima smisla raditi 1-opt korak
        opt_count=opt_count+1
        if (oneOpt(incidences,v,k,lmbd,b,r,sum_total).eqv..true.) then
           theEnd=.true.
           print *, "Total count of 1-opt steps:", opt_count
           return
        endif
     endif

     cell_liveness_factor=0

     i=int(generateRandomNumber()*(v))+1
     j=int(generateRandomNumber()*(b))+1
     sum_cols=sum(incidences(:,j))
     sum_rows=sum(incidences(i,:))

     ! If the cell is dormant
     if (incidences(i,j)==0) then
        do point=1,v
           if (point==i) cycle
           if (dot_product(incidences(i,:),incidences(point,:))<lmbd) cell_liveness_factor=cell_liveness_factor+1 ! Najviše v-1
        enddo
        if (sum_total<v*r) cell_liveness_factor=cell_liveness_factor+1 ! Najviše 1
        if (sum_cols<k) cell_liveness_factor=cell_liveness_factor+(k-sum_cols) ! Najviše k
        if (sum_rows<r) cell_liveness_factor=cell_liveness_factor+(r-sum_rows) ! Najviše r
        if(int(generateRandomNumber()*(v+k-1+r-1))<cell_liveness_factor) call flip(incidences(i,j),sum_total)
        ! If the cell is active
     else if (incidences(i,j)==1) then
        do point=1,v
           if (point==i) cycle
           if (dot_product(incidences(i,:),incidences(point,:))>lmbd) cell_liveness_factor=cell_liveness_factor+1 ! Najviše v-1
        enddo
        if (sum_total>v*r) cell_liveness_factor=cell_liveness_factor+1 ! Najviše 1
        if (sum_cols>k) cell_liveness_factor=cell_liveness_factor+(sum_cols-k) ! Najviše v-k+1?
        if (sum_rows>r) cell_liveness_factor=cell_liveness_factor+(sum_rows-r) ! Najviše b-r+1?
        if(int(generateRandomNumber()*(v+v-k+1+b-r+1))<cell_liveness_factor) call flip(incidences(i,j),sum_total)! Eksperimenti, dat prednost živim incidencijama
     endif
  enddo
end subroutine randomCA_BIBD

subroutine seedRandomGenerator()
  use mtmod

  external ZBQLINI

  call ZBQLINI(time())
  call sgrnd(time())
end subroutine seedRandomGenerator

real function generateRandomNumber()
  use mtmod

  double precision ZBQLU01,ZBQLUAB, ZBQLEXP
  integer i

  !generateRandomNumber=ZBQLU01(i)
  generateRandomNumber=grnd() !<- Ova implementacija Mersenne twistera je koma za ovu svrhu
  if (generateRandomNumber==1) generateRandomNumber=0 ! Računamo da se nikad ne dobije 1, jer bi to potrgalo sve
  return

end function generateRandomNumber

logical function oneOpt(incidences,v,k,lmbd,b,r,sum_total)
  implicit none

  integer v,k,lmbd,b,r,sum_total,i,j
  integer incidences(1:v,1:b)

  logical isBIBD

  oneOpt=.false.
  !print *, "Starting 1-opt..."
  do i=1,v
     do j=1,b
        if(&
             (incidences(i,j)==0 .and. sum_total==v*r-1 .and. sum(incidences(i,:))==r-1 .and. sum(incidences(:,j))==k-1 )&
             .or.&
             (incidences(i,j)==1 .and. sum_total==v*r+1 .and. sum(incidences(i,:))==r+1 .and. sum(incidences(:,j))==k+1) ) then
           call flip(incidences(i,j),sum_total)
           if (isBIBD(incidences,v,k,lmbd,b,r).eqv..True.) then
              !write (*,*) "1-opt yielded a BIBD"
              oneOpt=.true.
              return
           else
              !write (*,*) "1-opt yielded no BIBD"
              call flip(incidences(i,j),sum_total)
              oneOpt=.false.
           endif
        endif
     enddo
  enddo
  return
end function oneOpt


logical function twoOpt(incidences,v,k,lmbd,b,r,sum_total)
  implicit none

  integer v,k,lmbd,b,r,sum_total,i,j
  integer incidences(1:v,1:b)
  logical isBIBD, oneOpt

  twoOpt=.false.
  !print *, "2-opt!"
  do i=1,v
     do j=1,b
        if(&
             (incidences(i,j)==0 .and. sum_total==v*r-2 .and. (sum(incidences(i,:))==r-2 .or. sum(incidences(:,j))==k-2) )&
             .or.&
             (incidences(i,j)==1 .and. sum_total==v*r+2 .and. (sum(incidences(i,:))==r+2 .or. sum(incidences(:,j))==k+2) )&
             ) then
           call flip(incidences(i,j),sum_total)
           if (abs(sum_total-v*r)==1) then  !ako sada ima smisla raditi 1-opt korak
              if(oneOpt(incidences,v,k,lmbd,b,r,sum_total).eqv..true.) then
                 print *, "2-opt yielded a BIBD"
                 twoOpt=.true.
                 return
              endif
           else
              call flip(incidences(i,j),sum_total)
              twoOpt=.false.
           endif
        endif
     enddo
  enddo
  return
end function twoOpt

! Flips a CA's cell value, and updates the sum
subroutine flip (cell, sum_total)
  integer cell,sum_total

  cell=abs(cell-1)
  if (cell==0) sum_total=sum_total-1
  if (cell==1) sum_total=sum_total+1

end subroutine flip
