module incidence_structure
  use utils
  
  implicit none

  type IncidenceStructure
     integer, dimension(:,:), allocatable:: incidences
     integer, dimension(:,:), allocatable:: dp
     integer, dimension(:), allocatable:: sumInRow
     integer, dimension(:), allocatable:: sumInCol
     integer v,k,b,lmbd,r,sumTotal     
  end type IncidenceStructure

contains

  subroutine construct(incs,v,k,lmbd)
    type(IncidenceStructure) incs
    integer v,k,lmbd
    real r_r, b_r
    integer i,j

    call seedRandomGenerator()
    
    incs%v=v
    incs%k=k
    incs%lmbd=lmbd

    ! Test for neccessary BIBD conditions
    r_r=incs%lmbd*(incs%v-1)/(incs%k-1)
    b_r=r_r*incs%v/incs%k

    if (r_r/=int(r_r).or.b_r/=int(b_r)) then
       write (*,*) "Invalid BIBD parameters."       
       stop
    else
       incs%r=int(r_r)
       incs%b=int(b_r)
       print *, "v,k,λ,b,r=", incs%v,incs%k,incs%lmbd,incs%b,incs%r
    endif

    allocate(incs%incidences(1:incs%v,1:incs%b))
    allocate(incs%dp(1:incs%v,1:incs%v))
    allocate(incs%sumInRow(1:incs%v))
    allocate(incs%sumInCol(1:incs%b))
    
    do i=1,incs%v
       do j=1,incs%b
          incs%incidences(i,j)=int(generateRandomNumber())
          if(active(incs,i,j)) then
             incs%sumInRow(i)=incs%sumInRow(i)+1
             incs%sumInCol(j)=incs%sumInCol(j)+1
             incs%sumTotal=incs%sumTotal+1
          endif
          if (incs%incidences(i,j)/=1 .and. incs%incidences(i,j)/=0) then
             print *, "Random generator error, the incidence matrix is:"
             print *, incs%incidences(i,j)
          endif
       enddo
    enddo

    do i=1,incs%v
       do j=1,incs%v
          incs%dp(i,j)=dot_product(incs%incidences(i,:), incs%incidences(j,:))
       enddo
    enddo

    call writeMatrix(incs)
  end subroutine construct

  subroutine deconstruct(incs)
    type(IncidenceStructure) incs
    deallocate(incs%incidences)
    deallocate(incs%dp)
    deallocate(incs%sumInRow)
    deallocate(incs%sumInCol)
  end subroutine deconstruct

  logical function active(incs,row,col)
    type(IncidenceStructure) incs
    integer row,col

    if(incs%incidences(row,col)>0) then
       active=.True.
    else
       active=.False.
    endif

    return
  end function active

  logical function dormant(incs,row,col)
    type(IncidenceStructure) incs
    integer row,col

    if(incs%incidences(row,col)>0) then
       dormant=.False.
    else
       dormant=.true.
    endif

    return
  end function dormant

  ! Flips a CA's cell value, and updates the sums and products
  subroutine flip (incs,row,col)

    type(IncidenceStructure) incs
    integer row,col
    integer otherRow,newVal
    !!write (*,*) "DEBUG: flip"
    newVal=abs(incs%incidences(row,col)-1)
    incs%incidences(row,col)=newVal
    if (newVal==0) then
       incs%sumTotal=incs%sumTotal-1

       incs%sumInRow(row)=incs%sumInRow(row)-1
       incs%sumInCol(col)=incs%sumInCol(col)-1
       incs%dp(row,row)=incs%dp(row,row)-1
    endif
    if (newVal==1) then
       incs%sumTotal=incs%sumTotal+1
       incs%sumInRow(row)=incs%sumInRow(row)+1
       incs%sumInCol(col)=incs%sumInCol(col)+1
       incs%dp(row,row)=incs%dp(row,row)+1
    endif

    do otherRow=1,incs%v
       if(otherRow==row) cycle
       if(newVal==0 .and. incs%incidences(otherRow,col)==1) then
          incs%dp(otherRow,row)=incs%dp(otherRow,row)-1
          incs%dp(row,otherRow)=incs%dp(row,otherRow)-1
       endif
       if(newVal==1 .and. incs%incidences(otherRow,col)==1) then
          incs%dp(otherRow,row)=incs%dp(otherRow,row)+1
          incs%dp(row,otherRow)=incs%dp(row,otherRow)+1
       endif
    enddo
  end subroutine flip

  !!***
  !! Checks if the incidence matrix with the other parameters is a BIBD
  !! Returns .True. if this is a 2-(v,k,λ) BIBD
  !! .False. otherwise
  logical function isBIBD(incs)

    type(IncidenceStructure) incs

    integer i,j

    !write (*,*) "DEBUG: isBIBD"
    !  call writeMatrix(incs%incidences,v,k)

    ! If the incidence matrix has a row with sum non-equal to r, this is not a BIBD
    do i=1,incs%v
       if(incs%sumInRow(i)/=incs%r) then
          isBIBD=.False.
          return
       endif
    enddo

    ! If the incidence matrix has a row with sum non-equal to k, this is not a BIBD
    do i=1,incs%b
       if(incs%sumInCol(i)/=incs%r) then
          isBIBD=.False.
          return
       endif
    enddo

    ! If there are two incs%b (rows) whose intersection does not contain lambda incs%v,
    ! this is not a BIBD
    do i=1,(incs%v-1)
       do j=(i+1),incs%v
          if (incs%dp(i,j)/=incs%lmbd) then
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

  subroutine writeMatrix(incs)
    type(IncidenceStructure) incs
    integer i,j

    write (*,*) "writeMatrix"
    do i=1,incs%v
       do j=1,incs%b
          write (*,"(I1)",advance='no') incs%incidences(i,j)
       enddo
       print *
    enddo
  end subroutine writeMatrix

end module incidence_structure
