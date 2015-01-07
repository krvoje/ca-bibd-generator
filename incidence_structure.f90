module incidence_structure

  implicit none

  type IncidenceStructure
     integer, dimension(:,:), allocatable:: incidences
     integer, dimension(:,:), allocatable:: dp
     integer, dimension(:), allocatable:: sumInRow
     integer, dimension(:), allocatable:: sumInCol
     integer v,k,b,lmbd,r,sumTotal     
  end type IncidenceStructure

contains

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

end module incidence_structure
