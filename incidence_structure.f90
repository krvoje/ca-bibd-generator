module incidence_structure
  use utils
  use mtmod
  
  implicit none

  type IncidenceStructure
     integer, dimension(:,:), allocatable:: INCIDENCES ! incidences
     integer, dimension(:,:), allocatable:: BLOCK_INTERSECTION ! dp
     integer, dimension(:), allocatable:: SUM_IN_ROW ! SUM_IN_ROW
     integer, dimension(:), allocatable:: SUM_IN_COL ! SUM_IN_COL
     integer VERTICES !v
     integer INCIDENCES_PER_VERTICE!k
     integer BLOCKS !b
     integer LAMBDA!lmbd
     
     ! Heuristic helper vals
     integer VERTICES_PER_BLOCK!r
     integer SUM_TOTAL     
     integer VERTICES_X_BLOCKS     
     integer SUM_IDEAL ! Should be VERTICES_PER_BLOC*VERTICES
     logical SUM_TOTAL_LESS_THAN_IDEAL
     logical SUM_TOTAL_MORE_THAN_IDEAL
     logical SUM_TOTAL_NOT_IDEAL

     integer MAX_CHANGE_FACTOR
     integer MAX_CHANGE_FACTOR_DORMANT
     integer MAX_CHANGE_FACTOR_ACTIVE

     integer CHANGE_FACTOR_ALIVE_INCREMENT
  end type IncidenceStructure

contains

  subroutine construct(IS,v,k,lmbd)
    type(IncidenceStructure) IS    
    integer v,k,lmbd
    real r_r, b_r
    integer i,j

    IS%VERTICES=v
    IS%INCIDENCES_PER_VERTICE=k
    IS%LAMBDA=lmbd

    ! Test for neccessary BIBD conditions
    r_r=IS%LAMBDA*(IS%VERTICES-1)/(IS%INCIDENCES_PER_VERTICE-1)
    b_r=r_r*IS%VERTICES/IS%INCIDENCES_PER_VERTICE

    if (r_r/=int(r_r).or.b_r/=int(b_r)) then
       write (*,*) "Invalid BIBD parameters."       
       stop
    else
       IS%VERTICES_PER_BLOCK=int(r_r)
       IS%BLOCKS=int(b_r)
       print *, "v,k,λ,b,r=", IS%VERTICES,IS%INCIDENCES_PER_VERTICE,IS%LAMBDA,IS%BLOCKS,IS%VERTICES_PER_BLOCK
    endif

    allocate(IS%incidences(1:IS%VERTICES,1:IS%BLOCKS))
    allocate(IS%BLOCK_INTERSECTION(1:IS%VERTICES,1:IS%VERTICES))
    allocate(IS%SUM_IN_ROW(1:IS%VERTICES))
    allocate(IS%SUM_IN_COL(1:IS%BLOCKS))

    IS%incidences(1:IS%VERTICES,1:IS%BLOCKS)=0
    IS%BLOCK_INTERSECTION(1:IS%VERTICES,1:IS%VERTICES)=0
    IS%SUM_IN_ROW(1:IS%VERTICES)=0
    IS%SUM_IN_COL(1:IS%BLOCKS)=0
    IS%SUM_TOTAL=0

    IS%SUM_IDEAL=IS%VERTICES_PER_BLOCK*IS%VERTICES
    
    IS%SUM_TOTAL_LESS_THAN_IDEAL=IS%SUM_TOTAL<=IS%SUM_IDEAL
    IS%SUM_TOTAL_MORE_THAN_IDEAL=IS%SUM_TOTAL>=IS%SUM_IDEAL
    IS%SUM_TOTAL_NOT_IDEAL = (IS%SUM_TOTAL /= IS%SUM_IDEAL)

    IS%VERTICES_X_BLOCKS = IS%VERTICES * IS%BLOCKS

    IS%MAX_CHANGE_FACTOR_DORMANT=(IS%VERTICES-1) &
         + abs(IS%SUM_IDEAL - IS%SUM_TOTAL) &
         + IS%INCIDENCES_PER_VERTICE -1 &
         + (IS%VERTICES_PER_BLOCK-1)

    IS%MAX_CHANGE_FACTOR_ACTIVE=IS%CHANGE_FACTOR_ALIVE_INCREMENT * (IS%VERTICES-1) &
         + IS%BLOCKS &
         + IS%VERTICES
    IS%MAX_CHANGE_FACTOR = IS%MAX_CHANGE_FACTOR_DORMANT + IS%MAX_CHANGE_FACTOR_ACTIVE

    IS%CHANGE_FACTOR_ALIVE_INCREMENT = max(IS%INCIDENCES_PER_VERTICE, IS%VERTICES_PER_BLOCK)
  end subroutine construct

  subroutine updateCache(IS)
    type(IncidenceStructure) IS

    integer i,j

    do i=1,IS%VERTICES
       IS%SUM_IN_ROW(i)=sum(IS%incidences(i,:))
    enddo

    do j=1,IS%BLOCKS
       IS%SUM_IN_COL(i)=sum(IS%incidences(:,j))
    enddo
    
    do i=1,IS%VERTICES
       do j=1,IS%VERTICES
          IS%BLOCK_INTERSECTION(i,j)=dot_product(IS%incidences(i,:), IS%incidences(j,:))
       enddo
    enddo

    IS%SUM_TOTAL=sum(IS%SUM_IN_ROW(:))  
    
    IS%SUM_TOTAL_LESS_THAN_IDEAL = IS%SUM_TOTAL <= IS%SUM_IDEAL
    IS%SUM_TOTAL_MORE_THAN_IDEAL = IS%SUM_TOTAL >= IS%SUM_IDEAL
    IS%SUM_TOTAL_NOT_IDEAL = (IS%SUM_TOTAL /= IS%SUM_IDEAL)
  end subroutine updateCache
  
  subroutine deconstruct(IS)
    type(IncidenceStructure) IS
    deallocate(IS%incidences)
    deallocate(IS%BLOCK_INTERSECTION)
    deallocate(IS%SUM_IN_ROW)
    deallocate(IS%SUM_IN_COL)
  end subroutine deconstruct

  logical function active(IS,row,col)
    type(IncidenceStructure) IS
    integer row,col

    if(IS%incidences(row,col)>0) then
       active=.True.
    else
       active=.False.
    endif

    return
  end function active

  logical function dormant(IS,row,col)
    type(IncidenceStructure) IS
    integer row,col

    if(IS%incidences(row,col)>0) then
       dormant=.False.
    else
       dormant=.true.
    endif

    return
  end function dormant

  ! Flips a CA's cell value, and updates the sums and products
  subroutine flip (IS,row,col)

    type(IncidenceStructure) IS
    integer row,col
    integer otherRow,newVal
    !!write (*,*) "DEBUG: flip"
    newVal=abs(IS%incidences(row,col)-1)
    IS%incidences(row,col)=newVal
    if (newVal==0) then
       IS%SUM_TOTAL=IS%SUM_TOTAL-1

       IS%SUM_IN_ROW(row)=IS%SUM_IN_ROW(row)-1
       IS%SUM_IN_COL(col)=IS%SUM_IN_COL(col)-1
       IS%BLOCK_INTERSECTION(row,row)=IS%BLOCK_INTERSECTION(row,row)-1
    endif
    if (newVal==1) then
       IS%SUM_TOTAL=IS%SUM_TOTAL+1
       IS%SUM_IN_ROW(row)=IS%SUM_IN_ROW(row)+1
       IS%SUM_IN_COL(col)=IS%SUM_IN_COL(col)+1
       IS%BLOCK_INTERSECTION(row,row)=IS%BLOCK_INTERSECTION(row,row)+1
    endif

    do otherRow=1,IS%VERTICES
       if(otherRow==row) cycle
       if(newVal==0 .and. IS%incidences(otherRow,col)==1) then
          IS%BLOCK_INTERSECTION(otherRow,row)=IS%BLOCK_INTERSECTION(otherRow,row)-1
          IS%BLOCK_INTERSECTION(row,otherRow)=IS%BLOCK_INTERSECTION(row,otherRow)-1
       endif
       if(newVal==1 .and. IS%incidences(otherRow,col)==1) then
          IS%BLOCK_INTERSECTION(otherRow,row)=IS%BLOCK_INTERSECTION(otherRow,row)+1
          IS%BLOCK_INTERSECTION(row,otherRow)=IS%BLOCK_INTERSECTION(row,otherRow)+1
       endif
    enddo
  end subroutine flip

  !!***
  !! Checks if the incidence matrix with the other parameters is a BIBD
  !! Returns .True. if this is a 2-(v,k,λ) BIBD
  !! .False. otherwise
  logical function isBIBD(IS)

    type(IncidenceStructure) IS

    integer i,j
    
    !call writeMatrix(IS)

    if(IS%SUM_TOTAL /= IS%SUM_IDEAL) then
       isBIBD=.False.
       return
    endif
    
    ! If the incidence matrix has a row with sum non-equal to r, this is not a BIBD
    do i=1,IS%VERTICES
       if(IS%SUM_IN_ROW(i)/=IS%INCIDENCES_PER_VERTICE) then
          isBIBD=.False.
          return
       endif
    enddo

    ! If the incidence matrix has a row with sum non-equal to k, this is not a BIBD
    do i=1,IS%BLOCKS
       if(IS%SUM_IN_COL(i)/=IS%VERTICES_PER_BLOCK) then
          isBIBD=.False.
          return
       endif
    enddo

    ! If there are two IS%BLOCKS (rows) whose intersection does not contain lambda IS%VERTICES,
    ! this is not a BIBD
    do i=1,(IS%VERTICES-1)
       do j=(i+1),IS%VERTICES
          if (IS%BLOCK_INTERSECTION(i,j)/=IS%LAMBDA) then
             isBIBD=.False.
             return
          endif
       enddo
    end do

    ! Otherwise, we got a BIBD
    isBIBD=.True.
    return
  end function isBIBD

  subroutine writeMatrix(IS)
    type(IncidenceStructure) IS
    integer i,j

    write(*,*),"["
    do j=1,IS%BLOCKS
       do i=1,IS%VERTICES
          write (*,"(I1)",advance='no') IS%incidences(i,j)
       enddo
       write (*,"(A1)",advance='no') "-"
       write (*,"(I1)",advance='no') IS%SUM_IN_COL(j)
       print *
    enddo
    do i=1,IS%VERTICES
       write (*,"(A1)",advance='no') "|"
    enddo
    print *
    do i=1,IS%VERTICES
       write (*,"(I1)",advance='no') IS%SUM_IN_ROW(i)
    enddo
    print *
    write(*,*),"]"
  end subroutine writeMatrix

end module incidence_structure
