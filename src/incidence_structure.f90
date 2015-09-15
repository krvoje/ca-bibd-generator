module incidence_structure
  use utils
  use mtmod
  
  implicit none

  type IncidenceStructure
     integer, dimension(:,:), allocatable:: INCIDENCES ! incidences
     integer, dimension(:,:), allocatable:: ROW_INTERSECTION ! dp
     integer, dimension(:,:), allocatable:: COL_INTERSECTION ! dp
     integer, dimension(:), allocatable:: SUM_IN_ROW ! SUM_IN_ROW
     integer, dimension(:), allocatable:: SUM_IN_COL ! SUM_IN_COL
     integer v !v
     integer r!r
     integer k!k
     integer b !b
     integer LAMBDA!lmbd
     
     ! Heuristic helper vals
     integer SUM_TOTAL     
     integer SUM_IDEAL ! Should be VERTICES_PER_BLOC*VERTICES

     integer heuristic_distance
     integer max_hd_coef
     integer generations
  end type IncidenceStructure

contains

  subroutine construct(IS,v,k,lmbd)
    type(IncidenceStructure) IS    
    integer v,k,lmbd
    real r_r, b_r
    integer i,j

    IS%v=v
    IS%k=k
    IS%LAMBDA=lmbd

    ! Test for neccessary BIBD conditions
    r_r=IS%LAMBDA*(IS%v-1)/(IS%k-1)
    b_r=r_r*IS%v/IS%k

    if (r_r/=int(r_r).or.b_r/=int(b_r)) then
       write (*,*) "Invalid BIBD parameters."       
       stop
    else
       IS%r=int(r_r)
       IS%b=int(b_r)
       print *, "v,k,λ,b,r=", IS%v,IS%k,IS%LAMBDA,IS%b,IS%r       
    endif

    allocate(IS%incidences(1:IS%v,1:IS%b))
    allocate(IS%ROW_INTERSECTION(1:IS%v,1:IS%v))
    allocate(IS%COL_INTERSECTION(1:IS%b,1:IS%b))
    allocate(IS%SUM_IN_ROW(1:IS%v))
    allocate(IS%SUM_IN_COL(1:IS%b))

    IS%incidences(1:IS%v,1:IS%b)=0
    IS%ROW_INTERSECTION(1:IS%v,1:IS%v)=0
    IS%COL_INTERSECTION(1:IS%b,1:IS%b)=0
    IS%SUM_IN_ROW(1:IS%v)=0
    IS%SUM_IN_COL(1:IS%b)=0
    IS%SUM_TOTAL=0

    IS%SUM_IDEAL=IS%k*IS%b
    is%generations=0
    
    call updateCache(IS)
  end subroutine construct

  subroutine updateCache(IS)
    type(IncidenceStructure) IS

    integer i,j

    do i=1,IS%v
       IS%SUM_IN_ROW(i)=sum(IS%incidences(i,:))
    enddo

    do j=1,IS%b
       IS%SUM_IN_COL(j)=sum(IS%incidences(:,j))
    enddo
    
    do i=1,IS%v
       do j=1,IS%v
          IS%ROW_INTERSECTION(i,j)=dot_product(IS%incidences(i,:), IS%incidences(j,:))
       enddo
    enddo

    do i=1,IS%b
       do j=1,IS%b
          IS%COL_INTERSECTION(i,j)=dot_product(IS%incidences(:,i), IS%incidences(:,j))
       enddo
    enddo

    IS%SUM_TOTAL=sum(IS%SUM_IN_ROW(1:IS%v))
    
    call calculateHeuristicDistance(IS)
  end subroutine updateCache
  
  subroutine deconstruct(IS)
    type(IncidenceStructure) IS
    deallocate(IS%incidences)
    deallocate(IS%ROW_INTERSECTION)
    deallocate(IS%COL_INTERSECTION)
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
    integer otherRow,otherCol,newVal
    !!write (*,*) "DEBUG: flip"
    newVal=abs(IS%incidences(row,col)-1)
    IS%incidences(row,col)=newVal
    if (newVal==0) then
       call decrement(IS%SUM_TOTAL,1)
       call decrement(IS%SUM_IN_ROW(row),1)
       call decrement(IS%SUM_IN_COL(col),1)
       call decrement(IS%ROW_INTERSECTION(row,row),1)
       call decrement(IS%COL_INTERSECTION(col,col),1)
    endif
    if (newVal==1) then
       call increment(IS%SUM_TOTAL,1)
       call increment(IS%SUM_IN_ROW(row),1)
       call increment(IS%SUM_IN_COL(col),1)
       call increment(IS%ROW_INTERSECTION(row,row),1)
       call increment(IS%COL_INTERSECTION(col,col),1)
    endif

    do otherRow=1,IS%v
       if(IS%incidences(otherRow,col)==1) then 
          if(newVal==0) then
             call decrement(IS%ROW_INTERSECTION(otherRow,row),1)
             call decrement(IS%ROW_INTERSECTION(row,otherRow),1)
          endif
          if(newVal==1) then
             call increment(IS%ROW_INTERSECTION(otherRow,row),1)
             call increment(IS%ROW_INTERSECTION(row,otherRow),1)
          endif
       endif
    enddo

    do otherCol=1,IS%b
       if(IS%incidences(row,otherCol)==1) then 
          if(newVal==0) then
             call decrement(IS%COL_INTERSECTION(otherCol,col),1)
             call decrement(IS%COL_INTERSECTION(col,otherCol),1)
          endif
          if(newVal==1) then
             call increment(IS%COL_INTERSECTION(otherCol,col),1)
             call increment(IS%COL_INTERSECTION(col,otherCol),1)
          endif
       endif
    enddo

    call calculateHeuristicDistance(IS)
  end subroutine flip

  subroutine calculateHeuristicDistance(IS)
    type(IncidenceStructure) IS

    integer i,j

    IS%HEURISTIC_DISTANCE = (abs(is%sum_total - is%sum_ideal))
    do i=1,IS%v
       call increment(IS%HEURISTIC_DISTANCE, abs(IS%SUM_IN_ROW(i) - IS%r))
    enddo
    do i=1,IS%b
       call increment(IS%HEURISTIC_DISTANCE, abs(IS%SUM_IN_COL(i) - IS%k))
    enddo
    do i=1,IS%v
       do j=1,IS%v
          if(i==j) cycle
          call increment(IS%HEURISTIC_DISTANCE, abs(IS%ROW_INTERSECTION(i,j) - IS%LAMBDA))
       enddo
    enddo
    IS%max_hd_coef = is%SUM_TOTAL + (IS%v - 1) + 2
  end subroutine calculateHeuristicDistance
  
  !!***
  !! Checks if the incidence matrix with the other parameters is a BIBD
  !! Returns .True. if this is a 2-(v,k,λ) BIBD
  !! .False. otherwise
  logical function isBIBD(IS)

    type(IncidenceStructure) IS

    integer i,j

    ! if(is%heuristic_distance <= is%max_hd_coef ) then
    !     print *, "hd: ", is%max_hd_coef, "/", is%heuristic_distance
    !     call writeMatrix(IS)
    ! endif
    
    if(IS%SUM_TOTAL /= IS%SUM_IDEAL) then
       isBIBD=.False.
       !print *, "Sum is non-ideal", is%sum_total, is%sum_ideal
       return
    endif
    
    ! If the incidence matrix has a row with sum non-equal to r, this is not a BIBD
    do i=1,IS%v
       if(IS%SUM_IN_ROW(i)/=IS%r) then
          isBIBD=.False.
          !print *, "Sum in row is non-ideal", is%sum_in_row(i), is%r
          return
       endif
    enddo

    ! If the incidence matrix has a row with sum non-equal to k, this is not a BIBD
    do i=1,IS%b
       if(IS%SUM_IN_COL(i)/=IS%k) then
          isBIBD=.False.
          !print *, "Sum in col is non-ideal", is%sum_in_col(i), is%k
          return
       endif
    enddo

    ! If there are two IS%b (rows) whose intersection does not contain lambda IS%v,
    ! this is not a BIBD
    do i=1,(IS%v-1)
       do j=(i+1),IS%v
          if (IS%ROW_INTERSECTION(i,j)/=IS%LAMBDA) then
             isBIBD=.False.
             !print *, "Intersection is non-lambda", is%row_intersection(i,j)
             return
          endif
       enddo
    end do

    ! Otherwise, we got a BIBD
    isBIBD=.True.
    print *, "Generations: ", is%generations
    return
  end function isBIBD

  subroutine writeMatrix(IS)
    type(IncidenceStructure) IS
    integer row,col

    write(*,*),"["
    do row=1,IS%v
       do col=1,IS%b
          write (*,"(I1)",advance='no') IS%incidences(row,col)
       enddo
       write (*,"(A1)",advance='no') "-"
       write (*,"(I1)",advance='no') IS%SUM_IN_ROW(row)
       print *
    enddo
    do col=1,IS%b
       write (*,"(A1)",advance='no') "|"
    enddo
    print *
    do col=1,IS%b
       write (*,"(I1)",advance='no') IS%SUM_IN_col(col)
    enddo
    print *
    write(*,*),"]"
  end subroutine writeMatrix

end module incidence_structure
